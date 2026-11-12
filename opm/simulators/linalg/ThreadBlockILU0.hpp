/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_THREAD_BLOCK_ILU0_HEADER_INCLUDED
#define OPM_THREAD_BLOCK_ILU0_HEADER_INCLUDED

#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/istlexception.hh>
#include <dune/istl/solvercategory.hh>

#include <algorithm>
#include <cstddef>
#include <vector>

#if HAVE_OPENMP
#include <omp.h>
#endif

namespace Dune
{

/*! \brief Block-Jacobi ILU(0) over a fixed partition of the rows, one OpenMP task per part.
 *
 *  The rows of the matrix are split into a fixed number of parts. Each part gets the
 *  ILU(0) factorisation of its own diagonal block of the matrix (its rows in their
 *  natural order); couplings between parts are dropped. Factorisation and the triangular
 *  solves of the parts are independent, so update() and apply() run in parallel over them.
 *
 *  Two partitions: "rows" cuts the row numbers into contiguous ranges; "columns" first
 *  finds the largest forward offset of a matrix entry (the layer size of a structured
 *  grid numbered layer by layer), and cuts the position of a row within that stride into
 *  contiguous ranges, so every part holds whole vertical columns of cells and only
 *  lateral couplings are dropped.
 *
 *  The number of parts does not depend on the number of threads, so the preconditioner,
 *  and with it the linear iterations, are the same at any thread count. With one part
 *  it is the plain ILU(0) of the whole matrix.
 *
 *  \tparam M The matrix type to operate on
 *  \tparam X Type of the update
 *  \tparam Y Type of the defect
 */
template <class M, class X, class Y>
class ThreadBlockILU0 : public PreconditionerWithUpdate<X, Y>
{
public:
    using matrix_type = M;
    using domain_type = X;
    using range_type = Y;
    using field_type = typename X::field_type;
    using block_type = typename M::block_type;

    /*! \param A          the matrix (kept by reference; its values may change between updates)
     *  \param parts      number of parts (at least 1, at most the number of rows)
     *  \param w          relaxation factor applied to the result
     *  \param byColumns  false: contiguous row ranges; true: ranges of the position within
     *                    the largest forward offset (whole vertical columns, see above)
     */
    ThreadBlockILU0(const M& A, int parts, field_type w, bool byColumns = false)
        : A_(A)
        , w_(w)
    {
        const std::size_t n = A_.N();
        const std::size_t np = std::max<std::size_t>(1, std::min<std::size_t>(parts < 1 ? 1 : parts, n));
        std::size_t stride = n;
        if (byColumns) {
            stride = 1;
            for (std::size_t i = 0; i < n; ++i) {
                for (auto it = A_[i].begin(); it != A_[i].end(); ++it) {
                    if (it.index() > i) {
                        stride = std::max(stride, it.index() - i);
                    }
                }
            }
        }
        partOf_.resize(n);
        partRows_.assign(np, {});
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t pos = byColumns ? i % stride : i;
            const std::size_t p = std::min(np - 1, (pos * np) / std::min(stride, n));
            partOf_[i] = static_cast<int>(p);
            partRows_[p].push_back(i);   // ascending: the natural order within the part
        }

        // Pattern of each row restricted to its own part, columns in the order of A.
        rowStart_.assign(n + 1, 0);
        diag_.assign(n, 0);
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t count = 0;
            bool hasDiag = false;
            for (auto it = A_[i].begin(); it != A_[i].end(); ++it) {
                const std::size_t j = it.index();
                if (partOf_[j] == partOf_[i]) {
                    ++count;
                    hasDiag = hasDiag || (j == i);
                }
            }
            if (!hasDiag) {
                DUNE_THROW(ISTLError, "ThreadBlockILU0: diagonal entry missing");
            }
            rowStart_[i + 1] = count;
        }
        for (std::size_t i = 0; i < n; ++i) {
            rowStart_[i + 1] += rowStart_[i];
        }
        cols_.resize(rowStart_[n]);
        vals_.resize(rowStart_[n]);
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t k = rowStart_[i];
            for (auto it = A_[i].begin(); it != A_[i].end(); ++it) {
                const std::size_t j = it.index();
                if (partOf_[j] == partOf_[i]) {
                    if (j == i) {
                        diag_[i] = k;
                    }
                    cols_[k++] = j;
                }
            }
        }
        update();
    }

    void pre(X&, Y&) override {}
    void post(X&) override {}

    //! \brief Copies the current matrix values and factorises every part.
    void update() override
    {
        const int np = static_cast<int>(partRows_.size());
        bool ok = true;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) reduction(&& : ok)
#endif
        for (int p = 0; p < np; ++p) {
            ok = factorPart(p) && ok;
        }
        if (!ok) {
            DUNE_THROW(ISTLError, "ThreadBlockILU0: failed to invert a diagonal block");
        }
    }

    //! \brief Forward and backward solve of every part.
    void apply(X& v, const Y& d) override
    {
        const int np = static_cast<int>(partRows_.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1)
#endif
        for (int p = 0; p < np; ++p) {
            solvePart(p, v, d);
        }
    }

    SolverCategory::Category category() const override
    {
        return SolverCategory::sequential;
    }

    bool hasPerfectUpdate() const override
    {
        return true;
    }

private:
    bool factorPart(const int p)
    {
        const auto& rows = partRows_[p];
        // copy the values of the part's diagonal block of A
        for (const std::size_t i : rows) {
            std::size_t k = rowStart_[i];
            for (auto it = A_[i].begin(); it != A_[i].end(); ++it) {
                if (partOf_[it.index()] == p) {
                    vals_[k++] = *it;
                }
            }
        }
        // left-looking ILU(0) with stored inverse diagonal (as bilu0_decomposition)
        for (const std::size_t i : rows) {
            const std::size_t rowEnd = rowStart_[i + 1];
            for (std::size_t ij = rowStart_[i]; ij < diag_[i]; ++ij) {
                const std::size_t j = cols_[ij];
                // L_ij = A_ij * A_jj^-1
                vals_[ij].rightmultiply(vals_[diag_[j]]);
                std::size_t ik = ij + 1;
                std::size_t jk = diag_[j] + 1;
                const std::size_t jEnd = rowStart_[j + 1];
                while (ik < rowEnd && jk < jEnd) {
                    if (cols_[ik] == cols_[jk]) {
                        block_type B(vals_[jk]);
                        B.leftmultiply(vals_[ij]);
                        vals_[ik] -= B;
                        ++ik;
                        ++jk;
                    } else if (cols_[ik] < cols_[jk]) {
                        ++ik;
                    } else {
                        ++jk;
                    }
                }
            }
            try {
                vals_[diag_[i]].invert();
            }
            catch (const FMatrixError&) {
                return false;
            }
        }
        return true;
    }

    void solvePart(const int p, X& v, const Y& d) const
    {
        using dblock = typename Y::block_type;
        using vblock = typename X::block_type;
        const auto& rows = partRows_[p];
        // lower solve, unit diagonal
        for (const std::size_t i : rows) {
            dblock rhs(d[i]);
            for (std::size_t k = rowStart_[i]; k < diag_[i]; ++k) {
                vals_[k].mmv(v[cols_[k]], rhs);
            }
            v[i] = rhs;
        }
        // upper solve with the stored inverse diagonal
        for (auto r = rows.rbegin(); r != rows.rend(); ++r) {
            const std::size_t i = *r;
            vblock rhs(v[i]);
            for (std::size_t k = diag_[i] + 1; k < rowStart_[i + 1]; ++k) {
                vals_[k].mmv(v[cols_[k]], rhs);
            }
            vals_[diag_[i]].mv(rhs, v[i]);
        }
        if (w_ != 1.0) {
            for (const std::size_t i : rows) {
                v[i] *= w_;
            }
        }
    }

    const M& A_;
    field_type w_;
    std::vector<int> partOf_;
    std::vector<std::vector<std::size_t>> partRows_;
    std::vector<std::size_t> rowStart_;
    std::vector<std::size_t> diag_;
    std::vector<std::size_t> cols_;
    std::vector<block_type> vals_;
};

} // namespace Dune

#endif // OPM_THREAD_BLOCK_ILU0_HEADER_INCLUDED
