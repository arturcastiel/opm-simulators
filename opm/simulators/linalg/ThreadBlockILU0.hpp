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

/*! \brief Block-Jacobi ILU(0) over contiguous row ranges, one OpenMP task per range.
 *
 *  The rows of the matrix are split into a fixed number of contiguous ranges. Each range
 *  gets the ILU(0) factorisation of its own diagonal block of the matrix; couplings
 *  between ranges are dropped. Factorisation and the triangular solves of the ranges are
 *  independent, so update() and apply() run in parallel over the ranges.
 *
 *  The number of ranges does not depend on the number of threads, so the preconditioner,
 *  and with it the linear iterations, are the same at any thread count. With one range
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

    /*! \param A       the matrix (kept by reference; its values may change between updates)
     *  \param ranges  number of contiguous row ranges (at least 1, at most the number of rows)
     *  \param w       relaxation factor applied to the result
     */
    ThreadBlockILU0(const M& A, int ranges, field_type w)
        : A_(A)
        , w_(w)
    {
        const std::size_t n = A_.N();
        const std::size_t nr = std::max<std::size_t>(1, std::min<std::size_t>(ranges < 1 ? 1 : ranges, n));
        rangeStart_.resize(nr + 1);
        for (std::size_t r = 0; r <= nr; ++r) {
            rangeStart_[r] = (n * r) / nr;
        }

        // Pattern of each row restricted to its own range, columns in the order of A.
        rowStart_.assign(n + 1, 0);
        diag_.assign(n, 0);
        for (std::size_t r = 0; r < nr; ++r) {
            const std::size_t lo = rangeStart_[r];
            const std::size_t hi = rangeStart_[r + 1];
            for (std::size_t i = lo; i < hi; ++i) {
                std::size_t count = 0;
                bool hasDiag = false;
                for (auto it = A_[i].begin(); it != A_[i].end(); ++it) {
                    const std::size_t j = it.index();
                    if (j >= lo && j < hi) {
                        ++count;
                        hasDiag = hasDiag || (j == i);
                    }
                }
                if (!hasDiag) {
                    DUNE_THROW(ISTLError, "ThreadBlockILU0: diagonal entry missing");
                }
                rowStart_[i + 1] = count;
            }
        }
        for (std::size_t i = 0; i < n; ++i) {
            rowStart_[i + 1] += rowStart_[i];
        }
        cols_.resize(rowStart_[n]);
        vals_.resize(rowStart_[n]);
        for (std::size_t r = 0; r < nr; ++r) {
            const std::size_t lo = rangeStart_[r];
            const std::size_t hi = rangeStart_[r + 1];
            for (std::size_t i = lo; i < hi; ++i) {
                std::size_t k = rowStart_[i];
                for (auto it = A_[i].begin(); it != A_[i].end(); ++it) {
                    const std::size_t j = it.index();
                    if (j >= lo && j < hi) {
                        if (j == i) {
                            diag_[i] = k;
                        }
                        cols_[k++] = j;
                    }
                }
            }
        }
        update();
    }

    void pre(X&, Y&) override {}
    void post(X&) override {}

    //! \brief Copies the current matrix values and factorises every range.
    void update() override
    {
        const int nr = static_cast<int>(rangeStart_.size()) - 1;
        bool ok = true;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) reduction(&& : ok)
#endif
        for (int r = 0; r < nr; ++r) {
            ok = factorRange(rangeStart_[r], rangeStart_[r + 1]) && ok;
        }
        if (!ok) {
            DUNE_THROW(ISTLError, "ThreadBlockILU0: failed to invert a diagonal block");
        }
    }

    //! \brief Forward and backward solve of every range.
    void apply(X& v, const Y& d) override
    {
        const int nr = static_cast<int>(rangeStart_.size()) - 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1)
#endif
        for (int r = 0; r < nr; ++r) {
            solveRange(rangeStart_[r], rangeStart_[r + 1], v, d);
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
    bool factorRange(const std::size_t lo, const std::size_t hi)
    {
        // copy the values of the range's diagonal block of A
        for (std::size_t i = lo; i < hi; ++i) {
            std::size_t k = rowStart_[i];
            for (auto it = A_[i].begin(); it != A_[i].end(); ++it) {
                const std::size_t j = it.index();
                if (j >= lo && j < hi) {
                    vals_[k++] = *it;
                }
            }
        }
        // left-looking ILU(0) with stored inverse diagonal (as bilu0_decomposition)
        for (std::size_t i = lo; i < hi; ++i) {
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

    void solveRange(const std::size_t lo, const std::size_t hi, X& v, const Y& d) const
    {
        using dblock = typename Y::block_type;
        using vblock = typename X::block_type;
        // lower solve, unit diagonal
        for (std::size_t i = lo; i < hi; ++i) {
            dblock rhs(d[i]);
            for (std::size_t k = rowStart_[i]; k < diag_[i]; ++k) {
                vals_[k].mmv(v[cols_[k]], rhs);
            }
            v[i] = rhs;
        }
        // upper solve with the stored inverse diagonal
        for (std::size_t i = hi; i-- > lo;) {
            vblock rhs(v[i]);
            for (std::size_t k = diag_[i] + 1; k < rowStart_[i + 1]; ++k) {
                vals_[k].mmv(v[cols_[k]], rhs);
            }
            vals_[diag_[i]].mv(rhs, v[i]);
        }
        if (w_ != 1.0) {
            for (std::size_t i = lo; i < hi; ++i) {
                v[i] *= w_;
            }
        }
    }

    const M& A_;
    field_type w_;
    std::vector<std::size_t> rangeStart_;
    std::vector<std::size_t> rowStart_;
    std::vector<std::size_t> diag_;
    std::vector<std::size_t> cols_;
    std::vector<block_type> vals_;
};

} // namespace Dune

#endif // OPM_THREAD_BLOCK_ILU0_HEADER_INCLUDED
