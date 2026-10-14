/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS

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

#ifndef AQUIFER_GRID_UTILS_HEADER_INCLUDED
#define AQUIFER_GRID_UTILS_HEADER_INCLUDED

#include <opm/grid/CpGrid.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/Aquancon.hpp>
#include <opm/input/eclipse/EclipseState/Grid/LgrCollection.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <vector>

namespace Opm {

template<class Grid>
struct IsNumericalAquiferCell {
    explicit IsNumericalAquiferCell(const Grid&)
    {}

    template<class T>
    bool operator()(const T&) const { return false; }
};

template<>
struct IsNumericalAquiferCell<Dune::CpGrid> {
    explicit IsNumericalAquiferCell(const Dune::CpGrid& grid)
        : grid_(grid)
    {}

    template<class T>
    bool operator()(const T& elem) const
    {
        const auto& aquiferCells = grid_.sortedNumAquiferCells();
        if (aquiferCells.empty())
        {
          return false;
        }
        auto candidate = std::lower_bound(aquiferCells.begin(),
                                          aquiferCells.end(), elem.index());
        return candidate != aquiferCells.end() && *candidate == elem.index();
    }

private:
    const Dune::CpGrid& grid_;
};

/// \brief The connections of an analytical aquifer that reach a simulation cell.
///
/// AQUANCON names a cell by its level-zero index.  When a CARFIN box refines
/// that cell, the LGR cells take its place and the cell is not part of the
/// simulation grid: the connection is ignored, with a warning for each one.
/// Decided from the deck alone, so every rank keeps the same connections.
template<class Simulator>
std::vector<Aquancon::AquancCell>
aquiferConnectionsOutsideLgrs(const std::vector<Aquancon::AquancCell>& connections,
                              const Simulator& simulator,
                              const int aquiferID)
{
    const auto& vanguard = simulator.vanguard();
    const auto& lgrs = vanguard.eclState().getLgrs();
    if (lgrs.size() == 0) {
        return connections;
    }

    const auto& dims = vanguard.cartesianDimensions();
    const bool report = vanguard.grid().comm().rank() == 0;

    auto kept = std::vector<Aquancon::AquancCell>{};
    kept.reserve(connections.size());
    for (const auto& conn : connections) {
        const auto cell = static_cast<int>(conn.global_index);
        const auto ijk = std::array{cell % dims[0],
                                    (cell / dims[0]) % dims[1],
                                    cell / (dims[0] * dims[1])};
        bool refined = false;
        for (std::size_t l = 0; (l < lgrs.size()) && !refined; ++l) {
            // A nested box lies inside its parent's box.
            const auto& box = lgrs.getLgr(l);
            refined = (box.PARENT_NAME() == "GLOBAL")
                && (box.I1() <= ijk[0]) && (ijk[0] <= box.I2())
                && (box.J1() <= ijk[1]) && (ijk[1] <= box.J2())
                && (box.K1() <= ijk[2]) && (ijk[2] <= box.K2());
        }
        if (! refined) {
            kept.push_back(conn);
        }
        else if (report) {
            OpmLog::warning(fmt::format("Aquifer {}: the connection to cell ({},{},{}) is "
                                        "ignored. A CARFIN box refines that cell, so it is "
                                        "not part of the simulation grid.",
                                        aquiferID, ijk[0] + 1, ijk[1] + 1, ijk[2] + 1));
        }
    }

    return kept;
}

} // namespace Opm

#endif
