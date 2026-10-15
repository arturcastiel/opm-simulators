// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \brief Test a parallel run whose LGRs keep no active cell.
 *
 * MINPV removes every host cell of the only LGR, so the trilinear refinement
 * refines nothing.  The grid must still be set up on every process.
 */

#define BOOST_TEST_MODULE LgrUnrefinedRanks

#include "SimulatorFixture.hpp"

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <fstream>
#include <string>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

// A 4x4x1 grid.  LGR1 refines cells (1..2,1..2,1), whose small porosity puts
// them below MINPV.
const char* deck = R"(RUNSPEC
DIMENS
 4 4 1 /
OIL
WATER
METRIC
START
 1 'JAN' 2020 /
EQLDIMS
/
TABDIMS
/
GRID
CARFIN
 'LGR1' 1 2 1 2 1 1 6 6 1 /
ENDFIN
DXV
 4*100.0 /
DYV
 4*100.0 /
DZ
 16*10.0 /
TOPS
 16*2000 /
PORO
 16*0.25 /
PERMX
 16*500 /
COPY
 PERMX PERMY /
 PERMX PERMZ /
/
BOX
 1 2 1 2 1 1 /
PORO
 4*0.001 /
ENDBOX
MINPV
 1000.0 /
PROPS
PVTW
 250.0 1.02 4.0E-05 0.5 0.0 /
ROCK
 250.0 4.0E-05 /
DENSITY
 850 1014 1.0 /
PVDO
 100.0 1.05 2.0
 250.0 1.02 2.5
 400.0 1.00 3.0 /
SWOF
 0.2 0.0  1.0 0.0
 0.8 0.5  0.0 0.0 /
SOLUTION
EQUIL
 2000 250.0 2100 0.0 /
SCHEDULE
TSTEP
 10 /
END
)";

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(LgrsWithoutActiveCells)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    const auto& comm = Opm::FlowGenericVanguard::comm();
    const auto dir = std::filesystem::temp_directory_path() / ("test_lgr_unrefined_ranks_np" + std::to_string(comm.size()));
    const auto deck_file = (dir / "LGR_UNREFINED_RANKS.DATA").string();
    if (comm.rank() == 0) {
        std::filesystem::create_directories(dir);
        std::ofstream os(deck_file);
        os << deck;
    }
    comm.barrier();

    // The trilinear refinement is flow's default.
    std::unique_ptr<Opm::GetPropType<TypeTag, Opm::Properties::Simulator>> simulator;
    BOOST_CHECK_NO_THROW(simulator = Opm::initSimulator<TypeTag>(deck_file.c_str(),
                                                                 "test_lgr_unrefined_ranks",
                                                                 /*threads_per_process=*/1,
                                                                 "low"));
    BOOST_REQUIRE(simulator);

    // Nothing is refined, and the 12 cells MINPV keeps are all simulated.
    const auto& grid = simulator->vanguard().grid();
    BOOST_CHECK_EQUAL(comm.max(grid.maxLevel()), 0);
    int interior = 0;
    for (const auto& elem : elements(grid.leafGridView(), Dune::Partitions::interior)) {
        static_cast<void>(elem);
        ++interior;
    }
    BOOST_CHECK_EQUAL(comm.sum(interior), 12);

    comm.barrier();
    if (comm.rank() == 0) {
        std::filesystem::remove_all(dir);
    }
}
