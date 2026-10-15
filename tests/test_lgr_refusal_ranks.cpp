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
 * \brief Test that every process stops when the trilinear refinement refuses
 *        a box.
 *
 * The trilinear refinement cannot refine some cells, such as a cell on a fault
 * or a cell with a non-neighbour connection.  Each process refines its own
 * cells, and only the process that holds such a cell would see it: every
 * process must stop with the refusal.
 */

#define BOOST_TEST_MODULE LgrRefusalRanks

#include "SimulatorFixture.hpp"

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <fstream>
#include <string>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

// An 8x1x2 grid; LGR1 refines cell (1,1,1).  The GRID section ends with the
// lines the test case puts in.
std::string deck(const std::string& tops, const std::string& gridEnd)
{
    return R"(RUNSPEC
DIMENS
 8 1 2 /
OIL
WATER
METRIC
START
 1 'JAN' 2020 /
TABDIMS
/
EQLDIMS
/
GRID
CARFIN
'LGR1' 1 1 1 1 1 1 2 2 2 /
ENDFIN
DX
 16*100 /
DY
 16*100 /
DZ
 16*10 /
TOPS
)" + tops + R"( /
PORO
 16*0.2 /
PERMX
 16*100 /
PERMY
 16*100 /
PERMZ
 16*10 /
)" + gridEnd + R"(
PROPS
SWOF
 0.2 0.0 1.0 0
 1.0 1.0 0.0 0 /
PVTW
 200 1.0 4e-5 0.5 0 /
PVDO
 100 1.05 1.0
 300 1.03 1.0 /
DENSITY
 800 1000 1 /
ROCK
 200 1e-5 /
SOLUTION
EQUIL
 2000 200 2100 0 1000 0 /
SCHEDULE
TSTEP
 1 /
END
)";
}

// Every process stops: the deck cannot be refined by the trilinear refinement,
// which is flow's default.
void checkEveryProcessStops(const std::string& deckText)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    const auto& comm = Opm::FlowGenericVanguard::comm();
    const auto dir = std::filesystem::temp_directory_path() / ("test_lgr_refusal_ranks_np" + std::to_string(comm.size()));
    const auto deck_file = (dir / "LGR_REFUSAL_RANKS.DATA").string();
    if (comm.rank() == 0) {
        std::filesystem::create_directories(dir);
        std::ofstream os(deck_file);
        os << deckText;
    }
    comm.barrier();

    BOOST_CHECK_THROW(Opm::initSimulator<TypeTag>(deck_file.c_str(), "test_lgr_refusal_ranks",
                                                  /*threads_per_process=*/1, "low"),
                      std::exception);

    comm.barrier();
    if (comm.rank() == 0) {
        std::filesystem::remove_all(dir);
    }
}

} // Anonymous namespace

// The columns from I=2 on lie 5 m deeper, so the first column meets them across
// a fault and its cells have more than six faces.
BOOST_AUTO_TEST_CASE(BoxOnAFaultedCell)
{
    checkEveryProcessStops(deck(" 2000 7*2005", ""));
}

// An NNC joins the refined cell to the far end of the grid.
BOOST_AUTO_TEST_CASE(BoxOnACellWithAnNnc)
{
    checkEveryProcessStops(deck(" 8*2000", "NNC\n 1 1 1 8 1 2 10.0 /\n/"));
}
