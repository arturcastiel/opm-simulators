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
 * \brief Test the transmissibilities of a grid with an LGR and an NNC on many
 *        processes.
 *
 * With one cell per process, the NNC joins cells of two other processes than
 * the one that refines the LGR.  Every process must set up its
 * transmissibilities, and the NNC must keep its transmissibility.
 */

#define BOOST_TEST_MODULE LgrNncRanks

#include "SimulatorFixture.hpp"

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <fstream>
#include <string>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

// A 4x1x1 grid.  LGR1 refines cell (1,1,1) into two cells, and an NNC joins
// cells (2,1,1) and (4,1,1).
const char* deck = R"(RUNSPEC
DIMENS
  4 1 1 /
EQLDIMS
/
TABDIMS
/
OIL
WATER
FIELD
START
  1 'JAN' 2015 /
GRID
CARFIN
'LGR1'  1  1  1  1  1  1  2  1  1 /
ENDFIN
DX
  4*1000 /
DY
  4*1000 /
DZ
  4*50 /
TOPS
  4*8325 /
PORO
  4*0.3 /
PERMX
  4*500 /
PERMY
  4*500 /
PERMZ
  4*50 /
NNC
  2 1 1   4 1 1   10.0 /
/
PROPS
PVTW
  4017.55 1.038 3.22E-6 0.318 0.0 /
PVDO
  14.7  1.0  1.0
  9014.7  0.9  1.0 /
ROCK
  14.7 3E-6 /
SWOF
  0.12  0        1      0
  1     0.00001  0      0 /
DENSITY
  53.66 64.49 0.0533 /
SOLUTION
EQUIL
  8400 4800 8450 0 /
SCHEDULE
TSTEP
  1 /
END
)";

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(TransmissibilitiesOfAnLgrAndAnNnc)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    const auto& comm = Opm::FlowGenericVanguard::comm();
    const auto dir = std::filesystem::temp_directory_path() / ("test_lgr_nnc_ranks_np" + std::to_string(comm.size()));
    const auto deck_file = (dir / "LGR_NNC_RANKS.DATA").string();
    if (comm.rank() == 0) {
        std::filesystem::create_directories(dir);
        std::ofstream os(deck_file);
        os << deck;
    }
    comm.barrier();

    // The trilinear refinement.
    std::unique_ptr<Opm::GetPropType<TypeTag, Opm::Properties::Simulator>> simulator;
    BOOST_CHECK_NO_THROW(simulator = Opm::initSimulator<TypeTag>(deck_file.c_str(),
                                                                 "test_lgr_nnc_ranks",
                                                                 /*threads_per_process=*/1,
                                                                 "low", { "--lgr-backend=trilinear" }));
    BOOST_REQUIRE(simulator);

    // The process that holds both cells of the NNC as interior or overlap cells
    // has its transmissibility: 10 in the deck's units.
    const auto& vanguard = simulator->vanguard();
    const auto& trans = simulator->problem().eclTransmissibilities();
    const auto& cartMapper = vanguard.cartesianIndexMapper();
    int cell2 = -1;
    int cell4 = -1;
    for (const auto& elem : elements(vanguard.gridView())) {
        const int idx = vanguard.gridView().indexSet().index(elem);
        if (elem.hasFather()) {
            continue;
        }
        const int cart = cartMapper.cartesianIndex(idx);
        if (cart == 1) {
            cell2 = idx;
        }
        if (cart == 3) {
            cell4 = idx;
        }
    }
    const double nncTrans = (cell2 >= 0 && cell4 >= 0)
        ? trans.transmissibility(cell2, cell4) : 0.0;
    const double deckTrans = Opm::UnitSystem::newFIELD()
        .to_si(Opm::UnitSystem::measure::transmissibility, 10.0);
    BOOST_CHECK_CLOSE(comm.max(nncTrans), deckTrans, 1.0e-8);

    comm.barrier();
    if (comm.rank() == 0) {
        std::filesystem::remove_all(dir);
    }
}
