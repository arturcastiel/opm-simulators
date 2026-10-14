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
 * \brief Test that an aquifer connection on a refined cell is ignored
 *
 * LGR_AQUIFER_HOST.DATA refines cell (1,1,1) of a 3x1x1 grid into 2x1x1 and
 * connects aquifer 1 to (1,1,1) and to (3,1,1).  The refined cell is not part
 * of the simulation grid, so the aquifer keeps the connection to (3,1,1) only.
 */

#define BOOST_TEST_MODULE LgrAquiferHost

#include "SimulatorFixture.hpp"

#include <opm/input/eclipse/EclipseState/Aquifer/AquiferConfig.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/simulators/aquifers/AquiferGridUtils.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

BOOST_AUTO_TEST_CASE(ConnectionOnRefinedCellIsIgnored)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    auto simulator = Opm::initSimulator<TypeTag>("LGR_AQUIFER_HOST.DATA", "test_lgr_aquifer_host",
                                                 /*threads_per_process=*/1, "low");
    const auto& vanguard = simulator->vanguard();
    BOOST_REQUIRE_EQUAL(vanguard.grid().maxLevel(), 1);

    const auto& deckConnections = vanguard.eclState().aquifer().connections().getConnections(1);
    BOOST_REQUIRE_EQUAL(deckConnections.size(), 2U);

    const auto kept = Opm::aquiferConnectionsOutsideLgrs(deckConnections, *simulator, 1);
    BOOST_REQUIRE_EQUAL(kept.size(), 1U);
    BOOST_CHECK_EQUAL(kept.front().global_index, 2U);

    // The aquifer tracer finds its cells through the same map: no cell for
    // the refined cell either.
    BOOST_CHECK_EQUAL(vanguard.compressedIndex(0), -1);
}
