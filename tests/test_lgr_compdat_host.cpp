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
 * \brief Test that a COMPDAT connection in a refined cell stops the run
 *
 * LGR_COMPDAT_HOST.DATA refines cell (1,1,1) of a 3x1x1 grid into 2x1x1 and
 * completes well PROD in that cell with COMPDAT.  The refined cell is not part
 * of the simulation grid, so the simulator must stop and name WELSPECL and
 * COMPDATL, with which a well is completed in an LGR.
 */

#define BOOST_TEST_MODULE LgrCompdatHost

#include "SimulatorFixture.hpp"

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <string>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

BOOST_AUTO_TEST_CASE(CompdatInRefinedCellStops)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    const auto namesCompdatl = [](const std::invalid_argument& e)
    {
        const auto msg = std::string{e.what()};
        return (msg.find("COMPDATL") != std::string::npos)
            && (msg.find("well PROD in cell (1,1,1)") != std::string::npos);
    };

    BOOST_CHECK_EXCEPTION(Opm::initSimulator<TypeTag>("LGR_COMPDAT_HOST.DATA",
                                                      "test_lgr_compdat_host",
                                                      /*threads_per_process=*/1, "low"),
                          std::invalid_argument, namesCompdatl);
}
