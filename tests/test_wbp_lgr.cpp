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
*/

#define BOOST_TEST_MODULE WellBlockAveragePressuresLGR

#include "SimulatorFixture.hpp"

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/output/data/Wells.hpp>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/utils/readDeck.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

namespace Opm::Properties::TTag {
    struct TestWbpLgrTypeTag {
        using InheritsFrom = std::tuple<TestTypeTag>;
    };
}

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

// A 6x6x3 grid of 30 m x 30 m x 5 m cells, given either directly (first
// argument empty) or as an LGR that refines a 2x2x1 grid of 90 m x 90 m x
// 15 m cells into the same 6x6x3 cells.  The well is completed in the
// X direction in cell (1,1,1), a corner, so that its block average
// pressures depend on which of the neighbouring cells exist.
std::string deck(const bool lgr)
{
    const std::string grid = lgr
        ? R"(DIMENS
 2 2 1 /
LGR
 1 100 /
)" : R"(DIMENS
 6 6 3 /
)";

    const std::string geometry = lgr
        ? R"(CARFIN
 'LGR1' 1 2 1 2 1 1 6 6 3 /
ENDFIN
DXV
 2*90.0 /
DYV
 2*90.0 /
DZ
 4*15.0 /
TOPS
 4*2000 /
PORO
 4*0.25 /
PERMX
 4*500 /
)" : R"(DXV
 6*30.0 /
DYV
 6*30.0 /
DZ
 108*5.0 /
TOPS
 36*2000 /
PORO
 108*0.25 /
PERMX
 108*500 /
)";

    const std::string well = lgr
        ? R"(WELSPECL
 PROD G LGR1 1 1 1* OIL /
/
COMPDATL
 PROD LGR1 1 1 1 1 OPEN 1* 1* 0.2 1* 1* 1* X /
/
)" : R"(WELSPECS
 PROD G 1 1 1* OIL /
/
COMPDAT
 PROD 1 1 1 1 OPEN 1* 1* 0.2 1* 1* 1* X /
/
)";

    return "RUNSPEC\n" + grid + R"(OIL
WATER
METRIC
START
 1 'JAN' 2020 /
WELLDIMS
 2 10 2 2 /
EQLDIMS
/
TABDIMS
/
GRID
)" + geometry + R"(COPY
 PERMX PERMY /
 PERMX PERMZ /
/
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
 0.3 0.05 0.6 0.0
 0.5 0.15 0.3 0.0
 0.7 0.35 0.1 0.0
 0.8 0.5  0.0 0.0 /
SOLUTION
EQUIL
 2000 250.0 2100 0.0 /
SCHEDULE
)" + well + R"(WCONPROD
 PROD OPEN BHP 5* 150 /
/
TSTEP
 1 /
END
)";
}

// A simulator for the deck.  The deck is read with low parsing strictness,
// since the LGR keywords are flagged in normal mode.
template <class TypeTag>
std::unique_ptr<Opm::GetPropType<TypeTag, Opm::Properties::Simulator>>
initSimulator(const std::string& filename)
{
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

    std::vector<std::string> args { "test_wbp_lgr", "--ecl-deck-file-name=" + filename };
    std::vector<const char*> argv;
    for (const auto& arg : args) {
        argv.push_back(arg.c_str());
    }

    Opm::Parameters::reset();
    Opm::registerAllParameters_<TypeTag>(false);
    Opm::registerEclTimeSteppingParameters<double>();
    Opm::BlackoilModelParameters<double>::registerParameters();
    Opm::Parameters::Register<Opm::Parameters::EnableTerminalOutput>("Do *NOT* use!");
    Opm::Parameters::SetDefault<Opm::Parameters::ThreadsPerProcess>(1);
    Opm::Parameters::endRegistration();
    Opm::setupParameters_<TypeTag>(static_cast<int>(argv.size()), argv.data(),
                                   /*registerParams=*/false,
                                   /*allowUnused=*/false,
                                   /*handleHelp=*/true,
                                   /*myRank=*/0);

    auto& params = Opm::FlowGenericVanguard::modelParams_;
    Opm::readDeck(Opm::FlowGenericVanguard::comm(), filename,
                  params.eclState_, params.eclSchedule_, params.udqState_,
                  params.actionState_, params.wtestState_, params.eclSummaryConfig_,
                  nullptr, "low", "normal", "100", false, false, false, {},
                  /*slaveMode=*/false);

    return std::make_unique<Simulator>();
}

// Block average pressures of the well at the initial state.
Opm::data::WellBlockAvgPress initialWbp(const bool lgr)
{
    using TypeTag = Opm::Properties::TTag::TestWbpLgrTypeTag;

    const auto dir = std::filesystem::temp_directory_path() /
        (lgr ? "test_wbp_lgr_lgr" : "test_wbp_lgr_plain");
    std::filesystem::create_directories(dir);
    const auto deck_file = (dir / "WBP.DATA").string();
    {
        std::ofstream os(deck_file);
        os << deck(lgr);
    }

    auto simulator = initSimulator<TypeTag>(deck_file);

    simulator->model().applyInitialSolution();
    simulator->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(Opm::unit::day);
    simulator->problem().resetIterationForNewTimestep();

    auto& wellModel = simulator->problem().wellModel();
    wellModel.beginReportStep(0);

    const auto wbp = wellModel.wellBlockAveragePressures();
    BOOST_REQUIRE_EQUAL(wbp.values.count("PROD"), 1u);

    std::filesystem::remove_all(dir);
    return wbp.values.at("PROD");
}

} // Anonymous namespace

// A well in an LGR has the block average pressures of the same well in the
// same cells without an LGR.
BOOST_AUTO_TEST_CASE(WellInLgrAsWithoutLgr)
{
    using Quantity = Opm::data::WellBlockAvgPress::Quantity;

    const auto plain = initialWbp(false);
    const auto lgr = initialWbp(true);

    for (const auto q : { Quantity::WBP, Quantity::WBP4, Quantity::WBP5, Quantity::WBP9 }) {
        BOOST_CHECK_GT(plain[q], 0.0);
        BOOST_CHECK_CLOSE(lgr[q], plain[q], 1.0e-8);
    }
}
