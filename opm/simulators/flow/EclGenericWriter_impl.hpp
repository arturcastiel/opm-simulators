// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
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
#ifndef OPM_ECL_GENERIC_WRITER_IMPL_HPP
#define OPM_ECL_GENERIC_WRITER_IMPL_HPP

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/grid/GridHelpers.hpp>
#include <opm/grid/utility/cartesianToCompressed.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/RegionSetMatcher.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/RPTConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMatcher.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/output/eclipse/Summary.hpp>

#include <opm/simulators/flow/EclGenericWriter.hpp>

#if HAVE_MPI
#include <opm/simulators/utils/MPISerializer.hpp>
#endif

#if HAVE_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

/*!
 * \brief Detect whether two cells are direct vertical neighbours.
 *
 * I.e. have the same i and j index and all cartesian cells between them
 * along the vertical column are inactive.
 *
 * \tparam CM The type of the cartesian index mapper.
 * \param cartMapper The mapper onto cartesian indices.
 * \param cartesianToActive The mapping of cartesian indices to active indices.
 * \param smallGlobalIndex The cartesian cell index of the cell with smaller index
 * \param largeGlobalIndex The cartesian cell index of the cell with larger index
 * \return True if the cells have the same i and j indices and all cartesian cells
 *         between them are inactive.
 */
bool directVerticalNeighbors(const std::array<int, 3>& cartDims,
                             const std::unordered_map<int,int>& cartesianToActive,
                             int smallGlobalIndex, int largeGlobalIndex)
{
    assert(smallGlobalIndex <= largeGlobalIndex);
    std::array<int, 3> ijk1, ijk2;
    auto globalToIjk = [cartDims](int gc) {
                           std::array<int, 3> ijk;
                           ijk[0] = gc % cartDims[0];
                           gc /= cartDims[0];
                           ijk[1] = gc % cartDims[1];
                           ijk[2] = gc / cartDims[1];
                           return ijk;
                       };
    ijk1 = globalToIjk(smallGlobalIndex);
    ijk2 = globalToIjk(largeGlobalIndex);
    assert(ijk2[2]>=ijk1[2]);

    if ( ijk1[0] == ijk2[0] && ijk1[1] == ijk2[1] && (ijk2[2] - ijk1[2]) > 1)
    {
        assert((largeGlobalIndex-smallGlobalIndex)%(cartDims[0]*cartDims[1])==0);
        for ( int gi = smallGlobalIndex + cartDims[0] * cartDims[1]; gi < largeGlobalIndex;
              gi += cartDims[0] * cartDims[1] )
        {
            if ( cartesianToActive.find( gi ) != cartesianToActive.end() )
            {
                return false;
            }
        }
        return true;
    } else
        return false;
}

std::unordered_map<std::string, Opm::data::InterRegFlowMap>
getInterRegFlowsAsMap(const Opm::InterRegFlowMap& map)
{
    auto maps = std::unordered_map<std::string, Opm::data::InterRegFlowMap>{};

    const auto& regionNames = map.names();
    auto flows = map.getInterRegFlows();
    const auto nmap = regionNames.size();

    maps.reserve(nmap);
    for (auto mapID = 0*nmap; mapID < nmap; ++mapID) {
        maps.emplace(regionNames[mapID], std::move(flows[mapID]));
    }

    return maps;
}

struct EclWriteTasklet : public Opm::TaskletInterface
{
    Opm::Action::State actionState_;
    Opm::WellTestState wtestState_;
    Opm::SummaryState summaryState_;
    Opm::UDQState udqState_;
    Opm::EclipseIO& eclIO_;
    int reportStepNum_;
    std::optional<int> timeStepNum_;
    bool isSubStep_;
    double secondsElapsed_;
    Opm::RestartValue restartValue_;
    bool writeDoublePrecision_;

    explicit EclWriteTasklet(const Opm::Action::State& actionState,
                             const Opm::WellTestState& wtestState,
                             const Opm::SummaryState& summaryState,
                             const Opm::UDQState& udqState,
                             Opm::EclipseIO& eclIO,
                             int reportStepNum,
                             std::optional<int> timeStepNum,
                             bool isSubStep,
                             double secondsElapsed,
                             Opm::RestartValue restartValue,
                             bool writeDoublePrecision)
        : actionState_(actionState)
        , wtestState_(wtestState)
        , summaryState_(summaryState)
        , udqState_(udqState)
        , eclIO_(eclIO)
        , reportStepNum_(reportStepNum)
        , timeStepNum_(timeStepNum)
        , isSubStep_(isSubStep)
        , secondsElapsed_(secondsElapsed)
        , restartValue_(std::move(restartValue))
        , writeDoublePrecision_(writeDoublePrecision)
    {}

    // callback to eclIO serial writeTimeStep method
    void run() override
    {
        this->eclIO_.writeTimeStep(this->actionState_,
                                   this->wtestState_,
                                   this->summaryState_,
                                   this->udqState_,
                                   this->reportStepNum_,
                                   this->isSubStep_,
                                   this->secondsElapsed_,
                                   std::move(this->restartValue_),
                                   this->writeDoublePrecision_,
                                   this->timeStepNum_
);
    }
};

}

namespace Opm {

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
EclGenericWriter(const Schedule& schedule,
                 const EclipseState& eclState,
                 const SummaryConfig& summaryConfig,
                 const Grid& grid,
                 const EquilGrid* equilGrid,
                 const GridView& gridView,
                 const Dune::CartesianIndexMapper<Grid>& cartMapper,
                 const Dune::CartesianIndexMapper<EquilGrid>* equilCartMapper,
                 bool enableAsyncOutput,
                 bool enableEsmry )
    : collectOnIORank_(grid,
                       equilGrid,
                       gridView,
                       cartMapper,
                       equilCartMapper,
                       summaryConfig.fip_regions_interreg_flow())
    , grid_           (grid)
    , gridView_       (gridView)
    , schedule_       (schedule)
    , eclState_       (eclState)
    , cartMapper_     (cartMapper)
    , equilCartMapper_(equilCartMapper)
    , equilGrid_      (equilGrid)
{
    if (this->collectOnIORank_.isIORank()) {
        this->eclIO_ = std::make_unique<EclipseIO>
            (this->eclState_,
             UgGridHelpers::createEclipseGrid(*equilGrid, eclState_.getInputGrid()),
             this->schedule_, summaryConfig, "", enableEsmry);
    }

    // create output thread if enabled and rank is I/O rank
    // async output is enabled by default if pthread are enabled
    int numWorkerThreads = 0;
    if (enableAsyncOutput && collectOnIORank_.isIORank()) {
        numWorkerThreads = 1;
    }

    this->taskletRunner_.reset(new TaskletRunner(numWorkerThreads));
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
const EclipseIO& EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
eclIO() const
{
    assert(eclIO_);
    return *eclIO_;
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
writeInit()
{
    if (collectOnIORank_.isIORank()) {
        std::map<std::string, std::vector<int>> integerVectors;
        if (collectOnIORank_.isParallel()) {
            integerVectors.emplace("MPI_RANK", collectOnIORank_.globalRanks());
        }

        eclIO_->writeInitial(*this->outputTrans_,
                             integerVectors,
                             this->outputNnc_);
        this->outputTrans_.reset();
    }
}
template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
extractOutputTransAndNNC(const std::function<unsigned int(unsigned int)>& map)
{
    if (collectOnIORank_.isIORank()) {
        auto cartMap = cartesianToCompressed(equilGrid_->size(0), UgGridHelpers::globalCell(*equilGrid_));
        computeTrans_(cartMap, map);
        exportNncStructure_(cartMap, map);
    }

#if HAVE_MPI
    if (collectOnIORank_.isParallel()) {
        const auto& comm = grid_.comm();
        Parallel::MpiSerializer ser(comm);
        ser.broadcast(Parallel::RootRank{0}, outputNnc_);
    }
#endif
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
computeTrans_(const std::unordered_map<int,int>& cartesianToActive,
              const std::function<unsigned int(unsigned int)>& map) const
{
    if (!outputTrans_) {
        outputTrans_ = std::make_unique<data::Solution>();
    }

    const auto& cartMapper = *equilCartMapper_;
    const auto& cartDims = cartMapper.cartesianDimensions();

    auto createCellData = [&cartDims]() {
        return data::CellData{
                UnitSystem::measure::transmissibility,
                std::vector<double>(cartDims[0] * cartDims[1] * cartDims[2], 0.0),
                data::TargetType::INIT
        };
    };

    outputTrans_->clear();
    outputTrans_->emplace("TRANX", createCellData());
    outputTrans_->emplace("TRANY", createCellData());
    outputTrans_->emplace("TRANZ", createCellData());

    auto& tranx = this->outputTrans_->at("TRANX");
    auto& trany = this->outputTrans_->at("TRANY");
    auto& tranz = this->outputTrans_->at("TRANZ");

    using GlobalGridView = typename EquilGrid::LeafGridView;
    using GlobElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView>;
    const GlobalGridView& globalGridView = this->equilGrid_->leafGridView();
    const GlobElementMapper globalElemMapper { globalGridView, Dune::mcmgElementLayout() };

    auto isNumAquCell = [numAquCell = this->eclState_.aquifer().hasNumericalAquifer()
                         ? this->eclState_.aquifer().numericalAquifers().allAquiferCellIds()
                         : std::vector<std::size_t>{}]
        (const std::size_t cellIdx)
    {
        return std::binary_search(numAquCell.begin(), numAquCell.end(), cellIdx);
    };

    for (const auto& elem : elements(globalGridView)) {
        for (const auto& is : intersections(globalGridView, elem)) {
            if (!is.neighbor())
                continue; // intersection is on the domain boundary

            // Not 'const' because remapped if 'map' is non-null.
            unsigned c1 = globalElemMapper.index(is.inside());
            unsigned c2 = globalElemMapper.index(is.outside());

            if (c1 > c2)
                continue; // we only need to handle each connection once, thank you.

            // Ordering of compressed and uncompressed index should be the same
            const int cartIdx1 = cartMapper.cartesianIndex( c1 );
            const int cartIdx2 = cartMapper.cartesianIndex( c2 );

            if (isNumAquCell(cartIdx1) || isNumAquCell(cartIdx2)) {
                // Connections involving numerical aquifers are always NNCs
                // for the purpose of file output.  This holds even for
                // connections between cells like (I,J,K) and (I+1,J,K)
                // which are nominally neighbours in the Cartesian grid.
                continue;
            }

            // Ordering of compressed and uncompressed index should be the same
            assert(cartIdx1 <= cartIdx2);
            int gc1 = std::min(cartIdx1, cartIdx2);
            int gc2 = std::max(cartIdx1, cartIdx2);

            // Re-ordering in case of non-empty mapping between equilGrid to grid
            if (map) {
                c1 = map(c1); // equilGridToGrid map
                c2 = map(c2);
            }

            if (gc2 - gc1 == 1 && cartDims[0] > 1 ) {
                tranx.template data<double>()[gc1] = globalTrans().transmissibility(c1, c2);
                continue; // skip other if clauses as they are false, last one needs some computation
            }

            if (gc2 - gc1 == cartDims[0] && cartDims[1] > 1) {
                trany.template data<double>()[gc1] = globalTrans().transmissibility(c1, c2);
                continue; // skipt next if clause as it needs some computation
            }

            if ( gc2 - gc1 == cartDims[0]*cartDims[1] ||
                 directVerticalNeighbors(cartDims, cartesianToActive, gc1, gc2))
                tranz.template data<double>()[gc1] = globalTrans().transmissibility(c1, c2);
        }
    }
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
std::vector<NNCdata>
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
exportNncStructure_(const std::unordered_map<int,int>& cartesianToActive,
                    const std::function<unsigned int(unsigned int)>& map) const
{
    auto isNumAquCell = [numAquCell = this->eclState_.aquifer().hasNumericalAquifer()
                         ? this->eclState_.aquifer().numericalAquifers().allAquiferCellIds()
                         : std::vector<std::size_t>{}]
        (const std::size_t cellIdx)
    {
        return std::binary_search(numAquCell.begin(), numAquCell.end(), cellIdx);
    };

    auto isNumAquConn = [&isNumAquCell](const std::size_t cellIdx1,
                                        const std::size_t cellIdx2)
    {
        return isNumAquCell(cellIdx1) || isNumAquCell(cellIdx2);
    };

    auto isCartesianNeighbour = [nx = this->eclState_.getInputGrid().getNX(),
                                 ny = this->eclState_.getInputGrid().getNY()]
        (const std::size_t cellIdx1, const std::size_t cellIdx2)
    {
        const auto cellDiff = cellIdx2 - cellIdx1;

        return (cellDiff == 1)
            || (cellDiff == nx)
            || (cellDiff == nx * ny);
    };

    auto activeCell = [&cartesianToActive](const std::size_t cellIdx)
    {
        auto pos = cartesianToActive.find(cellIdx);
        return (pos == cartesianToActive.end()) ? -1 : pos->second;
    };

    const auto& nncData = this->eclState_.getInputNNC().input();
    const auto& unitSystem = this->eclState_.getDeckUnitSystem();

    for (const auto& entry : nncData) {
        // Ignore most explicit NNCs between otherwise neighbouring cells.
        // We keep NNCs that involve cells with numerical aquifers even if
        // these might be between neighbouring cells in the Cartesian
        // grid--e.g., between cells (I,J,K) and (I+1,J,K).  All such
        // connections should be written to NNC output arrays provided the
        // transmissibility value is sufficiently large.
        //
        // The condition cell2 >= cell1 holds by construction of nncData.
        assert (entry.cell2 >= entry.cell1);

        if (! isCartesianNeighbour(entry.cell1, entry.cell2) ||
            isNumAquConn(entry.cell1, entry.cell2))
        {
            // Pick up transmissibility value from 'globalTrans()' since
            // multiplier keywords like MULTREGT might have impacted the
            // values entered in primary sources like NNC/EDITNNC/EDITNNCR.
            const auto c1 = activeCell(entry.cell1);
            const auto c2 = activeCell(entry.cell2);

            if ((c1 < 0) || (c2 < 0)) {
                // Connection between inactive cells?  Unexpected at this
                // level.  Might consider 'throw'ing if this happens...
                continue;
            }

            const auto trans = this->globalTrans().transmissibility(c1, c2);
            const auto tt = unitSystem
                .from_si(UnitSystem::measure::transmissibility, trans);

            // ECLIPSE ignores NNCs (with EDITNNC/EDITNNCR applied) with
            // small transmissibility values.  Seems like the threshold is
            // 1.0e-6 in output units.
            if (std::isnormal(tt) && ! (tt < 1.0e-6)) {
                this->outputNnc_.emplace_back(entry.cell1, entry.cell2, trans);
            }
        }
    }

    auto isDirectNeighbours = [&isCartesianNeighbour, &cartesianToActive,
                               cartDims = &this->cartMapper_.cartesianDimensions()]
        (const std::size_t cellIdx1, const std::size_t cellIdx2)
    {
        return isCartesianNeighbour(cellIdx1, cellIdx2)
            || directVerticalNeighbors(*cartDims, cartesianToActive, cellIdx1, cellIdx2);
    };

    using GlobalGridView = typename EquilGrid::LeafGridView;
    using GlobElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView>;
    const GlobalGridView& globalGridView = this->equilGrid_->leafGridView();
    const GlobElementMapper globalElemMapper { globalGridView, Dune::mcmgElementLayout() };

    // Cartesian index mapper for the serial I/O grid
    const auto& equilCartMapper = *equilCartMapper_;
    for (const auto& elem : elements(globalGridView)) {
        for (const auto& is : intersections(globalGridView, elem)) {
            if (!is.neighbor())
                continue; // intersection is on the domain boundary

            // Not 'const' because remapped if 'map' is non-null.
            unsigned c1 = globalElemMapper.index(is.inside());
            unsigned c2 = globalElemMapper.index(is.outside());

            if (c1 > c2)
                continue; // we only need to handle each connection once, thank you.

            std::size_t cc1 = equilCartMapper.cartesianIndex( c1 );
            std::size_t cc2 = equilCartMapper.cartesianIndex( c2 );

            if ( cc2 < cc1 )
                std::swap(cc1, cc2);

            // Re-ordering in case of non-empty mapping between equilGrid to grid
            if (map) {
                c1 = map(c1); // equilGridToGrid map
                c2 = map(c2);
            }

            if (isNumAquConn(cc1, cc2) || ! isDirectNeighbours(cc1, cc2)) {
                // We need to check whether an NNC for this face was also
                // specified via the NNC keyword in the deck.
                auto t = this->globalTrans().transmissibility(c1, c2);
                auto candidate = std::lower_bound(nncData.begin(), nncData.end(),
                                                  NNCdata { cc1, cc2, 0.0 });

                while ((candidate != nncData.end()) &&
                       (candidate->cell1 == cc1) &&
                       (candidate->cell2 == cc2))
                {
                    t -= candidate->trans;
                    ++candidate;
                }

                // ECLIPSE ignores NNCs with zero transmissibility
                // (different threshold than for NNC with corresponding
                // EDITNNC above).  In addition we do set small
                // transmissibilities to zero when setting up the simulator.
                // These will be ignored here, too.
                const auto tt = unitSystem
                    .from_si(UnitSystem::measure::transmissibility, t);

                if (std::isnormal(tt) && (tt > 1.0e-12)) {
                    this->outputNnc_.emplace_back(cc1, cc2, t);
                }
            }
        }
    }

    return this->outputNnc_;
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
doWriteOutput(const int                          reportStepNum,
              const std::optional<int>           timeStepNum,
              const bool                         isSubStep,
              data::Solution&&                   localCellData,
              data::Wells&&                      localWellData,
              data::GroupAndNetworkValues&&      localGroupAndNetworkData,
              data::Aquifers&&                   localAquiferData,
              WellTestState&&                    localWTestState,
              const Action::State&               actionState,
              const UDQState&                    udqState,
              const SummaryState&                summaryState,
              const std::vector<Scalar>&         thresholdPressure,
              Scalar                             curTime,
              Scalar                             nextStepSize,
              bool                               doublePrecision,
              bool                               isFlowsn,
              std::array<FlowsData<double>, 3>&& flowsn,
              bool                               isFloresn,
              std::array<FlowsData<double>, 3>&& floresn)
{
    const auto isParallel = this->collectOnIORank_.isParallel();
    const bool needsReordering = this->collectOnIORank_.doesNeedReordering();

    RestartValue restartValue {
        (isParallel || needsReordering)
        ? this->collectOnIORank_.globalCellData()
        : std::move(localCellData),

        isParallel ? this->collectOnIORank_.globalWellData()
                   : std::move(localWellData),

        isParallel ? this->collectOnIORank_.globalGroupAndNetworkData()
                   : std::move(localGroupAndNetworkData),

        isParallel ? this->collectOnIORank_.globalAquiferData()
                   : std::move(localAquiferData)
    };

    if (eclState_.getSimulationConfig().useThresholdPressure()) {
        restartValue.addExtra("THRESHPR", UnitSystem::measure::pressure,
                              thresholdPressure);
    }

    // Add suggested next timestep to extra data.
    if (! isSubStep) {
        restartValue.addExtra("OPMEXTRA", std::vector<double>(1, nextStepSize));
    }

    // Add nnc flows and flores.
    if (isFlowsn) {
        const auto flowsn_global = isParallel ? this->collectOnIORank_.globalFlowsn() : std::move(flowsn);
        for (const auto& flows : flowsn_global) {
            if (flows.name.empty())
                continue;
            if (flows.name == "FLOGASN+") {
                restartValue.addExtra(flows.name, UnitSystem::measure::gas_surface_rate, flows.values);
            } else {
                restartValue.addExtra(flows.name, UnitSystem::measure::liquid_surface_rate, flows.values);
            }
        }
    }
    if (isFloresn) {
        const auto floresn_global = isParallel ? this->collectOnIORank_.globalFloresn() : std::move(floresn);
        for (const auto& flores : floresn_global) {
            if (flores.name.empty()) {
                continue;
            }
            restartValue.addExtra(flores.name, UnitSystem::measure::rate, flores.values);
        }
    }
    // make sure that the previous I/O request has been completed
    // and the number of incomplete tasklets does not increase between
    // time steps
    this->taskletRunner_->barrier();

    // check if there might have been a failure in the TaskletRunner
    if (this->taskletRunner_->failure()) {
        throw std::runtime_error("Failure in the TaskletRunner while writing output.");
    }

    // create a tasklet to write the data for the current time step to disk
    auto eclWriteTasklet = std::make_shared<EclWriteTasklet>(
        actionState,
        isParallel ? this->collectOnIORank_.globalWellTestState() : std::move(localWTestState),
        summaryState, udqState, *this->eclIO_,
        reportStepNum, timeStepNum, isSubStep, curTime, std::move(restartValue), doublePrecision);

    // finally, start a new output writing job
    this->taskletRunner_->dispatch(std::move(eclWriteTasklet));
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
evalSummary(const int                                            reportStepNum,
            const Scalar                                         curTime,
            const data::Wells&                                   localWellData,
            const data::WellBlockAveragePressures&               localWBPData,
            const data::GroupAndNetworkValues&                   localGroupAndNetworkData,
            const std::map<int,data::AquiferData>&               localAquiferData,
            const std::map<std::pair<std::string, int>, double>& blockData,
            const std::map<std::string, double>&                 miscSummaryData,
            const std::map<std::string, std::vector<double>>&    regionData,
            const Inplace&                                       inplace,
            const std::optional<Inplace>&                        initialInPlace,
            const InterRegFlowMap&                               interRegFlows,
            SummaryState&                                        summaryState,
            UDQState&                                            udqState)
{
    if (collectOnIORank_.isIORank()) {
        const auto& summary = eclIO_->summary();

        const auto& wellData = this->collectOnIORank_.isParallel()
            ? this->collectOnIORank_.globalWellData()
            : localWellData;

        const auto& wbpData = this->collectOnIORank_.isParallel()
            ? this->collectOnIORank_.globalWBPData()
            : localWBPData;

        const auto& groupAndNetworkData = this->collectOnIORank_.isParallel()
            ? this->collectOnIORank_.globalGroupAndNetworkData()
            : localGroupAndNetworkData;

        const auto& aquiferData = this->collectOnIORank_.isParallel()
            ? this->collectOnIORank_.globalAquiferData()
            : localAquiferData;

        summary.eval(summaryState,
                     reportStepNum,
                     curTime,
                     wellData,
                     wbpData,
                     groupAndNetworkData,
                     miscSummaryData,
                     initialInPlace,
                     inplace,
                     regionData,
                     blockData,
                     aquiferData,
                     getInterRegFlowsAsMap(interRegFlows));

        // Off-by-one-fun: The reportStepNum argument corresponds to the
        // report step these results will be written to, whereas the
        // argument to UDQ function evaluation corresponds to the report
        // step we are currently on.
        const auto udq_step = reportStepNum - 1;

        this->schedule_[udq_step].udq()
            .eval(udq_step,
                  this->schedule_.wellMatcher(udq_step),
                  this->schedule_[udq_step].group_order(),
                  this->schedule_.segmentMatcherFactory(udq_step),
                  [es = std::cref(this->eclState_)]() {
                      return std::make_unique<RegionSetMatcher>
                          (es.get().fipRegionStatistics());
                  },
                  summaryState, udqState);
    }

#if HAVE_MPI
    if (collectOnIORank_.isParallel()) {
        Parallel::MpiSerializer ser(grid_.comm());
        ser.append(summaryState);
    }
#endif
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
const typename EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::TransmissibilityType&
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
globalTrans() const
{
    assert (globalTrans_);
    return *globalTrans_;
}

} // namespace Opm

#endif // OPM_ECL_GENERIC_WRITER_IMPL_HPP
