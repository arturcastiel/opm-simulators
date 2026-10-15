/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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

#ifndef OPM_BLACKOILWELLMODEL_WBP_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_WBP_HEADER_INCLUDED

#include <opm/input/eclipse/EclipseState/Grid/GridDims.hpp>

#include <opm/output/data/Wells.hpp>

#include <opm/simulators/wells/ParallelWBPCalculation.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <cstddef>
#include <map>
#include <optional>
#include <vector>

namespace Opm {

template<typename Scalar, typename IndexTraits> class BlackoilWellModelGeneric;

/// Class for handling the blackoil well model.
template<typename Scalar, typename IndexTraits>
class BlackoilWellModelWBP
{
public:
    explicit BlackoilWellModelWBP(BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model);

    void initializeSources(typename ParallelWBPCalculation<Scalar>::GlobalToLocal index,
                           typename ParallelWBPCalculation<Scalar>::Evaluator eval);

    void registerOpenWellsForWBPCalculation();

    typename ParallelWBPCalculation<Scalar>::EvaluatorFactory
    makeWellSourceEvaluatorFactory(const std::vector<Well>::size_type wellIdx) const;

    void initializeWBPCalculationService();

    data::WellBlockAveragePressures
    computeWellBlockAveragePressures(const Scalar gravity) const;

private:
    /// Calculation service of the grid a well's connections are in: the
    /// global grid (0) or an LGR (its number in the deck, from 1).
    ParallelWBPCalculation<Scalar>& wbpCalculationService(int grid) const;

    BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model_;
    mutable ParallelWBPCalculation<Scalar> wbpCalculationService_;

    /// Dimensions of each LGR, keyed by the LGR's number in the deck.  The
    /// LGR's calculation service refers to them.
    std::map<int, GridDims> lgrDims_{};

    /// One calculation service per LGR, keyed by the LGR's number in the
    /// deck.  It identifies the LGR's cells by their Cartesian index in the
    /// LGR, so that a connection's neighbours are those of the LGR.
    mutable std::map<int, ParallelWBPCalculation<Scalar>> lgrWbpCalculationServices_{};

    struct WBPCalcID
    {
        std::optional<typename std::vector<WellInterfaceGeneric<Scalar, IndexTraits>*>::size_type> openWellIdx_{};
        std::size_t wbpCalcIdx_{};
        int grid_{0};
    };

    std::vector<WBPCalcID> wbpCalcMap_{};
};


} // namespace Opm

#endif
