/*
  Copyright 2026, Equinor ASA

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
#ifndef FLOWEXP_COMP_PH_HPP
#define FLOWEXP_COMP_PH_HPP

#include <opm/models/phflash/phflashmodel.hh>

#include <flowexperimental/comp/flowexp_comp.hpp>
#include <flowexperimental/comp/FlowProblemCompPh.hpp>

#include <tuple>

namespace Opm::Properties {
namespace TTag {

/*!
 * \brief The isenthalpic (P-H) variant of the experimental compositional
 *        simulator.
 *
 * The PT frontend comes first in the tuple so that everything it configures
 * (fluid system, well/aquifer/tracer models, discretization, deck problem
 * machinery) keeps its established precedence. The P-H specifics are
 * restated on this type tag below rather than relied upon from the
 * PHFlashModel parent, so they win independently of tuple order.
 */
template<int NumComp, bool EnableWater>
struct FlowExpCompPhProblem {
   using InheritsFrom = std::tuple<FlowExpCompProblem<NumComp, EnableWater>, PHFlashModel>;
};

}

// the P-H overrides, restated on this type tag (see the type tag brief):
// the isenthalpic solver, the intensive quantities that make temperature a
// flash output, and the model that registers the PhFlash* parameters

template <class TypeTag, int NumComp, bool EnableWater>
struct PhFlashSolver<TypeTag, TTag::FlowExpCompPhProblem<NumComp, EnableWater>>
{
    using type = PHFlash<GetPropType<TypeTag, Properties::Scalar>,
                         GetPropType<TypeTag, Properties::FluidSystem>>;
};

template <class TypeTag, int NumComp, bool EnableWater>
struct IntensiveQuantities<TypeTag, TTag::FlowExpCompPhProblem<NumComp, EnableWater>>
{ using type = PHFlashIntensiveQuantities<TypeTag>; };

template <class TypeTag, int NumComp, bool EnableWater>
struct Model<TypeTag, TTag::FlowExpCompPhProblem<NumComp, EnableWater>>
{ using type = PHFlashModel<TypeTag>; };

// the problem serving the P-H duck contract (specifiedEnthalpy(),
// phFlashConfig()) on top of the deck-driven compositional problem
template <class TypeTag, int NumComp, bool EnableWater>
struct Problem<TypeTag, TTag::FlowExpCompPhProblem<NumComp, EnableWater>>
{
    using type = FlowProblemCompPh<TypeTag>;
};

} // namespace Opm::Properties

#endif // FLOWEXP_COMP_PH_HPP
