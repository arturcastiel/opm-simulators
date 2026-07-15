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
/*!
 * \file
 *
 * \copydoc Opm::PHFlashModel
 */
#ifndef OPM_PHFLASH_MODEL_HH
#define OPM_PHFLASH_MODEL_HH

#include <opm/material/constraintsolvers/PHFlash.hpp>

#include <opm/models/io/vtkphflashmodule.hpp>

#include <opm/models/ptflash/flashmodel.hh>

#include <opm/models/phflash/phflashintensivequantities.hh>
#include <opm/models/phflash/phflashparameters.hh>

#include <memory>
#include <tuple>

namespace Opm {

template <class TypeTag>
class PHFlashModel;

}

namespace Opm::Properties {

namespace TTag {

/*!
 * \brief The type tag for the isenthalpic (P-H) compositional model.
 *
 * Inherits everything from the PT flash model type tag (indices, primary
 * variables, local residual, Newton method, parameters) and overrides only
 * what the P-H specification changes: the flash solver slot (a SEPARATE
 * PhFlashSolver property — the base FlashSolver keeps its PT meaning for the
 * inner solve), the intensive quantities (where the temperature becomes an
 * output of the flash), and the model class. The primary-variable set is
 * UNCHANGED: z, p, [Sw].
 */
struct PHFlashModel { using InheritsFrom = std::tuple<FlashModel>; };

} // namespace TTag

//! Bind the isenthalpic flash solver
template<class TypeTag>
struct PhFlashSolver<TypeTag, TTag::PHFlashModel>
{
    using type = PHFlash<GetPropType<TypeTag, Properties::Scalar>,
                         GetPropType<TypeTag, Properties::FluidSystem>>;
};

//! the IntensiveQuantities property: temperature comes from the P-H flash
template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::PHFlashModel>
{ using type = PHFlashIntensiveQuantities<TypeTag>; };

//! the Model property
template<class TypeTag>
struct Model<TypeTag, TTag::PHFlashModel>
{ using type = PHFlashModel<TypeTag>; };

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup PHFlashModel
 *
 * \brief A compositional multi-phase model where the per-cell temperature is
 *        determined by an isenthalpic (P-H) flash.
 *
 * Identical to the PT flash model (same primary variables and conservation
 * equations) except that the flash specification is (pressure, molar
 * enthalpy, composition): at every intensive-quantities update the
 * temperature is found such that the flashed mixture attains the enthalpy
 * the problem specifies for that degree of freedom. The specified enthalpy
 * and the solver configuration are served by the problem
 * (specifiedEnthalpy(), phFlashConfig()); the inner isothermal solve keeps
 * the PT flash runtime parameters.
 */
template <class TypeTag>
class PHFlashModel
    : public FlashModel<TypeTag>
{
    using ParentType = FlashModel<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    explicit PHFlashModel(GetPropType<TypeTag, Properties::Simulator>& simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Register the runtime parameters of the P-H model (on top of the
     *        PT flash model's).
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        // register runtime parameters of the P-H diagnostic VTK module
        VtkPhFlashModule<TypeTag>::registerParameters();

        Parameters::Register<Parameters::PhFlashEnthalpyModel>
            ("Enthalpy model inverted by the P-H flash: caloric, eos_departure");
        Parameters::Register<Parameters::PhFlashTempMin<Scalar>>
            ("Lower bound of the P-H flash temperature bracket [K]");
        Parameters::Register<Parameters::PhFlashTempMax<Scalar>>
            ("Upper bound of the P-H flash temperature bracket [K]");
        Parameters::Register<Parameters::PhFlashTolerance<Scalar>>
            ("Convergence tolerance of the P-H flash bracketing solver "
             "(temperature interval [K] and enthalpy residual [J/mol])");
        Parameters::Register<Parameters::PhFlashMaxIterations>
            ("Maximum iterations of the P-H flash bracketing solver");
    }

    // public like the PT flash model's: the discretization invokes this
    // through the CRTP implementation reference
    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        // the P-H diagnostics: reconstructed temperature, mixture enthalpy,
        // specified enthalpy, enthalpy residual (gated by --vtk-write-ph-flash)
        this->addOutputModule(std::make_unique<VtkPhFlashModule<TypeTag>>(this->simulator_));
    }
};

} // namespace Opm

#endif
