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
 * \copydoc Opm::VtkPhFlashModule
 */
#ifndef OPM_VTK_PHFLASH_MODULE_HPP
#define OPM_VTK_PHFLASH_MODULE_HPP

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/constraintsolvers/MixtureEnthalpy.hpp>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/io/vtkphflashparams.hpp>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the isenthalpic (P-H) flash diagnostics.
 *
 * This module writes the quantities a manual validation of a P-H run needs
 * per cell:
 *   - flashTemperature [K]: the temperature reconstructed by the P-H flash
 *   - mixtureEnthalpy [J/mol]: molar enthalpy of the flashed state under the
 *     configured enthalpy model
 *   - specifiedEnthalpy [J/mol]: the enthalpy the problem prescribes for the
 *     degree of freedom
 *   - enthalpyResidual [J/mol]: mixtureEnthalpy - specifiedEnthalpy (signed);
 *     values beyond the flash tolerance indicate a diagnosable defect
 *
 * It relies on the P-H problem duck methods (specifiedEnthalpy(),
 * phFlashConfig(), getEosType()), so it is only meaningful — and only
 * registered — for models deriving from the PHFlashModel type tag.
 */
template <class TypeTag>
class VtkPhFlashModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    using BufferType = typename ParentType::BufferType;
    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    explicit VtkPhFlashModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        params_.read();
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        VtkPhFlashParams::registerParameters();
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers() override
    {
        if (params_.phFlashOutput_) {
            this->resizeScalarBuffer_(flashTemperature_, BufferType::Dof);
            this->resizeScalarBuffer_(mixtureEnthalpy_, BufferType::Dof);
            this->resizeScalarBuffer_(specifiedEnthalpy_, BufferType::Dof);
            this->resizeScalarBuffer_(enthalpyResidual_, BufferType::Dof);
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant
     *        for an element
     */
    void processElement(const ElementContext& elemCtx) override
    {
        using Toolbox = MathToolbox<Evaluation>;

        if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
            return;
        }
        if (!params_.phFlashOutput_) {
            return;
        }

        const auto& problem = elemCtx.problem();
        const auto& cfg = problem.phFlashConfig();
        const auto& eosType = problem.getEosType();

        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            const unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            const Scalar hMix =
                Toolbox::value(MixtureEnthalpy<Scalar, FluidSystem>::mixtureEnthalpy(
                    fs, cfg.cpTable, cfg.refTemperature, eosType, cfg.model));
            const Scalar hSpec = problem.specifiedEnthalpy(elemCtx, i, /*timeIdx=*/0);

            flashTemperature_[I] = Toolbox::value(fs.temperature(/*phaseIdx=*/0));
            mixtureEnthalpy_[I] = hMix;
            specifiedEnthalpy_[I] = hSpec;
            enthalpyResidual_[I] = hMix - hSpec;
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter) override
    {
        if (!dynamic_cast<VtkMultiWriter*>(&baseWriter)) {
            return;
        }

        if (params_.phFlashOutput_) {
            this->commitScalarBuffer_(baseWriter, "flashTemperature", flashTemperature_, BufferType::Dof);
            this->commitScalarBuffer_(baseWriter, "mixtureEnthalpy", mixtureEnthalpy_, BufferType::Dof);
            this->commitScalarBuffer_(baseWriter, "specifiedEnthalpy", specifiedEnthalpy_, BufferType::Dof);
            this->commitScalarBuffer_(baseWriter, "enthalpyResidual", enthalpyResidual_, BufferType::Dof);
        }
    }

private:
    VtkPhFlashParams params_{};
    ScalarBuffer flashTemperature_{};
    ScalarBuffer mixtureEnthalpy_{};
    ScalarBuffer specifiedEnthalpy_{};
    ScalarBuffer enthalpyResidual_{};
};

} // namespace Opm

#endif // OPM_VTK_PHFLASH_MODULE_HPP
