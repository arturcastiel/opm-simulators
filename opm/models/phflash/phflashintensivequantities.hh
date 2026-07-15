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
 * \copydoc Opm::PHFlashIntensiveQuantities
 */
#ifndef OPM_PHFLASH_INTENSIVE_QUANTITIES_HH
#define OPM_PHFLASH_INTENSIVE_QUANTITIES_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <opm/common/Exceptions.hpp>

#include <opm/material/Constants.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <opm/models/common/energymodule.hh>
#include <opm/models/common/diffusionmodule.hh>

#include <opm/models/flash/flashproperties.hh>

#include <opm/models/ptflash/flashindices.hh>
#include <opm/models/phflash/phflashparameters.hh>

#include <array>
#include <iostream>
#include <string>

namespace Opm {

/*!
 * \ingroup PHFlashModel
 * \ingroup IntensiveQuantities
 *
 * \brief Intensive quantities of the isenthalpic (P-H) compositional model:
 *        per cell, the temperature is NOT taken from the primary variables
 *        or the problem — it is produced by the P-H flash, which inverts the
 *        mixture enthalpy for the temperature at which the flashed state
 *        attains the specified enthalpy.
 *
 * The specified enthalpy is per-dof state supplied by the problem through
 * the duck method
 *     Scalar specifiedEnthalpy(const ElementContext&, unsigned dofIdx,
 *                              unsigned timeIdx) const;   // molar [J/mol]
 * and the solver configuration through
 *     const typename PhFlashSolver::Config& phFlashConfig() const;
 * (built once from the PhFlash* runtime parameters — the same pattern as
 * problem.getEosType()).
 *
 * Value-level contract (deliberate, documented): the reconstructed
 * temperature carries NO derivatives w.r.t. the primary variables — the
 * P-H solver is values-only. This model stage targets explicit/diagnostic
 * and synthetic-case use, validated manually; wiring the temperature
 * sensitivities (dT/d{H,p,z} by the implicit-function theorem) into the AD
 * graph is the designated next stage and changes nothing in this class's
 * surface.
 */
template <class TypeTag>
class PHFlashIntensiveQuantities
    : public GetPropType<TypeTag, Properties::DiscIntensiveQuantities>
    , public DiffusionIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableDiffusion>() >
    , public EnergyIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableEnergy>() >
    , public GetPropType<TypeTag, Properties::FluxModule>::FluxIntensiveQuantities
{
    using ParentType = GetPropType<TypeTag, Properties::DiscIntensiveQuantities>;

    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    // primary variable indices (identical to the PT flash model: z, p, [Sw])
    enum { z0Idx = Indices::z0Idx };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    static constexpr bool enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>();
    static constexpr bool enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>();
    enum { dimWorld = GridView::dimensionworld };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { water0Idx = Indices::water0Idx};

    static constexpr bool waterEnabled = Indices::waterEnabled;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PhFlashSolver = GetPropType<TypeTag, Properties::PhFlashSolver>;

    using ComponentVector = Dune::FieldVector<Evaluation, numComponents>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using DiffusionIntensiveQuantities = ::Opm::DiffusionIntensiveQuantities<TypeTag, enableDiffusion>;
    using EnergyIntensiveQuantities = ::Opm::EnergyIntensiveQuantities<TypeTag, enableEnergy>;
    using FluxIntensiveQuantities = typename FluxModule::FluxIntensiveQuantities;

public:
    //! The type of the object returned by the fluidState() method
    using FluidState = CompositionalFluidState<Evaluation, FluidSystem, enableEnergy>;

    PHFlashIntensiveQuantities() = default;

    PHFlashIntensiveQuantities(const PHFlashIntensiveQuantities& other) = default;

    PHFlashIntensiveQuantities& operator=(const PHFlashIntensiveQuantities& other) = default;

    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);
        // NOTE: deliberately NO EnergyIntensiveQuantities::updateTemperatures_
        // here — in this model the temperature is an OUTPUT of the P-H flash,
        // not an input from the primary variables or the problem.

        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const auto& problem = elemCtx.problem();

        const Scalar ptTolerance = Parameters::Get<Parameters::FlashTolerance<Scalar>>();
        const int flashVerbosity = Parameters::Get<Parameters::FlashVerbosity>();
        const std::string twoPhaseMethod = Parameters::Get<Parameters::FlashTwoPhaseMethod>();

        // overall composition z from the primary variables (identical to the
        // PT flash model)
        ComponentVector z(0.);
        {
            Evaluation lastZ = 1.0;
            for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
                z[compIdx] = priVars.makeEvaluation(z0Idx + compIdx, timeIdx);
                lastZ -= z[compIdx];
            }
            z[numComponents - 1] = lastZ;

            Evaluation sumz = 0.0;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                z[compIdx] = max(z[compIdx], 1e-8);
                sumz += z[compIdx];
            }
            z /= sumz;
        }

        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            fluidState_.setMoleFraction(compIdx, z[compIdx]);
        }

        const Evaluation p = priVars.makeEvaluation(pressure0Idx, timeIdx);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fluidState_.setPressure(phaseIdx, p);
        }

        /////////////
        // The P-H flash: invert the mixture enthalpy for the temperature.
        // The solver seeds Wilson K-values and L internally at every trial
        // temperature (its documented input contract), so no K/L hints are
        // wired here; hint-based warm starts are a later optimization.
        /////////////
        const Scalar hSpec = problem.specifiedEnthalpy(elemCtx, dofIdx, timeIdx); // molar [J/mol]
        const auto& cfg = problem.phFlashConfig();
        const auto& eos_type = problem.getEosType();

        if (flashVerbosity >= 1) {
            const int spatialIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
            std::cout << " updating the P-H intensive quantities for Cell " << spatialIdx
                      << " (hSpec = " << hSpec << " J/mol)" << std::endl;
        }

        const bool ok = PhFlashSolver::solve(fluidState_, hSpec, cfg,
                                             twoPhaseMethod, ptTolerance, eos_type,
                                             flashVerbosity);
        if (!ok) {
            const int spatialIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
            throw NumericalProblem("P-H flash failed for cell " + std::to_string(spatialIdx)
                                   + ": specified enthalpy " + std::to_string(hSpec)
                                   + " J/mol has no solution on the temperature bracket");
        }

        /////////////
        // Everything below is identical to the PT flash model: the flashed,
        // temperature-consistent state feeds Z-factors, saturations,
        // viscosities and densities unchanged.
        /////////////
        typename FluidSystem::template ParameterCache<Evaluation> paramCache(eos_type);
        paramCache.updatePhase(fluidState_, FluidSystem::oilPhaseIdx);

        const Scalar R = Opm::Constants<Scalar>::R;
        const Evaluation Z_L = (paramCache.molarVolume(FluidSystem::oilPhaseIdx) *
                                fluidState_.pressure(FluidSystem::oilPhaseIdx)) /
                               (R * fluidState_.temperature(FluidSystem::oilPhaseIdx));
        paramCache.updatePhase(fluidState_, FluidSystem::gasPhaseIdx);
        const Evaluation Z_V = (paramCache.molarVolume(FluidSystem::gasPhaseIdx) *
                                fluidState_.pressure(FluidSystem::gasPhaseIdx)) /
                               (R * fluidState_.temperature(FluidSystem::gasPhaseIdx));

        Evaluation Sw = 0.0;
        if constexpr (waterEnabled) {
            Sw = priVars.makeEvaluation(water0Idx, timeIdx);
        }
        const Evaluation L = fluidState_.L();
        Evaluation So = max((1 - Sw) * (L * Z_L / ( L * Z_L + (1 - L) * Z_V)), 0.0);
        Evaluation Sg = max(1 - So - Sw, 0.0);
        const Scalar sumS = getValue(So) + getValue(Sg) + getValue(Sw);
        So /= sumS;
        Sg /= sumS;

        fluidState_.setSaturation(FluidSystem::oilPhaseIdx, So);
        fluidState_.setSaturation(FluidSystem::gasPhaseIdx, Sg);
        if constexpr (waterEnabled) {
            Sw /= sumS;
            fluidState_.setSaturation(FluidSystem::waterPhaseIdx, Sw);
        }

        fluidState_.setCompressFactor(FluidSystem::oilPhaseIdx, Z_L);
        fluidState_.setCompressFactor(FluidSystem::gasPhaseIdx, Z_V);

        const MaterialLawParams& materialParams = problem.materialLawParams(elemCtx, dofIdx, timeIdx);

        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams, fluidState_);
        Valgrind::CheckDefined(relativePermeability_);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (phaseIdx == static_cast<unsigned int>(FluidSystem::oilPhaseIdx) ||
                phaseIdx == static_cast<unsigned int>(FluidSystem::gasPhaseIdx))
            {
                paramCache.updatePhase(fluidState_, phaseIdx);
            }

            const Evaluation& mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);

            fluidState_.setViscosity(phaseIdx, mu);

            mobility_[phaseIdx] = relativePermeability_[phaseIdx] / mu;
            Valgrind::CheckDefined(mobility_[phaseIdx]);

            const Evaluation& rho = FluidSystem::density(fluidState_, paramCache, phaseIdx);
            fluidState_.setDensity(phaseIdx, rho);
        }

        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        Valgrind::CheckDefined(porosity_);

        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        EnergyIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        DiffusionIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::intrinsicPermeability
     */
    const DimMatrix& intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::relativePermeability
     */
    const Evaluation& relativePermeability(unsigned phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::mobility
     */
    const Evaluation& mobility(unsigned phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

private:
    DimMatrix intrinsicPerm_;
    FluidState fluidState_;
    Evaluation porosity_;
    std::array<Evaluation,numPhases> relativePermeability_;
    std::array<Evaluation,numPhases> mobility_;
};

} // namespace Opm

#endif
