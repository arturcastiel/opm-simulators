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
/*!
 * \file
 * \copydoc Opm::FlowProblemCompPh
 */
#ifndef OPM_FLOW_PROBLEM_COMP_PH_HPP
#define OPM_FLOW_PROBLEM_COMP_PH_HPP

#include <dune/common/fvector.hh>

#include <opm/common/Exceptions.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/constraintsolvers/IdealGasCaloricData.hpp>
#include <opm/material/constraintsolvers/MixtureEnthalpy.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <opm/models/ptflash/flashparameters.hh>
#include <opm/models/phflash/phflashparameters.hh>

#include <opm/simulators/flow/FlowProblemComp.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

namespace Opm {

/*!
 * \ingroup PHFlashModel
 *
 * \brief Compositional flow problem for the isenthalpic (P-H) model: serves
 *        the per-cell specified enthalpy and the P-H solver configuration
 *        that PHFlashIntensiveQuantities requires.
 *
 * The specified enthalpy is derived once, at the end of finishInit(), from
 * the deck-given initial state: each cell's (p, T, z) is flashed
 * isothermally and the molar mixture enthalpy of the flashed state becomes
 * that cell's enthalpy specification. The P-H model therefore starts from an
 * enthalpy field that is exactly consistent with the deck temperature — at
 * time zero the flash must reconstruct the deck temperature, which makes
 * |flashTemperature - TEMPI| (and the enthalpyResidual VTK field) a direct,
 * per-cell validation of the whole wiring. Prescribing an enthalpy field
 * that is NOT derived from a temperature field (e.g. from an energy-balance
 * source) is the designated next stage and only replaces the derivation in
 * computeSpecifiedEnthalpies_().
 *
 * The caloric data has a single derivation point here: one heat-capacity
 * table is built from the registered component names and handed to both the
 * fluid system (the energy path's enthalpy()) and the P-H flash config —
 * both paths invert the same enthalpy by construction.
 */
template <class TypeTag>
class FlowProblemCompPh : public FlowProblemComp<TypeTag>
{
    using FlowProblemCompType = FlowProblemComp<TypeTag>;

    // taken from the property system rather than re-exported from the base:
    // FlowProblemComp keeps its own aliases private
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    using PhFlashSolver = GetPropType<TypeTag, Properties::PhFlashSolver>;
    using PhFlashConfig = typename PhFlashSolver::Config;

public:
    explicit FlowProblemCompPh(Simulator& simulator)
        : FlowProblemCompType(simulator)
    {}

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        FlowProblemCompType::finishInit();

        initEnthalpyData_();
        computeSpecifiedEnthalpies_();
    }

    /*!
     * \brief The molar enthalpy [J/mol] the P-H flash must attain for a
     *        degree of freedom (duck method of PHFlashIntensiveQuantities).
     *
     * Expressed against the caloric reference datum of the P-H config
     * (IdealGasCaloricData: H(T0) = 0).
     */
    template <class Context>
    Scalar specifiedEnthalpy(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return specifiedEnthalpy_[context.globalSpaceIndex(spaceIdx, timeIdx)]; }

    /*!
     * \brief The P-H solver configuration, built once from the PhFlash*
     *        runtime parameters (duck method of PHFlashIntensiveQuantities;
     *        same pattern as getEosType()).
     */
    const PhFlashConfig& phFlashConfig() const
    { return phFlashConfig_; }

private:
    /*!
     * \brief Single derivation point of the caloric data: build the
     *        heat-capacity table from the registered component names and
     *        hand the SAME table to the fluid system and the P-H config.
     *
     * IdealGasCaloricData::byName() throws naming the component when no
     * preset exists — an unknown component is an input error, never a
     * silent fallback.
     */
    void initEnthalpyData_()
    {
        CpTable<Scalar, numComponents> table;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            table[compIdx] =
                IdealGasCaloricData<Scalar>::byName(std::string(FluidSystem::componentName(compIdx)));
        }
        FluidSystem::setEnthalpyData(table);

        phFlashConfig_.cpTable = table;
        phFlashConfig_.model =
            enthalpyModelFromString(Parameters::Get<Parameters::PhFlashEnthalpyModel>());
        phFlashConfig_.tempMin = Parameters::Get<Parameters::PhFlashTempMin<Scalar>>();
        phFlashConfig_.tempMax = Parameters::Get<Parameters::PhFlashTempMax<Scalar>>();
        phFlashConfig_.tolerance = Parameters::Get<Parameters::PhFlashTolerance<Scalar>>();
        phFlashConfig_.maxIterations = Parameters::Get<Parameters::PhFlashMaxIterations>();
        // refTemperature keeps the IdealGasCaloricData datum — the same
        // datum setEnthalpyData() defaults to above
    }

    /*!
     * \brief Derive each cell's enthalpy specification from its deck-given
     *        initial state: flash (p, T, z) isothermally, then evaluate the
     *        molar mixture enthalpy of the flashed state.
     *
     * A cell whose initial state cannot be flashed is an input error and
     * throws, naming the cell.
     */
    void computeSpecifiedEnthalpies_()
    {
        // PTFlash::solve()'s derivative update is written against an AD
        // value type (numVars/derivative()/setDerivative()) and does not
        // compile for a plain-Scalar fluid state, so the flash runs on the
        // model's own Evaluation type — the same instantiation the intensive
        // quantities use — and only the values are kept.
        using FlashFluidState = CompositionalFluidState<Evaluation, FluidSystem>;
        using PtFlash = PTFlash<Scalar, FluidSystem>;

        const auto& eosType = this->getEosType();
        const Scalar ptTolerance = Parameters::Get<Parameters::FlashTolerance<Scalar>>();
        const std::string twoPhaseMethod = Parameters::Get<Parameters::FlashTwoPhaseMethod>();

        const auto& initialStates = this->initialFluidStates();
        specifiedEnthalpy_.resize(initialStates.size());

        for (std::size_t dofIdx = 0; dofIdx < initialStates.size(); ++dofIdx) {
            const auto& initialFs = initialStates[dofIdx];
            const auto z = initialTotalMoleFractions_(initialFs, eosType);

            FlashFluidState fs;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                fs.setPressure(phaseIdx, initialFs.pressure(phaseIdx));
            }
            fs.setTemperature(initialFs.temperature(/*phaseIdx=*/0));
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                fs.setMoleFraction(compIdx, z[compIdx]);
            }

            // Wilson seed, as in FlowProblemComp::initial(): the isothermal
            // flash refines K and L from here
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                fs.setKvalue(compIdx, fs.wilsonK_(compIdx));
            }
            fs.setLvalue(-1.0);

            try {
                PtFlash::solve(fs, twoPhaseMethod, ptTolerance, eosType);
            }
            catch (const std::exception& e) {
                throw NumericalProblem(
                    "P-H initialization: isothermal flash of the initial state failed for cell "
                    + std::to_string(dofIdx) + " (" + e.what()
                    + "); check the EQUIL/TEMPI/ZMF input");
            }

            specifiedEnthalpy_[dofIdx] =
                getValue(MixtureEnthalpy<Scalar, FluidSystem>::mixtureEnthalpy(
                    fs, phFlashConfig_.cpTable, phFlashConfig_.refTemperature,
                    eosType, phFlashConfig_.model));
        }
    }

    /*!
     * \brief The total mole fractions z of an initial fluid state, mirroring
     *        the derivation FlowProblemComp::initial() uses for the primary
     *        variables.
     *
     * With ZMF initialization the deck gives z directly; otherwise z is
     * accumulated from the per-phase compositions weighted by saturation and
     * EoS molar density. All arithmetic is plain-Scalar: this feeds values,
     * not derivatives.
     */
    template <class InitialFluidState, class EosType>
    Dune::FieldVector<Scalar, numComponents>
    initialTotalMoleFractions_(const InitialFluidState& initialFs,
                               const EosType& eosType) const
    {
        Dune::FieldVector<Scalar, numComponents> z(0.0);

        if (this->zmfInitialization()) {
            // apply the IDENTICAL transformation the intensive quantities
            // apply when they rebuild z from the primary variables (last
            // component as the complement, clamp, renormalize): a deck ZMF
            // that does not sum exactly to one would otherwise flash a
            // different composition for H_spec than the model later
            // inverts, producing a systematic t=0 temperature offset with
            // a small residual — silently defeating the flashTemperature
            // vs. deck-temperature acceptance check
            Scalar lastZ = 1.0;
            for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
                z[compIdx] = initialFs.moleFraction(compIdx);
                lastZ -= z[compIdx];
            }
            z[numComponents - 1] = lastZ;

            Scalar sumz = 0.0;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                z[compIdx] = std::max(z[compIdx], Scalar{1e-8});
                sumz += z[compIdx];
            }
            z /= sumz;
            return z;
        }

        CompositionalFluidState<Scalar, FluidSystem> fs;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fs.setPressure(phaseIdx, initialFs.pressure(phaseIdx));
            fs.setSaturation(phaseIdx, initialFs.saturation(phaseIdx));
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                fs.setMoleFraction(phaseIdx, compIdx, initialFs.moleFraction(phaseIdx, compIdx));
            }
        }
        fs.setTemperature(initialFs.temperature(/*phaseIdx=*/0));

        typename FluidSystem::template ParameterCache<Scalar> paramCache(eosType);
        paramCache.updatePhase(fs, oilPhaseIdx);
        paramCache.updatePhase(fs, gasPhaseIdx);
        fs.setDensity(oilPhaseIdx, FluidSystem::density(fs, paramCache, oilPhaseIdx));
        fs.setDensity(gasPhaseIdx, FluidSystem::density(fs, paramCache, gasPhaseIdx));

        Scalar sumMoles = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (Indices::waterEnabled && phaseIdx == static_cast<unsigned>(waterPhaseIdx)) {
                continue;
            }
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar tmp = fs.molarity(phaseIdx, compIdx) * fs.saturation(phaseIdx);
                tmp = std::max(tmp, Scalar{1e-8});
                z[compIdx] += tmp;
                sumMoles += tmp;
            }
        }
        z /= sumMoles;

        return z;
    }

    PhFlashConfig phFlashConfig_{};
    std::vector<Scalar> specifiedEnthalpy_{};
};

} // namespace Opm

#endif // OPM_FLOW_PROBLEM_COMP_PH_HPP
