// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 *
 * \brief External-reference benchmark of the isenthalpic (P-H) flash
 *        machinery against CoolProp-generated states.
 *
 * Five methane/n-decane cases were evaluated OFFLINE with CoolProp (see the
 * per-case provenance blocks below); their results are embedded here as
 * constants. CoolProp is NOT a dependency of this test — only its numbers
 * appear, with full provenance.
 *
 * The test SHOWS the comparison (a per-case table on stdout; run the binary
 * directly for the full printout) and asserts only STRUCTURAL metrics.
 * Absolute enthalpy agreement is deliberately NOT asserted: OPM's
 * Peng-Robinson + polynomial-cp caloric model and CoolProp's reference
 * Helmholtz (HEOS) or PR formulations are different models, so a value
 * tolerance would either hide defects (too loose) or fail on honest model
 * differences (too tight). The printed table shows the magnitudes; the
 * CoolProp PR-vs-HEOS spread is printed alongside as the natural yardstick.
 *
 * Asserted metrics:
 *   M1 regime agreement   — OPM finds two-phase exactly where the reference
 *                           does (phase LABELS are not compared: at 410 K
 *                           methane is far supercritical and CoolProp's own
 *                           backends disagree on the label).
 *   M2 monotonicity       — OPM H(T) is strictly increasing around every
 *                           case point (cp > 0), as the reference data is.
 *   M3 rank agreement     — the five case enthalpies sort in the same order
 *                           on the OPM side as in the reference.
 *   M4 round-trip closure — flash at T, take H, invert with PHFlash: T is
 *                           recovered within the solver's own tolerance (an
 *                           internal-consistency bound, not a cross-model
 *                           tolerance).
 */
#include "config.h"

#define BOOST_TEST_MODULE PhFlashReference
#include <boost/test/unit_test.hpp>

#include <opm/material/components/C1.hpp>
#include <opm/material/components/C10.hpp>
#include <opm/material/constraintsolvers/IdealGasCaloricData.hpp>
#include <opm/material/constraintsolvers/MixtureEnthalpy.hpp>
#include <opm/material/constraintsolvers/PHFlash.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/eos/CubicEOS.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/BaseFluidSystem.hpp>
#include <opm/material/fluidsystems/PTFlashParameterCache.hpp>
#include <opm/material/viscositymodels/ViscosityModels.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <opm/material/common/MathToolbox.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

/*!
 * \brief The methane/n-decane fixture fluid system, mirroring the "F1"
 *        fixture of the opm-common ph_mvp suite (tests/material/
 *        ph_mvp_fixtures.hh) so both sides of the project benchmark the
 *        SAME parameterization: C1/C10 component classes, kij = 0.0411.
 *
 * Test-local by design; opm-common test headers are not installed, hence
 * the mirror rather than an include.
 */
template<class Scalar>
class TwoComponentFluidSystem
    : public Opm::BaseFluidSystem<Scalar, TwoComponentFluidSystem<Scalar>>
{
public:
    static constexpr int numPhases = 2;
    static constexpr int numComponents = 2;
    static constexpr int numMisciblePhases = 2;
    static constexpr int numMiscibleComponents = 2;
    static constexpr bool waterEnabled = false;
    static constexpr int oilPhaseIdx = 0;
    static constexpr int gasPhaseIdx = 1;
    static constexpr int waterPhaseIdx = -1;

    using Comp0 = Opm::C1<Scalar>;
    using Comp1 = Opm::C10<Scalar>;

    //! C1/nC10 binary interaction parameter (standard PR literature value;
    //! the unequal-index shortcut is valid only because the system is binary)
    static constexpr Scalar bipC1nC10 = 0.0411;

    template <class ValueType>
    using ParameterCache = Opm::PTFlashParameterCache<ValueType, TwoComponentFluidSystem<Scalar>>;
    using ViscosityModel = Opm::ViscosityModels<Scalar, TwoComponentFluidSystem<Scalar>>;
    using CubicEOS = Opm::CubicEOS<Scalar, TwoComponentFluidSystem<Scalar>>;

    static bool phaseIsActive(unsigned phaseIdx)
    { return phaseIdx == oilPhaseIdx || phaseIdx == gasPhaseIdx; }

    static Scalar acentricFactor(unsigned compIdx)
    {
        switch (compIdx) {
        case 0: return Comp0::acentricFactor();
        case 1: return Comp1::acentricFactor();
        default: throw std::runtime_error("Illegal component index for acentricFactor");
        }
    }

    static Scalar criticalTemperature(unsigned compIdx)
    {
        switch (compIdx) {
        case 0: return Comp0::criticalTemperature();
        case 1: return Comp1::criticalTemperature();
        default: throw std::runtime_error("Illegal component index for criticalTemperature");
        }
    }

    static Scalar criticalPressure(unsigned compIdx)
    {
        switch (compIdx) {
        case 0: return Comp0::criticalPressure();
        case 1: return Comp1::criticalPressure();
        default: throw std::runtime_error("Illegal component index for criticalPressure");
        }
    }

    static Scalar criticalVolume(unsigned compIdx)
    {
        switch (compIdx) {
        case 0: return Comp0::criticalVolume();
        case 1: return Comp1::criticalVolume();
        default: throw std::runtime_error("Illegal component index for criticalVolume");
        }
    }

    static Scalar molarMass(unsigned compIdx)
    {
        switch (compIdx) {
        case 0: return Comp0::molarMass();
        case 1: return Comp1::molarMass();
        default: throw std::runtime_error("Illegal component index for molarMass");
        }
    }

    static Scalar interactionCoefficient(unsigned comp1Idx, unsigned comp2Idx)
    { return comp1Idx != comp2Idx ? bipC1nC10 : 0.0; }

    static std::string_view phaseName(unsigned phaseIdx)
    {
        static const std::string_view name[] = {"o", "g"};
        assert(phaseIdx < 2);
        return name[phaseIdx];
    }

    static std::string_view componentName(unsigned compIdx)
    {
        static const std::string_view name[] = {Comp0::name(), Comp1::name()};
        assert(compIdx < 2);
        return name[compIdx];
    }

    template <class FluidState, class LhsEval = typename FluidState::ValueType, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& paramCache,
                           unsigned phaseIdx)
    {
        assert(phaseIdx == oilPhaseIdx || phaseIdx == gasPhaseIdx);
        return Opm::decay<LhsEval>(fluidState.averageMolarMass(phaseIdx) /
                                   paramCache.molarVolume(phaseIdx));
    }

    template <class FluidState, class LhsEval = typename FluidState::ValueType, class ParamCacheEval = LhsEval>
    static LhsEval viscosity(const FluidState& fluidState,
                             const ParameterCache<ParamCacheEval>& paramCache,
                             unsigned phaseIdx)
    { return Opm::decay<LhsEval>(ViscosityModel::LBC(fluidState, paramCache, phaseIdx)); }

    template <class FluidState, class LhsEval = typename FluidState::ValueType, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& paramCache,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(phaseIdx < numPhases);
        assert(compIdx < numComponents);
        return Opm::decay<LhsEval>(CubicEOS::computeFugacityCoefficient(
            fluidState, paramCache, phaseIdx, compIdx));
    }

    static bool isCompressible([[maybe_unused]] unsigned phaseIdx)
    { return true; }

    static bool isIdealMixture([[maybe_unused]] unsigned phaseIdx)
    { return false; }

    static bool isLiquid(unsigned phaseIdx)
    { return phaseIdx == 0; }

    static bool isIdealGas(unsigned phaseIdx)
    { return phaseIdx == 1; }
};

using FluidSystem = TwoComponentFluidSystem<double>;
//! arity bound to the PTFlash slot layout, as in the opm-common suite
using FlashEval = Opm::DenseAd::Evaluation<double, FluidSystem::numComponents + 1>;
using FluidState = Opm::CompositionalFluidState<FlashEval, FluidSystem>;
using PtFlash = Opm::PTFlash<double, FluidSystem>;
using PhFlash = Opm::PHFlash<double, FluidSystem>;
using Enthalpy = Opm::MixtureEnthalpy<double, FluidSystem>;
using EOSType = Opm::CompositionalConfig::EOSType;

constexpr double T0 = 298.15;             // OPM caloric datum [K]
constexpr const char* TWO_PHASE_METHOD = "ssi";
constexpr double PT_TOLERANCE = 1e-8;

/*!
 * \brief One embedded reference case.
 *
 * PROVENANCE (all cases): CoolProp — Bell, Wronski, Quoilin & Lemort (2014),
 * Ind. Eng. Chem. Res. 53(6), 2498-2508, doi:10.1021/ie4033999. Generated
 * with the CoolProp 8.0.0 python wheel (source clone 23b7323a for citation);
 * mixture "Methane&n-Decane". Two backends per case: HEOS (reference
 * multiparameter Helmholtz mixture model, Kunz & Wagner GERG-2008 form) and
 * PR (Peng-Robinson with kij = 0.0411 — the SAME value this fixture uses,
 * isolating implementation differences from the EoS model gap).
 * Fluid EoS papers (from CoolProp's BibTeX library at generation time):
 * methane — Setzmann & Wagner, J. Phys. Chem. Ref. Data 20 (1991) 1061;
 * n-decane — Lemmon & Span, J. Chem. Eng. Data 51 (2006) 785.
 *
 * DATUM: CoolProp enthalpies were shifted onto the OPM caloric datum
 * h(T0 = 298.15 K) = 0 by subtracting the mixture enthalpy at the near-ideal
 * anchor state (100 Pa, T0) of the same backend and composition (the two
 * backends agree there to 0.002%). The shifted values below are therefore
 * directly comparable to this test's mixtureEnthalpy() output.
 */
struct ReferenceCase {
    std::string_view name;
    std::string_view regime;
    double p;                 // [Pa]
    double t;                 // [K]
    std::array<double, 2> z;  // mole fractions (C1, nC10)
    bool twoPhase;            // reference phase regime
    double hHeos;             // shifted molar enthalpy, HEOS backend [J/mol]
    double hPr;               // shifted molar enthalpy, PR backend   [J/mol]
    double liqFracHeos;       // liquid molar fraction 1-Q, HEOS (two-phase only)
    double liqFracPr;         // liquid molar fraction 1-Q, PR   (two-phase only)
};

// generated 2026-07-15 (project record: isenthalpic-flash-imp/phases/pr/
// validation — generator, sweeps and citations live there, outside OPM)
constexpr std::array<ReferenceCase, 5> REFERENCE_CASES{{
    {.name = "case1", .regime = "two-phase (F1 anchor)",
     .p = 50e5, .t = 300.0, .z = {0.50, 0.50}, .twoPhase = true,
     .hHeos = -25332.042991, .hPr = -24989.869056,
     .liqFracHeos = 1.0 - 0.382560, .liqFracPr = 1.0 - 0.362444},
    {.name = "case2", .regime = "single-liquid",
     .p = 50e5, .t = 280.0, .z = {0.10, 0.90}, .twoPhase = false,
     .hHeos = -50638.687785, .hPr = -49118.923273,
     .liqFracHeos = 1.0, .liqFracPr = 1.0},
    // case3: single dense phase far above methane's critical temperature —
    // CoolProp's own backends disagree on the LABEL (HEOS "vapor",
    // PR "liquid"); only single-phase-ness is compared.
    {.name = "case3", .regime = "single-vapor (supercritical C1)",
     .p = 50e5, .t = 410.0, .z = {0.99, 0.01}, .twoPhase = false,
     .hHeos = 4060.067716, .hPr = 3968.274922,
     .liqFracHeos = 0.0, .liqFracPr = 0.0},
    {.name = "case4", .regime = "near-datum (two-phase)",
     .p = 1e5, .t = 298.15, .z = {0.50, 0.50}, .twoPhase = true,
     .hHeos = -25637.710968, .hPr = -24761.290139,
     .liqFracHeos = 1.0 - 0.499395, .liqFracPr = 1.0 - 0.498542},
    // case5: a single compressed-liquid state — at 200 bar the 50/50 feed
    // stays below its bubble point at 350 K on both CoolProp backends.
    {.name = "case5", .regime = "high-P single-liquid (round-trip stress)",
     .p = 200e5, .t = 350.0, .z = {0.50, 0.50}, .twoPhase = false,
     .hHeos = -15571.205136, .hPr = -15730.281124,
     .liqFracHeos = 1.0, .liqFracPr = 1.0},
}};

struct OpmResult {
    bool singlePhase;
    double L;           // liquid molar fraction from the flash
    double hCaloric;    // [J/mol]
    double hDeparture;  // [J/mol]
};

FluidState makeState(double p, double t, const std::array<double, 2>& z)
{
    FluidState fs;
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        fs.setPressure(phaseIdx, p);
    }
    fs.setTemperature(t);
    for (unsigned compIdx = 0; compIdx < FluidSystem::numComponents; ++compIdx) {
        fs.setMoleFraction(compIdx, z[compIdx]);
    }
    for (unsigned compIdx = 0; compIdx < FluidSystem::numComponents; ++compIdx) {
        fs.setKvalue(compIdx, fs.wilsonK_(compIdx));
    }
    fs.setLvalue(-1.0);
    return fs;
}

Opm::CpTable<double, 2> cpTable()
{
    return {Opm::IdealGasCaloricData<double>::methane(),
            Opm::IdealGasCaloricData<double>::decane()};
}

OpmResult flashPT(const ReferenceCase& rc)
{
    FluidState fs = makeState(rc.p, rc.t, rc.z);
    const bool singlePhase =
        PtFlash::solve(fs, TWO_PHASE_METHOD, PT_TOLERANCE, EOSType::PR);

    const auto table = cpTable();
    return {
        .singlePhase = singlePhase,
        .L = Opm::getValue(fs.L()),
        .hCaloric = Opm::getValue(Enthalpy::mixtureEnthalpy(
            fs, table, T0, EOSType::PR, Opm::EnthalpyModel::caloric)),
        .hDeparture = Opm::getValue(Enthalpy::mixtureEnthalpy(
            fs, table, T0, EOSType::PR, Opm::EnthalpyModel::eos_departure)),
    };
}

typename PhFlash::Config phConfig(Opm::EnthalpyModel model)
{
    typename PhFlash::Config cfg;
    cfg.cpTable = cpTable();
    cfg.model = model;
    // The default bracket [200, 600] K fails for the 99%-methane case: at
    // the 200 K endpoint (50 bar) that feed is near-critical (methane
    // Tc = 190.6 K) and the isothermal endpoint flash does not converge, so
    // PHFlash's bracket pre-check correctly reports no solution. Narrow to
    // the window the monotonicity probes of this test already flash
    // successfully — every case temperature (280..410 K) lies well inside.
    // (PHFlashConfig's own documentation advises narrowing the default
    // bracket for quantitative work.)
    cfg.tempMin = 270.0;
    cfg.tempMax = 460.0;
    return cfg;
}

//! PHFlash inversion; returns the recovered temperature, NaN if the solver
//! reports no solution on the bracket (displayed honestly, never fabricated)
double invertPH(const ReferenceCase& rc, double hSpec, Opm::EnthalpyModel model)
{
    FluidState fs = makeState(rc.p, rc.t, rc.z);
    const bool ok = PhFlash::solve(fs, hSpec, phConfig(model),
                                   TWO_PHASE_METHOD, PT_TOLERANCE, EOSType::PR);
    return ok ? Opm::getValue(fs.temperature(0))
              : std::numeric_limits<double>::quiet_NaN();
}

} // anonymous namespace

// ── The shown comparison ───────────────────────────────────────────────────
// Never fails on magnitudes: this is the human-facing table. Run the binary
// directly (--log_level=message not required; the table goes to stdout) to
// see OPM against both CoolProp backends, with the PR-vs-HEOS spread printed
// as the yardstick OPM cannot be expected to beat.
BOOST_AUTO_TEST_CASE(ShownComparison)
{
    std::cout << "\n=== P-H flash vs CoolProp references (J/mol, K) ===\n"
              << "    yardstick = |CoolProp/PR - CoolProp/HEOS|: the EoS model gap\n\n"
              << std::fixed << std::setprecision(1);

    for (const auto& rc : REFERENCE_CASES) {
        const OpmResult opm = flashPT(rc);
        const double yardstick = std::abs(rc.hPr - rc.hHeos);
        const double tBack = invertPH(rc, rc.hPr, Opm::EnthalpyModel::eos_departure);

        std::cout << rc.name << "  [" << rc.regime << "]  "
                  << "P = " << rc.p / 1e5 << " bar, T = " << rc.t
                  << " K, z(C1) = " << std::setprecision(2) << rc.z[0]
                  << std::setprecision(1) << "\n"
                  << "  H  opm/caloric " << std::setw(11) << opm.hCaloric
                  << "   opm/departure " << std::setw(11) << opm.hDeparture
                  << "   cp/PR " << std::setw(11) << rc.hPr
                  << "   cp/HEOS " << std::setw(11) << rc.hHeos << "\n"
                  << "  dH opm/departure - cp/PR = "
                  << std::setw(9) << opm.hDeparture - rc.hPr
                  << "   vs yardstick " << std::setw(8) << yardstick << "\n"
                  << "  L  opm " << std::setprecision(4) << opm.L
                  << "   cp/PR " << rc.liqFracPr
                  << "   cp/HEOS " << rc.liqFracHeos << std::setprecision(1) << "\n"
                  << "  T(cp/PR enthalpy) via PHFlash = " << std::setprecision(3)
                  << tBack << "  (case T = " << rc.t << ", dT = "
                  << tBack - rc.t << ")" << std::setprecision(1) << "\n\n";
    }
    BOOST_CHECK(true);  // the table itself is the deliverable of this case
}

// ── M1: regime agreement ───────────────────────────────────────────────────
BOOST_AUTO_TEST_CASE(M1_RegimeAgreement)
{
    for (const auto& rc : REFERENCE_CASES) {
        const OpmResult opm = flashPT(rc);
        BOOST_TEST_CONTEXT(rc.name << " [" << rc.regime << "]") {
            BOOST_CHECK_EQUAL(!opm.singlePhase, rc.twoPhase);
        }
    }
}

// ── M2: monotonicity of H(T) around every case point ───────────────────────
BOOST_AUTO_TEST_CASE(M2_Monotonicity)
{
    for (const auto& rc : REFERENCE_CASES) {
        // probe temperatures around the case point, clamped to the validity
        // window and deduplicated (clamping may collapse neighbors — strict
        // monotonicity needs strictly increasing probes)
        std::vector<double> probes;
        for (const double dT : {-30.0, -10.0, 0.0, +10.0, +30.0}) {
            const double t = std::clamp(rc.t + dT, 275.0, 460.0);
            if (probes.empty() || t > probes.back()) {
                probes.push_back(t);
            }
        }

        double prevCal = -std::numeric_limits<double>::infinity();
        double prevDep = -std::numeric_limits<double>::infinity();
        for (const double t : probes) {
            ReferenceCase probe = rc;
            probe.t = t;
            const OpmResult opm = flashPT(probe);
            BOOST_TEST_CONTEXT(rc.name << " at T = " << t) {
                BOOST_CHECK_GT(opm.hCaloric, prevCal);
                BOOST_CHECK_GT(opm.hDeparture, prevDep);
            }
            prevCal = opm.hCaloric;
            prevDep = opm.hDeparture;
        }
    }
}

// ── M3: rank agreement of the case enthalpies ──────────────────────────────
BOOST_AUTO_TEST_CASE(M3_RankAgreement)
{
    std::array<std::size_t, REFERENCE_CASES.size()> byOpm;
    std::array<std::size_t, REFERENCE_CASES.size()> byRef;
    std::iota(byOpm.begin(), byOpm.end(), 0);
    byRef = byOpm;

    std::array<double, REFERENCE_CASES.size()> hOpm;
    std::transform(REFERENCE_CASES.begin(), REFERENCE_CASES.end(), hOpm.begin(),
                   [](const auto& rc) { return flashPT(rc).hDeparture; });

    std::sort(byOpm.begin(), byOpm.end(),
              [&](auto a, auto b) { return hOpm[a] < hOpm[b]; });
    std::sort(byRef.begin(), byRef.end(),
              [](auto a, auto b)
              { return REFERENCE_CASES[a].hPr < REFERENCE_CASES[b].hPr; });

    BOOST_CHECK_EQUAL_COLLECTIONS(byOpm.begin(), byOpm.end(),
                                  byRef.begin(), byRef.end());
}

// ── M4: internal round-trip closure ────────────────────────────────────────
// Flash at the case temperature, take the OPM enthalpy, invert it with
// PHFlash: the case temperature must come back within the solver's own
// bracketing tolerance (1e-6 on the T interval; asserted with a 10x margin
// for bracket-end rounding). This is PR-quality internal consistency — it
// involves no cross-model comparison at all.
BOOST_AUTO_TEST_CASE(M4_RoundTripClosure)
{
    for (const auto& rc : REFERENCE_CASES) {
        for (const auto model : {Opm::EnthalpyModel::caloric,
                                 Opm::EnthalpyModel::eos_departure}) {
            const OpmResult opm = flashPT(rc);
            const double h = (model == Opm::EnthalpyModel::caloric)
                             ? opm.hCaloric : opm.hDeparture;
            const double tBack = invertPH(rc, h, model);
            BOOST_TEST_CONTEXT(rc.name << " model="
                               << Opm::enthalpyModelToString(model)) {
                BOOST_CHECK(!std::isnan(tBack));
                BOOST_CHECK_SMALL(tBack - rc.t, 1e-5);
            }
        }
    }
}
