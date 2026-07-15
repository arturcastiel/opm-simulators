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
 * \ingroup PHFlashModel
 *
 * \brief Declares the runtime parameters of the isenthalpic (P-H)
 *        compositional model.
 *
 * The defaults are copied verbatim from PHFlashConfig
 * (opm/material/constraintsolvers/PHFlash.hpp) so that building the config
 * from these parameters is a plain value copy. The inner isothermal flash
 * keeps its own parameters (FlashTolerance, FlashTwoPhaseMethod,
 * FlashVerbosity) — those govern the P-T solve that runs inside every trial
 * temperature.
 */
#ifndef OPM_PHFLASH_PARAMETERS_HH
#define OPM_PHFLASH_PARAMETERS_HH

#include <opm/models/ptflash/flashparameters.hh>

namespace Opm::Parameters {

//! Which enthalpy model the P-H flash inverts: "caloric" or "eos_departure"
struct PhFlashEnthalpyModel { static constexpr auto value = "caloric"; };

//! [K] lower bound of the temperature search bracket. Both bracket ends
//! are full inner-flash evaluations, so the defaults are the
//! field-validated window inside the isothermal flash's robust envelope
//! and the cp-correlation validity range — near-critical (methane-rich)
//! feeds fail to flash at colder bounds, killing cells with a misleading
//! "no solution" although the root is interior. Widen deliberately.
template<class Scalar>
struct PhFlashTempMin { static constexpr Scalar value = 270.0; };

//! [K] upper bound of the temperature search bracket
template<class Scalar>
struct PhFlashTempMax { static constexpr Scalar value = 460.0; };

//! convergence tolerance of the outer bracketing solver (applied to both the
//! temperature interval [K] and the enthalpy residual [J/mol])
template<class Scalar>
struct PhFlashTolerance { static constexpr Scalar value = 1e-6; };

//! maximum iterations of the outer bracketing solver
struct PhFlashMaxIterations { static constexpr int value = 100; };

} // namespace Opm::Parameters

#endif
