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
 * \ingroup FlashModel
 *
 * \brief Declares the properties required by the compositional
 *        multi-phase model based on flash calculations.
 */
#ifndef EWOMS_FLASH_PROPERTIES_HH
#define EWOMS_FLASH_PROPERTIES_HH

#include <opm/models/utils/propertysystem.hh>

namespace Opm::Properties {

/*!
 * \brief The type of the flash constraint solver.
 *
 * NOTE: this property is duck-typed — its call contract is defined by the
 * consuming IntensiveQuantities, and two distinct contracts already exist:
 * - the NCP-based flash model calls
 *   FlashSolver::guessInitial(fluidState, globalMolarities) and
 *   FlashSolver::solve<MaterialLaw>(fluidState, matParams, paramCache,
 *   globalMolarities, tolerance)
 *   (see opm/models/flash/flashintensivequantities.hh);
 * - the PT flash model calls
 *   FlashSolver::solve(fluidState, twoPhaseMethod, tolerance, eosType,
 *   verbosity)
 *   (see opm/models/ptflash/flashintensivequantities.hh).
 * A binding must match the contract of the model consuming it; do not add a
 * third meaning to this property.
 */
template<class TypeTag, class MyTypeTag>
struct FlashSolver { using type = UndefinedProperty; };

/*!
 * \brief The type of the isenthalpic (P-H) flash solver used by the PH flash
 *        model.
 *
 * Deliberately a separate property from FlashSolver: the P-H solver's
 * contract differs (the specified quantity is enthalpy, not temperature).
 * Expected contract: a nested Config type and
 * solve(fluidState, hSpec, config, twoPhaseMethod, ptTolerance, eosType,
 * verbosity) returning bool (false = no solution on the bracket).
 */
template<class TypeTag, class MyTypeTag>
struct PhFlashSolver { using type = UndefinedProperty; };

} // namespace Opm::Properties

#endif
