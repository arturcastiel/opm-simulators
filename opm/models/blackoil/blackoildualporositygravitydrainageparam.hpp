/*
  Copyright 2026 TNO

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
#ifndef OPM_BLACKOIL_DUAL_POROSITY_GRAVITY_DRAINAGE_PARAM_HPP
#define OPM_BLACKOIL_DUAL_POROSITY_GRAVITY_DRAINAGE_PARAM_HPP

#include <vector>

namespace Opm {

/*!
 * Per-cell static data for the dual-porosity gravity-drainage flux terms:
 * the mobile-fraction end points, the initial saturations and fractions,
 * and the typical matrix block height.  Filled once by the problem when a
 * gravity-drainage model is requested; empty otherwise.  CPU-side only —
 * dual-porosity runs are serial.
 */
template<class Scalar>
struct DualPorosityGravityDrainageParam
{
    struct CellData
    {
        // water-fraction end points
        Scalar swco = 0.0;
        Scalar swcr = 0.0;
        Scalar scohy = 0.0;
        Scalar scrhy = 0.0;
        // gas-fraction end points
        Scalar sgco = 0.0;
        Scalar sgcr = 0.0;
        Scalar slco = 0.0;
        Scalar slcr = 0.0;
        // initial saturations and mobile fractions
        Scalar swi = 0.0;
        Scalar sgi = 0.0;
        Scalar xwi = 0.0;
        Scalar xgi = 0.0;
        // typical matrix block height (zero = the keyword's zero-effect default)
        Scalar dzMatrix = 0.0;
    };

    bool active = false;
    Scalar gravity = 0.0;
    std::vector<CellData> cell;
};

} // namespace Opm

#endif // OPM_BLACKOIL_DUAL_POROSITY_GRAVITY_DRAINAGE_PARAM_HPP
