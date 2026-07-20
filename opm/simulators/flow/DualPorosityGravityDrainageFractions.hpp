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
#ifndef OPM_DUAL_POROSITY_GRAVITY_DRAINAGE_FRACTIONS_HPP
#define OPM_DUAL_POROSITY_GRAVITY_DRAINAGE_FRACTIONS_HPP

namespace Opm::DualPorosityFractions {

/*!
 * The gravity-drainage flux terms of the dual-porosity models are driven by
 * the fractional heights of mobile water and mobile gas in a cell.  These
 * fractions are estimated from the cell saturation with a vertical
 * equilibrium picture: the initial fraction locates the initial contact
 * inside the cell, and the current fraction moves between the connate and
 * critical end points on either side of it.  All functions are pure and
 * clamp their result to [0, 1].
 */

struct WaterFractionEndPoints
{
    double swco;    //!< connate water saturation
    double swcr;    //!< critical water saturation
    double scohy;   //!< connate hydrocarbon saturation
    double scrhy;   //!< critical hydrocarbon saturation
};

struct GasFractionEndPoints
{
    double sgco;    //!< connate gas saturation
    double sgcr;    //!< critical gas saturation
    double slco;    //!< connate liquid saturation
    double slcr;    //!< critical liquid saturation
};

//! The current-fraction functions are templated on the saturation value
//! type so automatic-differentiation types flow through into the flux
//! terms; the initial fractions and end points are plain state.
template <class Value>
constexpr Value clampFraction(const Value& x) noexcept
{
    if (x < 0.0) {
        return Value(0.0);
    }
    if (x > 1.0) {
        return Value(1.0);
    }
    return x;
}

//! Fraction of the cell initially below the water contact.
constexpr double initialWaterFraction(const double swi,
                                      const WaterFractionEndPoints& ep) noexcept
{
    const double denom = 1.0 - ep.scohy - ep.swco;
    if (denom <= 0.0) {
        return 0.0;
    }
    return clampFraction((swi - ep.swco) / denom);
}

//! Current fraction of the cell containing mobile water.
template <class Value>
constexpr Value waterFraction(const Value& sw,
                              const double swi,
                              const double xwi,
                              const WaterFractionEndPoints& ep) noexcept
{
    const double denom = (sw >= swi)
        ? 1.0 - ep.scrhy - ep.swco
        : 1.0 - ep.scohy - ep.swcr;
    if (denom <= 0.0) {
        return Value(xwi);
    }

    const double offset = (sw >= swi)
        ? xwi * (ep.scrhy - ep.scohy)
        : xwi * (ep.swcr - ep.swco);

    return clampFraction(Value((sw - offset - ep.swco) / denom));
}

//! Fraction of the cell initially above the gas contact.
constexpr double initialGasFraction(const double sgi,
                                    const GasFractionEndPoints& ep) noexcept
{
    const double denom = 1.0 - ep.slco - ep.sgco;
    if (denom <= 0.0) {
        return 0.0;
    }
    return clampFraction((sgi - ep.sgco) / denom);
}

//! Current fraction of the cell containing mobile gas.
template <class Value>
constexpr Value gasFraction(const Value& sg,
                            const double sgi,
                            const double xgi,
                            const GasFractionEndPoints& ep) noexcept
{
    const double denom = (sg >= sgi)
        ? 1.0 - ep.slcr - ep.sgco
        : 1.0 - ep.slco - ep.sgcr;
    if (denom <= 0.0) {
        return Value(xgi);
    }

    const double offset = (sg >= sgi)
        ? xgi * (ep.slcr - ep.slco)
        : xgi * (ep.sgcr - ep.sgco);

    return clampFraction(Value((sg - offset - ep.sgco) / denom));
}

} // namespace Opm::DualPorosityFractions

#endif // OPM_DUAL_POROSITY_GRAVITY_DRAINAGE_FRACTIONS_HPP
