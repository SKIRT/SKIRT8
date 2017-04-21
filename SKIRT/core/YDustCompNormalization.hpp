/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef YDUSTCOMPNORMALIZATION_HPP
#define YDUSTCOMPNORMALIZATION_HPP

#include "DustCompNormalization.hpp"

//////////////////////////////////////////////////////////////////////

/** YDustCompNormalization is a class that sets the normalization of a general 3D dust component by
    defining the Y-axis optical depth at an arbitrary wavelength. The Y-axis optical depth is
    defined as the integral of the opacity along the entire Y-axis, \f[
    \tau_\lambda^{\text{Y}} = \int_{-\infty}^\infty k_\lambda(0,y,0)\, {\text{d}}y. \f] */
class YDustCompNormalization : public DustCompNormalization
{
    ITEM_CONCRETE(YDustCompNormalization, DustCompNormalization,
                  "normalization by defining the Y-axis optical depth at some wavelength")

    PROPERTY_DOUBLE(wavelength, "the wavelength at which to specify the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

    PROPERTY_DOUBLE(opticalDepth, "the Y-axis optical depth at this wavelength")
        ATTRIBUTE_MIN_VALUE(opticalDepth, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the appropriate normalization factor for the specified
        geometry and dust mixture. */
    double normalizationFactor(const Geometry* geom, const DustMix* mix) const override;
};

//////////////////////////////////////////////////////////////////////

#endif
