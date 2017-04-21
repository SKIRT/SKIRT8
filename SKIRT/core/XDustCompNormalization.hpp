/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XDUSTCOMPNORMALIZATION_HPP
#define XDUSTCOMPNORMALIZATION_HPP

#include "DustCompNormalization.hpp"

//////////////////////////////////////////////////////////////////////

/** XDustCompNormalization is a class that sets the normalization of a general 3D dust component by
    defining the X-axis optical depth at an arbitrary wavelength. The X-axis optical depth is
    defined as the integral of the opacity along the entire X-axis, \f[
    \tau_\lambda^{\text{X}} = \int_{-\infty}^\infty k_\lambda(x,0,0)\, {\text{d}}x. \f] */
class XDustCompNormalization : public DustCompNormalization
{
    ITEM_CONCRETE(XDustCompNormalization, DustCompNormalization,
                  "normalization by defining the X-axis optical depth at some wavelength")

    PROPERTY_DOUBLE(wavelength, "the wavelength at which to specify the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

    PROPERTY_DOUBLE(opticalDepth, "the X-axis optical depth at this wavelength")
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
