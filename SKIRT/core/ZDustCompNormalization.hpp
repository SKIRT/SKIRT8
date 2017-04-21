/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ZDUSTCOMPNORMALIZATION_HPP
#define ZDUSTCOMPNORMALIZATION_HPP

#include "DustCompNormalization.hpp"

//////////////////////////////////////////////////////////////////////

/** ZDustCompNormalization is a class that sets the normalization of a general 3D dust component by
    defining the Z-axis optical depth at an arbitrary wavelength. The Z-axis optical depth is
    defined as the integral of the opacity along the entire Z-axis, \f[
    \tau_\lambda^{\text{Z}} = \int_{-\infty}^\infty k_\lambda(0,0,z)\, {\text{d}}z. \f] */
class ZDustCompNormalization : public DustCompNormalization
{
    ITEM_CONCRETE(ZDustCompNormalization, DustCompNormalization,
                  "normalization by defining the Z-axis optical depth at some wavelength")

    PROPERTY_DOUBLE(wavelength, "the wavelength at which to specify the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

    PROPERTY_DOUBLE(opticalDepth, "the Z-axis optical depth at this wavelength")
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
