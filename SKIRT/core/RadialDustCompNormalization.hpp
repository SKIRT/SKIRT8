/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RADIALDUSTCOMPNORMALIZATION_HPP
#define RADIALDUSTCOMPNORMALIZATION_HPP

#include "DustCompNormalization.hpp"

//////////////////////////////////////////////////////////////////////

/** RadialDustCompNormalization is a class that sets the normalization of a spherically
    symmetric dust component by defining the radial optical depth at an arbitrary wavelength.
    The radial optical depth is defined as the integral of the opacity along the line from the
    centre to the edge, \f[
    \tau_\lambda^{\text{rad}} = \int_{0}^\infty k_\lambda(r)\, {\text{d}}r. \f] */
class RadialDustCompNormalization : public DustCompNormalization
{
    ITEM_CONCRETE(RadialDustCompNormalization, DustCompNormalization,
                  "normalization by defining the radial optical depth at some wavelength")

    PROPERTY_DOUBLE(wavelength, "the wavelength at which to specify the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

    PROPERTY_DOUBLE(opticalDepth, "the radial optical depth at this wavelength")
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
