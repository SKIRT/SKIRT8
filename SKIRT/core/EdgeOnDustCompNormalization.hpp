/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EDGEONDUSTCOMPNORMALIZATION_HPP
#define EDGEONDUSTCOMPNORMALIZATION_HPP

#include "DustCompNormalization.hpp"

//////////////////////////////////////////////////////////////////////

/** EdgeOnDustCompNormalization is a class that sets the normalization of an axisymmetric
    dust component by defining the edge-on optical depth at an arbitrary wavelength. The edge-on
    optical depth is defined as the integral of the opacity along a radial line in the equatorial
    plane \f$z=0\f$, \f[ \tau_\lambda^{\text{edge-on}} = \int_0^\infty k_\lambda(R,0)\,
    {\text{d}}R. \f] */
class EdgeOnDustCompNormalization : public DustCompNormalization
{
    ITEM_CONCRETE(EdgeOnDustCompNormalization, DustCompNormalization,
                  "normalization by defining the edge-on optical depth at some wavelength")

    PROPERTY_DOUBLE(wavelength, "the wavelength at which to specify the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

    PROPERTY_DOUBLE(opticalDepth, "the edge-on optical depth at this wavelength")
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
