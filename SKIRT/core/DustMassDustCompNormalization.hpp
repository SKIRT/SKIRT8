/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTMASSDUSTCOMPNORMALIZATION_HPP
#define DUSTMASSDUSTCOMPNORMALIZATION_HPP

#include "DustCompNormalization.hpp"

//////////////////////////////////////////////////////////////////////

/** DustMassDustCompNormalization is a class that sets the normalization of a dust component by
    defining the total dust mass. */
class DustMassDustCompNormalization : public DustCompNormalization
{
    ITEM_CONCRETE(DustMassDustCompNormalization, DustCompNormalization, "normalization by defining the total dust mass")

    PROPERTY_DOUBLE(dustMass, "the total dust mass of the dust component")
        ATTRIBUTE_QUANTITY(dustMass, "mass")
        ATTRIBUTE_MIN_VALUE(dustMass, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the appropriate normalization factor for the specified
        geometry and dust mixture. */
    double normalizationFactor(const Geometry* geom, const DustMix* mix) const override;
};

//////////////////////////////////////////////////////////////////////

#endif
