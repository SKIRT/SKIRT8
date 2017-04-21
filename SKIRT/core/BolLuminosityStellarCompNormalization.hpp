/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOLLUMINOSITYSTELLARCOMPNORMALIZATION_HPP
#define BOLLUMINOSITYSTELLARCOMPNORMALIZATION_HPP

#include "StellarCompNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** BolLuminosityStellarCompNormalization is a class that sets the normalization of a stellar
    component by defining the total bolometric luminosity. */
class BolLuminosityStellarCompNormalization : public StellarCompNormalization
{
    ITEM_CONCRETE(BolLuminosityStellarCompNormalization, StellarCompNormalization,
                  "stellar component normalization through the bolometric luminosity")

    PROPERTY_DOUBLE(luminosity, "the bolometric luminosity for this component")
        ATTRIBUTE_QUANTITY(luminosity, "bolluminosity")
        ATTRIBUTE_MIN_VALUE(luminosity, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the bolometric luminosity of a (virtual) stellar component that
        would have a given %SED. For the present type of normalization, this function is trivial as
        the bolometric luminosity is a data member. */
    double totalLuminosity(SED* sed) const override;
};

////////////////////////////////////////////////////////////////////

#endif
