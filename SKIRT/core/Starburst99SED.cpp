/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Starburst99SED.hpp"
#include "Starburst99SEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

void Starburst99SED::setupSelfBefore()
{
    StellarSED::setupSelfBefore();

    // construct the library of SED models
    Starburst99SEDFamily sb(this);

    // get the luminosities for arbitrary mass and for the appropriate Z and t
    setLuminosities(sb.luminosities(1., _metallicity, _age/Constants::year()));
}

////////////////////////////////////////////////////////////////////
