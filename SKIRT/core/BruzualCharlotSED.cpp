/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BruzualCharlotSED.hpp"
#include "BruzualCharlotSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

void BruzualCharlotSED::setupSelfBefore()
{
    StellarSED::setupSelfBefore();

    // construct the library of SED models
    BruzualCharlotSEDFamily bc(this);

    // get the luminosities for arbitrary mass and for the appropriate Z and t
    setLuminosities(bc.luminosities(1., _metallicity, _age/Constants::year()));
}

////////////////////////////////////////////////////////////////////
