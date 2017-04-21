/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MappingsSED.hpp"
#include "MappingsSEDFamily.hpp"

////////////////////////////////////////////////////////////////////

void MappingsSED::setupSelfBefore()
{
    StellarSED::setupSelfBefore();

    // construct the library of SED models
    MappingsSEDFamily mappings(this);

    // get the luminosities for arbitrary SFR and for the appropriate Z, logC, pressure and fPDR
    setLuminosities(mappings.luminosities(1., _metallicity, _compactness, _pressure, _coveringFactor));
}

////////////////////////////////////////////////////////////////////
