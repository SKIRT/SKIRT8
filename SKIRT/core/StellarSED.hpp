/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STELLARSED_HPP
#define STELLARSED_HPP

#include "SED.hpp"

////////////////////////////////////////////////////////////////////

/** StellarSED is an abstract subclass of the SED class and represents spectral energy
    distributions corresponding to stellar sources. */
class StellarSED : public SED
{
    ITEM_ABSTRACT(StellarSED, SED, "a stellar SED")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
