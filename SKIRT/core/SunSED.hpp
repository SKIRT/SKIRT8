/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SUNSED_HPP
#define SUNSED_HPP

#include "StellarSED.hpp"

////////////////////////////////////////////////////////////////////

/** The SunSED class represents the spectral energy distribution from the Sun. */
class SunSED : public StellarSED
{
    ITEM_CONCRETE(SunSED, StellarSED, "the SED of the Sun")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the solar fluxes from a resource file into a vector, which is
        then resampled on the global wavelength grid. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
