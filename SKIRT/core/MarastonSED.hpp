/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MARASTONSED_HPP
#define MARASTONSED_HPP

#include "StellarSED.hpp"

////////////////////////////////////////////////////////////////////

/** MarastonSED is a class that represents spectral energy distributions of simple stellar
    populations (SSPs) according to the model of Maraston (1998, MNRAS, 300, 872–892). SSPs with
    different ages and metallicities can be chosen. */
class MarastonSED : public StellarSED
{
    ITEM_CONCRETE(MarastonSED, StellarSED, "a Maraston simple stellar population SED")

    PROPERTY_DOUBLE(metallicity, "the metallicity of the SSP")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.0001")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.07]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

    PROPERTY_DOUBLE(age, "the age of the SSP")
        ATTRIBUTE_QUANTITY(age, "time")
        ATTRIBUTE_MIN_VALUE(age, "[1e-6 Gyr")
        ATTRIBUTE_MAX_VALUE(age, "15 Gyr]")
        ATTRIBUTE_DEFAULT_VALUE(age, "5 Gyr")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads fluxes from a set of resource files and calculates a vector with the
        %SED by interpolating between the four SSPs with ages and metallicities that bracket the
        desired age \f$\tau\f$ and the desired metallicity \f$Z\f$. This vector is regridded on the
        global wavelength grid. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
