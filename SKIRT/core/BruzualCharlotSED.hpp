/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BRUZUALCHARLOTSED_HPP
#define BRUZUALCHARLOTSED_HPP

#include "StellarSED.hpp"

////////////////////////////////////////////////////////////////////

/** BruzualCharlotSED is a class that represents spectral energy distributions of simple stellar
    populations (SSPs), parameterized on metallicity and age according to the model of Bruzual &
    Charlot (2003). See the BruzualCharlotSEDFamily class for more information. */
class BruzualCharlotSED : public StellarSED
{
    ITEM_CONCRETE(BruzualCharlotSED, StellarSED, "a Bruzual & Charlot simple stellar population SED")

    PROPERTY_DOUBLE(metallicity, "the metallicity of the SSP")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.0001")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.05]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

    PROPERTY_DOUBLE(age, "the age of the SSP")
        ATTRIBUTE_QUANTITY(age, "time")
        ATTRIBUTE_MIN_VALUE(age, "[0 Gyr")
        ATTRIBUTE_MAX_VALUE(age, "20 Gyr]")
        ATTRIBUTE_DEFAULT_VALUE(age, "5 Gyr")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a temporary instance of the BruzualCharlotSEDFamily class to
        obtain an SED that corresponds to the values of the metallicity and age attributes. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
