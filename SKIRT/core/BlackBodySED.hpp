/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BLACKBODYSED_HPP
#define BLACKBODYSED_HPP

#include "StellarSED.hpp"

////////////////////////////////////////////////////////////////////

/** BlackBodySED is a class that describes black-body spectral energy distributions, i.e. the
    emission spectra of perfect absorbers which are in thermal equilibrium. Such an %SED is
    characterized by the temperature of the object, and its spectrum is the Planck spectrum. */
class BlackBodySED : public StellarSED
{
    ITEM_CONCRETE(BlackBodySED, StellarSED, "a black body SED")

    PROPERTY_DOUBLE(temperature, "the black body temperature")
        ATTRIBUTE_QUANTITY(temperature, "temperature")
        ATTRIBUTE_MIN_VALUE(temperature, "]0 K")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates a vector with the Planck function \f$ B_\lambda(T) \f$ sampled at
        all grid points of the global wavelength grid. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
