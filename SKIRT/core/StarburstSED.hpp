/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STARBURSTSED_HPP
#define STARBURSTSED_HPP

#include "StellarSED.hpp"

////////////////////////////////////////////////////////////////////

/** StarburstSED is a class that represents spectral energy distributions of starbursting stellar
    populations. The SEDs are generated from data of the Starburst99 library; see Leitherer et al.
    (1999, ApJS, 123, 3–40). More specifically, they represent stellar populations with a constant,
    continuous star formation rate that have evolved for 100 Myr. The initial mass function has is
    a simple Salpeter IMF (a power law with \f$\alpha=2.35\f$) with \f$1~M_\odot\f$ and
    \f$100~M_\odot\f$ as lower and upper masses. Populations with different metallicities can be
    chosen (with \f$Z\f$ ranging between 0.001 and 0.040). */
class StarburstSED : public StellarSED
{
    ITEM_CONCRETE(StarburstSED, StellarSED, "a Starburst stellar population SED")

    PROPERTY_DOUBLE(metallicity, "the metallicity of the population")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.001")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.040]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads fluxes from a resource file. It calculates a vector with the
        SED by interpolating between the two populations with metallicities that bracket the
        desired metallicity \f$Z\f$. This vector is regridded on the global wavelength grid. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
