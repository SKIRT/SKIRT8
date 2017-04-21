/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef KURUCZSED_HPP
#define KURUCZSED_HPP

#include "StellarSED.hpp"

////////////////////////////////////////////////////////////////////

/** KuruczSED is a class that represents spectral energy distributions of stars according to the
    model of Kurucz (1993). SSPs with different metallicities, effective temperatures and surface
    gravities can be chosen. */
class KuruczSED : public StellarSED
{
    ITEM_CONCRETE(KuruczSED, StellarSED, "a Kurucz SED")

    PROPERTY_DOUBLE(metallicity, "the metallicity")
        ATTRIBUTE_MIN_VALUE(metallicity, "[-5")
        ATTRIBUTE_MAX_VALUE(metallicity, "1]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

    PROPERTY_DOUBLE(temperature, "the effective temperature")
        ATTRIBUTE_QUANTITY(temperature, "temperature")
        ATTRIBUTE_MIN_VALUE(temperature, "[3500 K")
        ATTRIBUTE_MAX_VALUE(temperature, "10000 K]")
        ATTRIBUTE_DEFAULT_VALUE(temperature, "6000 K")

    PROPERTY_DOUBLE(gravity, "the surface gravity")
        ATTRIBUTE_MIN_VALUE(gravity, "[0")
        ATTRIBUTE_MAX_VALUE(gravity, "5]")
        ATTRIBUTE_DEFAULT_VALUE(gravity, "2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** The SKIRT resource data includes a huge library with Kurucz model atmospheres on a 3D
        parameter space grid in metallicity, effective temperature and surface gravity. This
        function first determines the metallicity from the metallicity grid that is closest to the
        desired metallicity \f$\mu\f$. The same goes for the surface gravity \f$g\f$. Once these
        values are known, the correct file is searched and the constructor calculates a vector with
        the %SED by interpolating between the two SEDs with effective temperatures that bracket the
        desired effective temperature \f$T_{\text{eff}}\f$. This vector is then regridded on the
        global wavelength grid. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
