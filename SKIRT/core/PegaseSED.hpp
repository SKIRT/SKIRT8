/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PEGASESED_HPP
#define PEGASESED_HPP

#include "StellarSED.hpp"

////////////////////////////////////////////////////////////////////

/** PegaseSED is a class that represents template galaxy %SEDs from the Pegase library. The library
    offers templates for elliptical, lenticular and spiral galaxies. See Fioc & Rocca-Volmerange
    (1997, A&A, 326, 950–962). */
class PegaseSED : public StellarSED
{
    /** The enumeration type indicating the spectral type of the %SED. */
    ENUM_DEF(SpectralType, E, S0, Sa, Sb, Sc)
    ENUM_VAL(SpectralType, E, "elliptical galaxy (E)")
    ENUM_VAL(SpectralType, S0, "lenticular galaxy (S0)")
    ENUM_VAL(SpectralType, Sa, "early-type spiral galaxy (Sa)")
    ENUM_VAL(SpectralType, Sb, "intermediate-type spiral galaxy (Sb)")
    ENUM_VAL(SpectralType, Sc, "late-type spiral galaxy (Sc)")
    ENUM_END()

    ITEM_CONCRETE(PegaseSED, StellarSED, "a Pegase galaxy SED")

    PROPERTY_ENUM(spectralType, SpectralType, "the spectral type of the SED")
        ATTRIBUTE_DEFAULT_VALUE(spectralType, "E")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads fluxes from a resource file corresponding to the type of the desired
        Pegase template. The flux vector is regridded on the global wavelength grid. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
