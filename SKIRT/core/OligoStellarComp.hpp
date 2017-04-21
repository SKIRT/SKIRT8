/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OLIGOSTELLARCOMP_HPP
#define OLIGOSTELLARCOMP_HPP

#include "GeometricStellarComp.hpp"

//////////////////////////////////////////////////////////////////////

/** The OligoStellarComp class represents a stellar component that uses a built-in geometry in an
    oligochromatic simulation. The spectral energy distribution over the small number of
    wavelengths is managed internally rather than through an instance of the StellarSED class. */
class OligoStellarComp : public GeometricStellarComp
{
    ITEM_CONCRETE(OligoStellarComp, GeometricStellarComp,
                  "a stellar component with a built-in geometry (in an oligochromatic simulation)")
        ATTRIBUTE_ALLOWED_IF(OligoStellarComp, "OligoMonteCarloSimulation")

    PROPERTY_DOUBLE_LIST(luminosities, "the luminosities, one for each wavelength")
        ATTRIBUTE_QUANTITY(luminosities, "monluminosity")
        ATTRIBUTE_MIN_VALUE(luminosities, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that this item resides in an oligochromatic simulation, and that the
        number of luminosities specified by the user matches the number of wavelengths. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the luminosity \f$L_\ell\f$ of the stellar component at the
        wavelength index \f$\ell\f$. Specifically, it returns the appropriate element from the \em
        luminosities property multiplied by the corresponding wavelength bin width. */
    double luminosity(int ell) const override;

    //======================== Data Members ========================

private:
    Array _Lv;
};

//////////////////////////////////////////////////////////////////////

#endif
