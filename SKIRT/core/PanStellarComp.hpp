/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PANSTELLARCOMP_HPP
#define PANSTELLARCOMP_HPP

#include "GeometricStellarComp.hpp"
#include "StellarCompNormalization.hpp"
#include "StellarSED.hpp"

//////////////////////////////////////////////////////////////////////

/** The PanStellarComp class represents a stellar component that uses a built-in geometry in a
    panchromatic simulation. It uses a spectral energy distribution (an instance of the StellarSED
    class) and a normalization method (an instance of the StellarCompNormalization class) to
    calculate the luminosity vector maintained by the GeometricStellarComp base class. */
class PanStellarComp : public GeometricStellarComp
{
    ITEM_CONCRETE(PanStellarComp, GeometricStellarComp,
                  "a stellar component with a built-in geometry (in an panchromatic simulation)")
        ATTRIBUTE_ALLOWED_IF(PanStellarComp, "PanMonteCarloSimulation")

    PROPERTY_ITEM(sed, StellarSED, "the spectral energy distribution for the stellar component")

    PROPERTY_ITEM(normalization, StellarCompNormalization, "the type of normalization for the stellar component")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that this item resides in a panchromatic simulation. */
    void setupSelfBefore() override;

    /** This function calculates the luminosity vector maintained by the GeometricStellarComp base
        class using the spectral energy distribution and normalization method provided as
        attributes. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the luminosity \f$L_\ell\f$ of the stellar component at the
        wavelength index \f$\ell\f$. It is found as \f$L_\ell = L_{\text{tot}}\, S_\ell\f$, with
        \f$L_{\text{tot}}\f$ the total bolometric luminosity and \f$S_\ell\f$ the value of the SED
        at the wavelength index \f$\ell\f$. */
    double luminosity(int ell) const override;

    //======================== Data Members ========================

private:
    Array _Lv;
};

//////////////////////////////////////////////////////////////////////

#endif
