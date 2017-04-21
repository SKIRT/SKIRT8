/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPECTRALLUMINOSITYSTELLARCOMPNORMALIZATION_HPP
#define SPECTRALLUMINOSITYSTELLARCOMPNORMALIZATION_HPP

#include "StellarCompNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** SpectralLuminosityStellarCompNormalization is a class that sets the normalization of a stellar
    component by defining the spectral luminosity (radiative power per wavelength) at a certain
    wavelength. */
class SpectralLuminosityStellarCompNormalization : public StellarCompNormalization
{
    ITEM_CONCRETE(SpectralLuminosityStellarCompNormalization, StellarCompNormalization,
                  "stellar component normalization through the luminosity at a given wavelength")

    PROPERTY_DOUBLE(wavelength, "the wavelength at which to specify the luminosity")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

    PROPERTY_DOUBLE(luminosity, "this component's luminosity at the specified wavelength")
        ATTRIBUTE_QUANTITY(luminosity, "monluminosity")
        ATTRIBUTE_MIN_VALUE(luminosity, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies the validity of the luminosity value and caches the wavelength bin
        corresponding to the specified wavelength. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the total, bolometric luminosity of a (virtual) stellar component
        that would have a given %SED. For the present type of normalization, the bolometric
        luminosity is \f[ L_{\text{bol}} = \frac{L_{\lambda,\ell}\Delta\lambda_\ell}{S_\ell} \f]
        with \f$\ell\f$ the wavelength bin corresponding to the specified wavelength,
        \f$\Delta\lambda_\ell\f$ the width of that wavelength bin, \f$L_{\lambda,\ell}\f$ the
        specified spectral luminosity (for the wavelength corresponding to the bin), and
        \f$S_\ell\f$ the value of the (normalized) %SED for the same wavelength bin. */
    double totalLuminosity(SED* sed) const override;

    //======================== Data Members ========================

private:
    // initialized during setup
    int _ell{0};
    double _dlambda{0.};
};

////////////////////////////////////////////////////////////////////

#endif
