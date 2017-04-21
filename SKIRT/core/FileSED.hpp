/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILESED_HPP
#define FILESED_HPP

#include "StellarSED.hpp"

////////////////////////////////////////////////////////////////////

/** FileSED is a simple class that represents spectral energy distributions read in directly from a
    file provided by the user. */
class FileSED : public StellarSED
{
    ITEM_CONCRETE(FileSED, StellarSED, "a tabulated SED from a file")

    PROPERTY_STRING(filename, "the name of the file that contains the SED")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the stellar fluxes from a file provided by the user. This file should
        first contain a single line with the number of data points and subsequently contain lines
        with two columns: wavelength \f$\lambda\f$ in micron and flux density \f$F_\lambda\f$ in
        arbitrary units. This vector is regridded on the global wavelength grid. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
