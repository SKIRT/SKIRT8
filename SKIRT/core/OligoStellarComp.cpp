/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OligoStellarComp.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "OligoWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void OligoStellarComp::setupSelfBefore()
{
    GeometricStellarComp::setupSelfBefore();

    // Verify that the wavelength grid (and thus the simulation) is of the Oligo type
    OligoWavelengthGrid* lambdagrid = find<OligoWavelengthGrid>();

    // Verify that the number of luminosities equals the number of wavelengths
    if (_luminosities.size() != static_cast<size_t>(lambdagrid->numWavelengths()))
        throw FATALERROR("The number of stellar component luminosities differs from the number of wavelengths");

    // Convert spectral luminosities (W/m) to luminosities in each wavelength bin (W)
    NR::assign(_Lv, _luminosities);
    _Lv *= lambdagrid->dlambdav();
}

////////////////////////////////////////////////////////////////////

double OligoStellarComp::luminosity(int ell) const
{
    return _Lv[ell];
}

//////////////////////////////////////////////////////////////////////
