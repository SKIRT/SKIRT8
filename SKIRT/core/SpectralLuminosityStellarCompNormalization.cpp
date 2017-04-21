/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpectralLuminosityStellarCompNormalization.hpp"
#include "FatalError.hpp"
#include "SED.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void SpectralLuminosityStellarCompNormalization::setupSelfBefore()
{
    StellarCompNormalization::setupSelfBefore();

    // remember the wavelength index and bin width corresponding to the specified wavelength
    WavelengthGrid* lambdagrid = find<WavelengthGrid>();
    _ell = lambdagrid->nearest(_wavelength);
    if (_ell < 0) throw FATALERROR("The given wavelength is outside of the simulation's wavelength grid");
    _dlambda = lambdagrid->dlambda(_ell);
}

////////////////////////////////////////////////////////////////////

double SpectralLuminosityStellarCompNormalization::totalLuminosity(SED* sed) const
{
    return _luminosity * _dlambda / sed->luminosity(_ell);
}

//////////////////////////////////////////////////////////////////////
