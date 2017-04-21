/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OligoWavelengthGrid.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void OligoWavelengthGrid::setupSelfBefore()
{
    WavelengthGrid::setupSelfBefore();

    // copy the wavelengths and sort them in ascending order
    Array lambdav;
    NR::assign(lambdav, _wavelengths);
    NR::sort(lambdav);

    // set the wavelength grid using artifical bin widths
    setWavelengthBins(lambdav, 0.001*lambdav);
}

//////////////////////////////////////////////////////////////////////

bool OligoWavelengthGrid::isSampledRange() const
{
    return false;
}

//////////////////////////////////////////////////////////////////////
