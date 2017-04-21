/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PanStellarComp.hpp"
#include "PanWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void PanStellarComp::setupSelfBefore()
{
    GeometricStellarComp::setupSelfBefore();

    // verify that the wavelength grid (and thus the simulation) is of the Pan type
    find<PanWavelengthGrid>();
}

////////////////////////////////////////////////////////////////////

void PanStellarComp::setupSelfAfter()
{
    GeometricStellarComp::setupSelfAfter();

    // calculate the luminosities (we need our children to be setup for this)
    _Lv = _normalization->totalLuminosity(_sed) * _sed->luminosities();
}

////////////////////////////////////////////////////////////////////

double PanStellarComp::luminosity(int ell) const
{
    return _Lv[ell];
}

////////////////////////////////////////////////////////////////////
