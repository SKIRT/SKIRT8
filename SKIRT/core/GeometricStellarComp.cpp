/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GeometricStellarComp.hpp"
#include "FatalError.hpp"
#include "PhotonPackage.hpp"
#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

int GeometricStellarComp::dimension() const
{
    return _geometry->dimension();
}

//////////////////////////////////////////////////////////////////////

void GeometricStellarComp::launch(PhotonPackage* pp, int ell, double L) const
{
    Position bfr = _geometry->generatePosition();
    Direction bfk = _geometry->generateDirection(ell,bfr);
    pp->launch(L,ell,bfr,bfk);
    pp->setAngularDistribution(_geometry);
}

//////////////////////////////////////////////////////////////////////
