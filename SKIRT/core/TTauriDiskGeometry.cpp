/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TTauriDiskGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void TTauriDiskGeometry::setupSelfBefore()
{
    AxGeometry::setupSelfBefore();

    // verify property values
    if (_Rout <= _Rinn) throw FATALERROR("the outer radius of the disk must be larger than the inner radius");

    // calculate cached values
    _rho0 = 17.0/32.0/M_PI / (_Rd*_Rd*_zd) /
            (pow(_Rout/_Rd,17.0/8.0) - pow(_Rinn/_Rd,17.0/8.0));
}

//////////////////////////////////////////////////////////////////////

double TTauriDiskGeometry::density(double R, double z) const
{
    if (R<_Rinn || R>_Rout) return 0.0;
    double h = _zd*pow(R/_Rd,1.125);
    return _rho0 / (R/_Rd) * exp(-M_PI/4.0*pow(z/h,2));
}

//////////////////////////////////////////////////////////////////////

Position TTauriDiskGeometry::generatePosition() const
{
    double phi = 2.0 * M_PI * random()->uniform();
    double tinn = pow(_Rinn,2.125);
    double tout = pow(_Rout,2.125);
    double R = pow(tinn+random()->uniform()*(tout-tinn),1.0/2.125);
    double h = _zd*pow(R/_Rd,1.125);
    double sigma = sqrt(2.0/M_PI)*h;
    double z = random()->gauss() * sigma;
    return Position(R,phi,z,Position::CoordinateSystem::CYLINDRICAL);
}

//////////////////////////////////////////////////////////////////////

double TTauriDiskGeometry::SigmaR() const
{
    return _rho0 * _Rd * log(_Rout/_Rinn);
}

//////////////////////////////////////////////////////////////////////

double TTauriDiskGeometry::SigmaZ() const
{
    return 0.0;
}

//////////////////////////////////////////////////////////////////////
