/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PinteBenchmarkGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

void PinteBenchmarkGeometry::setupSelfBefore()
{
    AxGeometry::setupSelfBefore();

    // set property values
    _Rinn = 0.1 * Constants::AU();
    _Rout = 400.0 * Constants::AU();
    _Rd = 100.0 * Constants::AU();
    _zd = 10.0 * Constants::AU();

    // calculate normalization factor
    _rho0 = 1. / ( 8./5. * pow(2.*M_PI,3./2.) * _Rd * _Rd *_zd * (pow(_Rout/_Rd,5./8.) - pow(_Rinn/_Rd,5./8.)) );
}

//////////////////////////////////////////////////////////////////////

double PinteBenchmarkGeometry::density(double R, double z) const
{
    if (R<_Rinn || R>_Rout) return 0.;
    double x = (z/_zd) / pow(R/_Rd,9./8.);
    return _rho0 * pow(R/_Rd,-5./2.) * exp(-0.5*x*x);
}

//////////////////////////////////////////////////////////////////////

Position PinteBenchmarkGeometry::generatePosition() const
{
    double phi = 2. * M_PI * random()->uniform();
    double tinn = pow(_Rinn,5./8.);
    double tout = pow(_Rout,5./8.);
    double R = pow(tinn+random()->uniform()*(tout-tinn),8./5.);
    double sigma = _zd*pow(R/_Rd,9./8.);
    double z = random()->gauss() * sigma;
    return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
}

//////////////////////////////////////////////////////////////////////

double PinteBenchmarkGeometry::SigmaR() const
{
    return 2./3. * _rho0 * _Rd * (pow(_Rinn/_Rd,-3./2.)-pow(_Rout/_Rd,-3./2.));
}

//////////////////////////////////////////////////////////////////////

double PinteBenchmarkGeometry::SigmaZ() const
{
    return 0.;
}

//////////////////////////////////////////////////////////////////////
