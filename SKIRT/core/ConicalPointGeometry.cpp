/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ConicalPointGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void ConicalPointGeometry::setupSelfBefore()
{
    Geometry::setupSelfBefore();
    _cosDelta = cos(_Delta);
}

//////////////////////////////////////////////////////////////////////

int ConicalPointGeometry::dimension() const
{
    return 2;
}

//////////////////////////////////////////////////////////////////////

double ConicalPointGeometry::density(Position bfr) const
{
    return bfr.radius() == 0 ? std::numeric_limits<double>::infinity() : 0;
}

//////////////////////////////////////////////////////////////////////

Position ConicalPointGeometry::generatePosition() const
{
    return Position();
}

//////////////////////////////////////////////////////////////////////

double ConicalPointGeometry::SigmaX() const
{
    return std::numeric_limits<double>::infinity();
}

//////////////////////////////////////////////////////////////////////

double ConicalPointGeometry::SigmaY() const
{
    return std::numeric_limits<double>::infinity();
}

//////////////////////////////////////////////////////////////////////

double ConicalPointGeometry::SigmaZ() const
{
    return std::numeric_limits<double>::infinity();
}

//////////////////////////////////////////////////////////////////////

double ConicalPointGeometry::probabilityForDirection(int /*ell*/, Position bfr, Direction bfk) const
{
    if (bfr.radius()>0.0)
        throw FATALERROR("the angular probability function is not defined for positions besides the origin");
    double theta, phi;
    bfk.spherical(theta, phi);
    if (theta > _Delta && theta < M_PI-_Delta)
        return 0.0;
    else
        return 1.0/(1.0-_cosDelta);
}

//////////////////////////////////////////////////////////////////////

Direction ConicalPointGeometry::generateDirection(int /*ell*/, Position bfr) const
{
    if (bfr.radius()>0.0)
        throw FATALERROR("no directions should be generated at positions besides the origin");
    double theta = 0.0;
    double X = random()->uniform();
    if (X<0.5)
        theta = acos(1.0-2.0*X*(1.0-_cosDelta));
    else
        theta = acos(1.0-2.0*_cosDelta-2.0*X*(1.0-_cosDelta));
    double phi = 2.0*M_PI*random()->uniform();
    return Direction(theta,phi);
}

//////////////////////////////////////////////////////////////////////
