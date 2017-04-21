/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpheBackgroundGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

int SpheBackgroundGeometry::dimension() const
{
    return 1;
}

//////////////////////////////////////////////////////////////////////

double SpheBackgroundGeometry::density(Position bfr) const
{
    return (bfr.radius() == _rbg) ?
                std::numeric_limits<double>::infinity() : 0.0;
}

//////////////////////////////////////////////////////////////////////

Position SpheBackgroundGeometry::generatePosition() const
{
    return Position(_rbg,random()->direction());
}

//////////////////////////////////////////////////////////////////////

double SpheBackgroundGeometry::SigmaX() const
{
    return 1.0/(2.0*M_PI*_rbg*_rbg);
}

//////////////////////////////////////////////////////////////////////

double SpheBackgroundGeometry::SigmaY() const
{
    return 1.0/(2.0*M_PI*_rbg*_rbg);
}

//////////////////////////////////////////////////////////////////////

double SpheBackgroundGeometry::SigmaZ() const
{
    return 1.0/(2.0*M_PI*_rbg*_rbg);
}

//////////////////////////////////////////////////////////////////////

double SpheBackgroundGeometry::probabilityForDirection(int /*ell*/, Position bfr, Direction bfk) const
{
    if (fabs(bfr.radius()/_rbg-1.0) > 1e-8)
        throw FATALERROR("the directional probability function is not defined for positions not on the background sphere");
    double costhetap = Vec::dot(bfr,bfk)/_rbg;
    if (costhetap > 0.0)
        return 0.0;
    else
        return -4.0*costhetap;
}

//////////////////////////////////////////////////////////////////////

Direction SpheBackgroundGeometry::generateDirection(int /*ell*/, Position bfr) const
{
    if (fabs(bfr.radius()/_rbg-1.0) > 1e-8)
        throw FATALERROR("cannot generate directions for positions not on the SpheBackground sphere");

    // picking a random (theta',phi')
    double thetap = M_PI-acos(sqrt(random()->uniform()));
    double phip = 2.0*M_PI*random()->uniform();
    Direction bfkp(thetap,phip);
    double kpx, kpy, kpz;
    bfkp.cartesian(kpx,kpy,kpz);

    // conversion to the regular coordinate system
    double r, theta, phi;
    bfr.spherical(r,theta,phi);
    double costheta = cos(theta);
    double sintheta = sin(theta);
    double cosphi = cos(phi);
    double sinphi = sin(phi);
    double kx = costheta*cosphi*kpx - sinphi*kpy + sintheta*cosphi*kpz;
    double ky = costheta*sinphi*kpx + cosphi*kpy + sintheta*sinphi*kpz;
    double kz = -sintheta*kpx + costheta*kpz;
    return Direction(kx,ky,kz);
}

//////////////////////////////////////////////////////////////////////
