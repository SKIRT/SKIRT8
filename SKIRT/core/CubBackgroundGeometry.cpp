/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CubBackgroundGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void CubBackgroundGeometry::setupSelfBefore()
{
    Geometry::setupSelfBefore();

     _h = 0.5 * _edgeLength;
}

//////////////////////////////////////////////////////////////////////

int CubBackgroundGeometry::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

double CubBackgroundGeometry::density(Position bfr) const
{
    double x, y, z;
    bfr.cartesian(x,y,z);
    double ax = fabs(x);
    double ay = fabs(y);
    double az = fabs(z);
    return ( (ax==_h && ay<=_h && az<=_h)
          || (ay==_h && ax<=_h && az<=_h)
          || (az==_h && ax<=_h && ay<=_h) ) ?
                std::numeric_limits<double>::infinity() : 0.0;
}

//////////////////////////////////////////////////////////////////////

Position CubBackgroundGeometry::generatePosition() const
{
    double t1 = _h*(2.0*random()->uniform()-1.0);
    double t2 = _h*(2.0*random()->uniform()-1.0);
    double X = random()->uniform();
    if (X<1.0/6.0)
        return Position(-_h,t1,t2);
    else if (X<1.0/3.0)
        return Position(_h,t1,t2);
    else if (X<0.5)
        return Position(t1,-_h,t2);
    else if (X<2.0/3.0)
        return Position(t1,_h,t2);
    else if (X<5.0/6.0)
        return Position(t1,t2,-_h);
    else
        return Position(t1,t2,_h);
}

//////////////////////////////////////////////////////////////////////

double CubBackgroundGeometry::SigmaX() const
{
    return 1.0/(12.0*_h*_h);
}

//////////////////////////////////////////////////////////////////////

double CubBackgroundGeometry::SigmaY() const
{
    return 1.0/(12.0*_h*_h);
}

//////////////////////////////////////////////////////////////////////

double CubBackgroundGeometry::SigmaZ() const
{
    return 1.0/(12.0*_h*_h);
}

//////////////////////////////////////////////////////////////////////

double CubBackgroundGeometry::probabilityForDirection(int /*ell*/, Position bfr, Direction bfk) const
{
    double eps = 1e-8;
    double x, y, z;
    bfr.cartesian(x,y,z);
    double ax = fabs(x);
    double ay = fabs(y);
    double az = fabs(z);
    Vec norm;
    if (fabs(x/_h+1.0)<eps && ay<=_h && az<=_h)
        norm = Vec(-1.,0.,0.);
    else if (fabs(x/_h-1.0)<eps && ay<=_h && az<=_h)
        norm = Vec(1.,0.,0.);
    else if (fabs(y/_h+1.0)<eps && ax<=_h && az<=_h)
        norm = Vec(0.,-1.,0.);
    else if (fabs(y/_h-1.0)<eps && ax<=_h && az<=_h)
        norm = Vec(0.,1.,0.);
    else if (fabs(z/_h+1.0)<eps && ax<=_h && ay<=_h)
        norm = Vec(0.,0.,-1.);
    else if (fabs(z/_h-1.0)<eps && ax<=_h && ay<=_h)
        norm = Vec(0.,0.,1.);
    else
        throw FATALERROR("the directional probability function is not defined for positions not on the background cube");
    double costhetap = Vec::dot(norm,bfk);
    if (costhetap > 0.0)
        return 0.0;
    else
        return -4.0*costhetap;
}

//////////////////////////////////////////////////////////////////////

Direction CubBackgroundGeometry::generateDirection(int /*ell*/, Position bfr) const
{
    // picking a random (theta',phi')
    double thetap = M_PI-acos(sqrt(random()->uniform()));
    double phip = 2.0*M_PI*random()->uniform();
    Direction bfkp(thetap,phip);
    double kpx, kpy, kpz;
    bfkp.cartesian(kpx,kpy,kpz);

    // conversion to the regular coordinate system
    double eps = 1e-8;
    double x, y, z;
    bfr.cartesian(x,y,z);
    double ax = fabs(x);
    double ay = fabs(y);
    double az = fabs(z);
    if (fabs(x/_h+1.0)<eps && ay<=_h && az<=_h)
        return Direction(-kpz,-kpy,-kpx);
    else if (fabs(x/_h-1.0)<eps && ay<=_h && az<=_h)
        return Direction(kpz,kpy,-kpx);
    else if (fabs(y/_h+1.0)<eps && ax<=_h && az<=_h)
        return Direction(kpy,-kpz,-kpx);
    else if (fabs(y/_h-1.0)<eps && ax<=_h && az<=_h)
        return Direction(-kpy,kpz,-kpx);
    else if (fabs(z/_h+1.0)<eps && ax<=_h && ay<=_h)
        return Direction(-kpx,kpy,-kpz);
    else if (fabs(z/_h-1.0)<eps && ax<=_h && ay<=_h)
        return Direction(kpx,kpy,kpz);
    else
        throw FATALERROR("cannot generate directions for positions not on the background cube");
}

//////////////////////////////////////////////////////////////////////
