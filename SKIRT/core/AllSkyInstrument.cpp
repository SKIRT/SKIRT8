/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AllSkyInstrument.hpp"
#include "FatalError.hpp"
#include "FITSInOut.hpp"
#include "LockFree.hpp"
#include "Log.hpp"
#include "PhotonPackage.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void AllSkyInstrument::setupSelfBefore()
{
    Instrument::setupSelfBefore();

    // set number of pixels using fixed aspect ratio
    _Nx = 2 * _numPixelsY;
    _Ny = _numPixelsY;

    // unit vectors along x and y axes in the plane normal to the line crosshair-observer
    //              with y-axis oriented along projection of up direction in that plane
    Vec kn(_Ox-_Cx, _Oy-_Cy, _Oz-_Cz);
    if (kn.norm() < 1e-20) throw FATALERROR("Crosshair is too close to observer");
    Vec ku(_Ux, _Uy, _Uz);
    if (ku.norm() < 1e-20) throw FATALERROR("Upwards direction cannot be null vector");
    Vec ky = Vec::cross(kn,Vec::cross(ku,kn));
    Vec kx = Vec::cross(ky,kn);
    _bfkx = Direction(kx/kx.norm());
    _bfky = Direction(ky/ky.norm());

    // setup the transformation from world to observer coordinates

    // translate to observer position
    _transform.translate(-_Ox, -_Oy, -_Oz);

    // unit vector in direction from crosshair to observer
    kn /= kn.norm();
    double a = kn.x();
    double b = kn.y();
    double c = kn.z();

    // rotate from world to observer coordinates just as for the perspective transformation
    double v = sqrt(b*b+c*c);
    if (v > 0.3)
    {
        _transform.rotateX(c/v, -b/v);
        _transform.rotateY(v, -a);
        double k = (b*b+c*c)*_Ux - a*b*_Uy - a*c*_Uz;
        double l = c*_Uy - b*_Uz;
        double u = sqrt(k*k+l*l);
        _transform.rotateZ(l/u, -k/u);
    }
    else
    {
        v = sqrt(a*a+c*c);
        _transform.rotateY(c/v, -a/v);
        _transform.rotateX(v, -b);
        double k = c*_Ux - a*_Uz;
        double l = (a*a+c*c)*_Uy - a*b*_Ux - b*c*_Uz;
        double u = sqrt(k*k+l*l);
        _transform.rotateZ(l/u, -k/u);
    }
    // rather than flipping the z-axis as is done for the perspective transformation,
    // rotate the axes into the alignment appropriate for our purposes (z-axis up, x-axis towards crosshair)
    _transform.rotateX(0., 1.);
    _transform.rotateZ(0., -1.);

    // the data cube
    _ftotv.initialize("Instrument " + instrumentName() + " total flux", _Nx*_Ny, this);
}

////////////////////////////////////////////////////////////////////

Direction AllSkyInstrument::bfkobs(const Position& bfr) const
{
    // vector and distance from launch to observer
    Vec k = Vec(_Ox,_Oy,_Oz) - bfr;
    double d = k.norm();

    // if the distance is very small, return something arbitrary - the photon package will not be detected anyway
    if (d < _radius) return Direction();

    // otherwise return a unit vector in the direction from launch to observer
    return Direction(k/d);
}

////////////////////////////////////////////////////////////////////

Direction AllSkyInstrument::bfkx() const
{
    return _bfkx;
}

////////////////////////////////////////////////////////////////////

Direction AllSkyInstrument::bfky() const
{
    return _bfky;
}

////////////////////////////////////////////////////////////////////

void AllSkyInstrument::detect(PhotonPackage* pp)
{
    // transform launch position from world to observer coordinates
    Vec p = _transform.transform(pp->position());

    // get distance from launch position to observer
    double d = p.norm();

    // if the distance is very small, ignore the photon package
    if (d < _radius) return;

    // otherwise get the spherical coordinates from the normalized launch position
    double inc, azi;
    Direction(p/d).spherical(inc, azi);

    // convert to longitude and latitude:  -pi < lam < pi  and  -pi/2 < phi < pi/2
    double lam = -azi;  // flip east-west
    double phi = inc - M_PI_2;

    // convert longitude and latitude to viewport coordinates:  -1 < x < 1  and -1 < y < 1
    // using the Hammer-Aitoff projection
    double t = 1/sqrt(1+cos(phi)*cos(lam/2));
    double x = t*cos(phi)*sin(lam/2);
    double y = t*sin(phi);

    // convert viewport coordinates to pixel indices
    int i = max(0, min(static_cast<int>((x+1)*_Nx/2.), _Nx-1));
    int j = max(0, min(static_cast<int>((y+1)*_Ny/2.), _Ny-1));

    // determine the photon package's luminosity, attenuated for the absorption along its path to the instrument
    double taupath = opticalDepth(pp,d);
    double L = pp->luminosity() * exp(-taupath);

    // adjust the luminosity for the distance from the launch position to the instrument
    L /= d*d;

    // add the adjusted luminosity to the appropriate pixel in the data cube
    int ell = pp->ell();
    int l = i + _Nx*j;
    LockFree::add(_ftotv(ell,l), L);
}

////////////////////////////////////////////////////////////////////

void AllSkyInstrument::write()
{
    Units* units = find<Units>();
    WavelengthGrid* lambdagrid = find<WavelengthGrid>();
    int Nlambda = find<WavelengthGrid>()->numWavelengths();

    // Collect the partial data cubes into one big cube at process 0
    std::shared_ptr<Array> completeCube = _ftotv.constructCompleteCube();

    // Determine the solid angle corresponding to each pixel, assuming an area preserving projection
    // and a usage fraction in the output rectangle of pi/4
    double omega = 16. / (_Nx * _Ny);

    // Multiply each sample by lambda/dlamdba and by the constant factor 1/(4 pi Omega)
    // to obtain the surface brightness and convert to output units (such as W/m2/arcsec2)
    double front = 1. / (4.*M_PI*omega);
    for (size_t m = 0; m < completeCube->size(); m++)
    {
        size_t ell = m / (_Nx*_Ny);
        double lambda = lambdagrid->lambda(ell);
        double dlambda = lambdagrid->dlambda(ell);
        (*completeCube)[m] = units->osurfacebrightness(lambda, (*completeCube)[m]*front/dlambda);
    }

    // Determine the scale of the output map axes, i.e. in geographical coordinates
    double inc = M_PI / _Ny;

    // Write a FITS file containing the data cube
    string filename = instrumentName() + "_total";
    string description = "total flux";
    FITSInOut::write(this, description, filename, *completeCube, _Nx, _Ny, Nlambda,
                     units->oposangle(inc), units->oposangle(inc), 0., 0.,
                     units->usurfacebrightness(), units->uposangle());
}

////////////////////////////////////////////////////////////////////
