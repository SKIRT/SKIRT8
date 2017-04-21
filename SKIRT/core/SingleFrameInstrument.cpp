/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SingleFrameInstrument.hpp"
#include "FITSInOut.hpp"
#include "Log.hpp"
#include "PhotonPackage.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void SingleFrameInstrument::setupSelfBefore()
{
    DistantInstrument::setupSelfBefore();

    // calculate derived values
    _Nframep = _Nxp * _Nyp;
    _xpmin = _xpc - 0.5*_fovxp;
    _xpmax = _xpc + 0.5*_fovxp;
    _xpsiz = _fovxp/_Nxp;
    _ypmin = _ypc - 0.5*_fovyp;
    _ypmax = _ypc + 0.5*_fovyp;
    _ypsiz = _fovyp/_Nyp;
}

////////////////////////////////////////////////////////////////////

int SingleFrameInstrument::pixelOnDetector(const PhotonPackage* pp) const
{
    // get the position
    double x, y, z;
    pp->position().cartesian(x,y,z);

    // transform to detector coordinates using inclination, azimuth, and position angle
    double xpp = - _sinphi*x + _cosphi*y;
    double ypp = - _cosphi*_costheta*x - _sinphi*_costheta*y + _sintheta*z;
    double xp = _cospa * xpp - _sinpa * ypp;
    double yp = _sinpa * xpp + _cospa * ypp;

    // scale and round to pixel index
    int i = static_cast<int>(floor((xp-_xpmin)/_xpsiz));
    int j = static_cast<int>(floor((yp-_ypmin)/_ypsiz));
    if (i<0 || i>=_Nxp || j<0 || j>=_Nyp) return -1;
    else return i + _Nxp*j;
}

////////////////////////////////////////////////////////////////////

void SingleFrameInstrument::calibrateAndWriteDataCubes(const vector<Array*>& farrays, const vector<string>& fnames)
{
    WavelengthGrid* lambdagrid = find<WavelengthGrid>();
    int Nlambda = lambdagrid->numWavelengths();

    // calibration step 1: conversion from bolometric luminosities (units W) to monochromatic luminosities (units W/m)
    for (int ell=0; ell<Nlambda; ell++)
    {
        double dlambda = lambdagrid->dlambda(ell);
        for (int i=0; i<_Nxp; i++)
        {
            for (int j=0; j<_Nyp; j++)
            {
                size_t m = i + _Nxp*j + _Nframep*ell;
                for (Array* farr : farrays)
                {
                    if (farr->size()) (*farr)[m] /= dlambda;
                }
            }
        }
    }

    // calibration step 2: correction for the area of the pixels of the images; the units are now W/m/sr
    double xpsizang = 2.0*atan(_xpsiz/(2.0*distance()));
    double ypsizang = 2.0*atan(_ypsiz/(2.0*distance()));
    double area = xpsizang*ypsizang;
    for (Array* farr : farrays)
    {
        (*farr) /= area;
    }

    // calibration step 3: conversion of the flux per pixel from monochromatic luminosity units (W/m/sr)
    // to flux density units (W/m3/sr) by taking into account the distance
    double fourpid2 = 4.0*M_PI*distance()*distance();
    for (Array* farr : farrays)
    {
        (*farr) /= fourpid2;
    }

    // conversion from program SI units (at this moment W/m3/sr) to the correct output units;
    // we use lambda*flambda for the surface brightness (in units like W/m2/arcsec2)
    Units* units = find<Units>();
    for (int ell=0; ell<Nlambda; ell++)
    {
        // for performance reasons, determine the units scaling factor only once for each wavelength
        double factor = units->osurfacebrightness(lambdagrid->lambda(ell), 1.);
        for (int i=0; i<_Nxp; i++)
        {
            for (int j=0; j<_Nyp; j++)
            {
                size_t m = i + _Nxp*j + _Nframep*ell;
                for (Array* farr : farrays)
                {
                    if (farr->size()) (*farr)[m] *= factor;
                }
            }
        }
    }

    // Write a FITS file for each array
    for (size_t q = 0; q < farrays.size(); q++)
    {
        if (farrays[q]->size())
        {
            string filename = instrumentName() + "_" + fnames[q];
            string description = fnames[q] + " flux";
            FITSInOut::write(this, description, filename, *(farrays[q]), _Nxp, _Nyp, Nlambda,
                             units->olength(_xpsiz), units->olength(_ypsiz), units->olength(_xpc), units->olength(_ypc),
                             units->usurfacebrightness(), units->ulength());
        }
    }
}

////////////////////////////////////////////////////////////////////
