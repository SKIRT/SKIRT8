/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InstrumentFrame.hpp"
#include "FITSInOut.hpp"
#include "LockFree.hpp"
#include "MultiFrameInstrument.hpp"
#include "PhotonPackage.hpp"
#include "StellarSystem.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void InstrumentFrame::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // calculate derived values
    _Nframep = _Nxp * _Nyp;
    _xpmin = _xpc - 0.5*_fovxp;
    _xpmax = _xpc + 0.5*_fovxp;
    _xpsiz = _fovxp/_Nxp;
    _ypmin = _ypc - 0.5*_fovyp;
    _ypmax = _ypc + 0.5*_fovyp;
    _ypsiz = _fovyp/_Nyp;

    // copy information from parent instrument
    _instrument = find<MultiFrameInstrument>();
    _writeTotal = _instrument->writeTotal();
    _writeStellarComps = _instrument->writeStellarComps();
    _distance = _instrument->distance();
    double inclination = _instrument->inclination();
    double azimuth = _instrument->azimuth();
    double positionangle = _instrument->positionAngle();
    _costheta = cos(inclination);
    _sintheta = sin(inclination);
    _cosphi = cos(azimuth);
    _sinphi = sin(azimuth);
    _cospa = cos(positionangle);
    _sinpa = sin(positionangle);

    // initialize pixel frame(s)
    if (_writeTotal) _ftotv.resize(_Nframep);
    if (_writeStellarComps) _fcompvv.resize(find<StellarSystem>()->numComponents(), _Nframep);
}

////////////////////////////////////////////////////////////////////

int InstrumentFrame::pixelOnDetector(const PhotonPackage* pp) const
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

void InstrumentFrame::detect(PhotonPackage* pp)
{
    int l = pixelOnDetector(pp);
    if (l >= 0)
    {
        double L = pp->luminosity();
        double taupath = _instrument->opticalDepth(pp);
        double extf = exp(-taupath);
        double Lextf = L*extf;

        if (_writeTotal) LockFree::add(_ftotv[l], Lextf);
        if (_writeStellarComps && pp->isStellar()) LockFree::add(_fcompvv(pp->stellarCompIndex(),l), Lextf);
    }
}

////////////////////////////////////////////////////////////////////

void InstrumentFrame::calibrateAndWriteData(int ell)
{
    // construct list of data cube pointers and the corresponding file names
    vector<Array*> farrays;
    vector<string> fnames;

    if (_writeTotal)
    {
        farrays.push_back(&_ftotv);
        fnames.push_back("total");
    }
    if (_writeStellarComps)
    {
        int Ncomp = find<StellarSystem>()->numComponents();
        for (int k=0; k<Ncomp; k++)
        {
            farrays.push_back( &(_fcompvv[k]) );
            fnames.push_back("stellar_" + std::to_string(k));
        }
    }

    // Sum the flux arrays element-wise across the different processes
    _instrument->sumResults(farrays);

    // calibrate and output the arrays
    calibrateAndWriteDataFrames(ell, farrays, fnames);
}

////////////////////////////////////////////////////////////////////

void InstrumentFrame::calibrateAndWriteDataFrames(int ell, const vector<Array*>& farrays, const vector<string>& fnames)
{
    Units* units = find<Units>();
    WavelengthGrid* lambdagrid = find<WavelengthGrid>();

    // conversion from bolometric luminosities (units W) to monochromatic luminosities (units W/m)
    // --> divide by delta-lambda
    double dlambda = lambdagrid->dlambda(ell);

    // correction for the area of the pixels of the images; the units are now W/m/sr
    // --> divide by area
    double xpsizang = 2.0*atan(_xpsiz/(2.0*_distance));
    double ypsizang = 2.0*atan(_ypsiz/(2.0*_distance));
    double area = xpsizang*ypsizang;

    // calibration step 3: conversion of the flux per pixel from monochromatic luminosity units (W/m/sr)
    // to flux density units (W/m3/sr) by taking into account the distance
    // --> divide by fourpid2
    double fourpid2 = 4.0*M_PI*_distance*_distance;

    // conversion from program SI units (at this moment W/m3/sr) to the correct output units
    // --> multiply by unit conversion factor
    double unitfactor = units->osurfacebrightness(lambdagrid->lambda(ell), 1.);

    // Perform the conversion, in place
    for (Array* farr : farrays)
    {
        (*farr) *= (unitfactor / (dlambda * area * fourpid2));
    }

    // Write a FITS file for each array
    for (size_t q = 0; q < farrays.size(); q++)
    {
        string filename = _instrument->instrumentName() + "_" + fnames[q] + "_" + std::to_string(ell);
        string description = fnames[q] + " flux " + std::to_string(ell);
        FITSInOut::write(this, description, filename, *(farrays[q]), _Nxp, _Nyp, 1,
                         units->olength(_xpsiz), units->olength(_ypsiz), units->olength(_xpc), units->olength(_ypc),
                         units->usurfacebrightness(), units->ulength());
    }
}

////////////////////////////////////////////////////////////////////
