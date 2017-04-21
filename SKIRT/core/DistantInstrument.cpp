/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DistantInstrument.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "PhotonPackage.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void DistantInstrument::setupSelfBefore()
{
    Instrument::setupSelfBefore();

    // calculate sine and cosine for our angles
    _costheta = cos(_inclination);
    _sintheta = sin(_inclination);
    _cosphi = cos(_azimuth);
    _sinphi = sin(_azimuth);
    _cospa = cos(_positionAngle);
    _sinpa = sin(_positionAngle);

    // calculate relevant directions
    _bfkobs = Direction(_inclination,_azimuth);
    _bfkx = Direction( + _cosphi*_costheta*_sinpa - _sinphi*_cospa,
                       + _sinphi*_costheta*_sinpa + _cosphi*_cospa,
                       - _sintheta*_sinpa );
    _bfky = Direction( - _cosphi*_costheta*_cospa - _sinphi*_sinpa,
                       - _sinphi*_costheta*_cospa + _cosphi*_sinpa,
                       + _sintheta*_cospa );
}

////////////////////////////////////////////////////////////////////

Direction DistantInstrument::bfkobs(const Position& /*bfr*/) const
{
    return _bfkobs;
}

////////////////////////////////////////////////////////////////////

Direction DistantInstrument::bfkx() const
{
    return _bfkx;
}

////////////////////////////////////////////////////////////////////

Direction DistantInstrument::bfky() const
{
    return _bfky;
}

////////////////////////////////////////////////////////////////////

void DistantInstrument::calibrateAndWriteSEDs(const vector<Array*>& Farrays, const vector<string>& Fnames)
{
    PeerToPeerCommunicator* comm = find<PeerToPeerCommunicator>();
    if (comm->rank()) return;

    WavelengthGrid* lambdagrid = find<WavelengthGrid>();
    int Nlambda = find<WavelengthGrid>()->numWavelengths();

    // calibration step 1: conversion from bolometric luminosities (units W) to monochromatic luminosities (units W/m)

    for (int ell=0; ell<Nlambda; ell++)
    {
        double dlambda = lambdagrid->dlambda(ell);
        for (Array* Farr : Farrays)
        {
            if (Farr->size()) (*Farr)[ell] /= dlambda;
        }
    }

    // calibration step 2: conversion of the integrated flux from monochromatic luminosity units (W/m) to
    //                     flux density units (W/m3) by taking into account the distance

    double fourpid2 = 4.0*M_PI*_distance*_distance;
    for (Array* Farr : Farrays)
    {
        (*Farr) /= fourpid2;
    }

    // write a text file for easy SED plotting

    Units* units = find<Units>();
    TextOutFile sedfile(this, instrumentName() + "_sed", "SED");

    // Write the header
    sedfile.addColumn("lambda (" + units->uwavelength() + ")", 'e', 8);
    for (size_t q = 0; q < Farrays.size(); q++)
    {
        sedfile.addColumn(Fnames[q] + "; " + units->sfluxdensity() + " " + "(" + units->ufluxdensity() + ")", 'e', 8);
    }

    // Write the body
    for (int ell=0; ell<Nlambda; ell++)
    {
        double lambda = lambdagrid->lambda(ell);
        vector<double> values;
        values.push_back(units->owavelength(lambda));
        for (Array* Farr : Farrays)
        {
            values.push_back(Farr->size() ? units->ofluxdensity(lambda, (*Farr)[ell]) : 0.);
        }
        sedfile.writeRow(values);
    }
}

////////////////////////////////////////////////////////////////////
