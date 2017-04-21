/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "KuruczSED.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

void KuruczSED::setupSelfBefore()
{
    StellarSED::setupSelfBefore();

    // base directory for Kurucz resource library
    string filename = FilePaths::resource("SED/Kurucz/");

    // determine first portion of filenames based on desired metallicity
    if (_metallicity < -4.75) filename += "km50/km50_";
    else if (_metallicity < -4.25) filename += "km45/km45_";
    else if (_metallicity < -3.75) filename += "km40/km40_";
    else if (_metallicity < -3.25) filename += "km35/km35_";
    else if (_metallicity < -2.75) filename += "km30/km30_";
    else if (_metallicity < -2.25) filename += "km25/km25_";
    else if (_metallicity < -1.75) filename += "km20/km20_";
    else if (_metallicity < -1.25) filename += "km15/km15_";
    else if (_metallicity < -0.75) filename += "km10/km10_";
    else if (_metallicity < -0.40) filename += "km05/km05_";
    else if (_metallicity < -0.25) filename += "km03/km03_";
    else if (_metallicity < -0.15) filename += "km02/km02_";
    else if (_metallicity < -0.05) filename += "km01/km01_";
    else if (_metallicity < 0.05) filename += "kp00/kp00_";
    else if (_metallicity < 0.15) filename += "kp01/kp01_";
    else if (_metallicity < 0.25) filename += "kp02/kp02_";
    else if (_metallicity < 0.40) filename += "kp03/kp03_";
    else if (_metallicity < 0.75) filename += "kp05/kp05_";
    else filename += "kp10/kp10_";

    // determine full filenames bracketing the desired effective temperature
    double TeffL = floor(_temperature/250.)*250.;
    if (TeffL==10000.) TeffL -= 250;  // include the rightmost margin in the last bin
    double TeffR = TeffL + 250.;
    string filenameL = filename + StringUtils::toString(TeffL, 'd') + ".dat";
    string filenameR = filename + StringUtils::toString(TeffR, 'd') + ".dat";

    // open both files
    std::ifstream fileL = System::ifstream(filenameL);
    if (! fileL.is_open()) throw FATALERROR("Could not open the data file " + filenameL);
    find<Log>()->info("Reading SED data from file " + filenameL + "...");
    std::ifstream fileR = System::ifstream(filenameR);
    if (! fileR.is_open()) throw FATALERROR("Could not open the data file " + filenameR);
    find<Log>()->info("Reading SED data from file " + filenameR + "...");

    // determine the flux choice index within each file depending on desired gravity
    int mchoice = static_cast<int>(floor(2.0*_gravity+0.5));

    // construct two bracketing SEDs from the files
    int Nlambda = 1221;
    Array lambdav(Nlambda);
    Array jRv(Nlambda);
    Array jLv(Nlambda);
    int numberL, numberR;
    double lambdaL, lambdaR;
    Array fluxLgv(11);
    Array fluxRgv(11);
    for (int k=0; k<1221; k++)
    {
        fileL >> numberL >> lambdaL;
        fileR >> numberR >> lambdaR;
        for (int m=0; m<11; m++)
        {
            fileL >> fluxLgv[m];
            fileR >> fluxRgv[m];
        }
        if (lambdaL != lambdaR) throw FATALERROR("Values for lambdaL and lambdaR should be equal");
        lambdav[k] = lambdaL/1e10; // conversion from Angstrom to m
        jLv[k] = fluxLgv[mchoice];
        jRv[k] = fluxRgv[mchoice];
    }

    // close the files
    fileL.close();
    find<Log>()->info("File " + filenameL + " closed.");
    fileR.close();
    find<Log>()->info("File " + filenameR + " closed.");

    // determine the jv[k] vector by linear interpolation
    Array jv(Nlambda);
    for (int k=0; k<Nlambda; k++)
        jv[k] = jLv[k] + (_temperature-TeffL)/(TeffR-TeffL)*(jRv[k]-jLv[k]);

    // finish up
    setEmissivities(lambdav,jv);
}

////////////////////////////////////////////////////////////////////
