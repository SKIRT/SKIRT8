/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MarastonSED.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "System.hpp"

//////////////////////////////////////////////////////////////////////

void MarastonSED::setupSelfBefore()
{
    StellarSED::setupSelfBefore();

    // age in Gyr
    double ageGyr = _age/Constants::year()/1e9;

    // verify that property values are in the grid
    if (_metallicity<0.001 && ageGyr<1.0) throw FATALERROR("For metallicity Z<0.001, the age should be above 1 Gyr");
    if (_metallicity>0.04 && ageGyr<1.0) throw FATALERROR("For metallicity Z>0.04, the age should be above 1 Gyr");

    // fill vector with metallicities
    Array Zv(6);
    Zv[0] = 0.0001;
    Zv[1] = 0.001;
    Zv[2] = 0.01;
    Zv[3] = 0.02;
    Zv[4] = 0.04;
    Zv[5] = 0.07;

    // determine the bracketing sed filenames based on desired metallicity
    int mL = NR::locateClip(Zv,_metallicity);
    double ZL = Zv[mL];
    double ZR = Zv[mL+1];
    string fileLname = FilePaths::resource("SED/Maraston/sed.ssz");
    string fileRname = fileLname;
    int NlinesL = 0, NlinesR = 0;
    if (mL==0) {
        fileLname += "10m4.rhb";
        fileRname += "0001.rhb";
        NlinesL = 19536;
        NlinesR = 81807;
    }
    else if (mL==1) {
        fileLname += "0001.rhb";
        fileRname += "001.rhb";
        NlinesL = 81807;
        NlinesR = 81807;
    }
    else if (mL==2) {
        fileLname += "001.rhb";
        fileRname += "002.rhb";
        NlinesL = 81807;
        NlinesR = 81807;
    }
    else if (mL==3) {
        fileLname += "002.rhb";
        fileRname += "004.rhb";
        NlinesL = 81807;
        NlinesR = 81807;
    }
    else if (mL==4) {
        fileLname += "004.rhb";
        fileRname += "007.rhb";
        NlinesL = 81807;
        NlinesR = 19536;
    }

    // fill vector with ages
    string filename = FilePaths::resource("SED/Maraston/ages.dat");
    std::ifstream file = System::ifstream(filename);
    if (! file.is_open()) throw FATALERROR("Could not open the data file " + filename);
    find<Log>()->info("Reading SED data from file " + filename + "...");
    int Ntau = 67;
    Array tauv(Ntau);
    for (int l=0; l<Ntau; l++) file >> tauv[l];
    file.close();
    find<Log>()->info("File " + filename + " closed.");

    // determine the bracketing ages
    int lL = NR::locateClip(tauv,ageGyr);
    double tauL = tauv[lL];
    double tauR = tauv[lL+1];

    // initialize the wavelength and flux vectors
    const int Nlambda = 1221;
    Array lambdav(Nlambda);
    Array jLLv(Nlambda);
    Array jLRv(Nlambda);
    Array jRLv(Nlambda);
    Array jRRv(Nlambda);
    double age, ZH, lambda, j;

    // read the fluxes from the left sed file
    std::ifstream fileL = System::ifstream(fileLname);
    if (! fileL.is_open()) throw FATALERROR("Could not open the data file " + fileLname);
    find<Log>()->info("Reading SED data from file " + fileLname + "...");
    for (int k=0; k<NlinesL; k++)
    {
        fileL >> age >> ZH >> lambda >> j;
        if (k<Nlambda)
            lambdav[k] = lambda*1e-10; // conversion from Angstrom to m
        if (age==tauL)
            jLLv[k%Nlambda] = j;
        else if (age==tauR)
            jLRv[k%Nlambda] = j;
    }
    fileL.close();
    find<Log>()->info("File " + fileLname + " closed.");

    // read the fluxes from the right sed file
    std::ifstream fileR = System::ifstream(fileRname);
    if (! fileR.is_open()) throw FATALERROR("Could not open the data file " + fileRname);
    find<Log>()->info("Reading SED data from file " + fileRname + "...");
    for (int k=0; k<NlinesR; k++)
    {
        fileR >> age >> ZH >> lambda >> j;
        if (age==tauL)
            jRLv[k%Nlambda] = j;
        else if (age==tauR)
            jRRv[k%Nlambda] = j;
    }
    fileR.close();
    find<Log>()->info("File " + fileRname + " closed.");

    // interpolate
    Array jv(Nlambda);
    double p = (_metallicity-ZL)/(ZR-ZL);
    double q = (ageGyr-tauL)/(tauR-tauL);
    for (int k=0; k<Nlambda; k++)
        jv[k] = (1.0-p)*(1.0-q)*jLLv[k]
                + p*(1.0-q)*jRLv[k]
                + (1.0-p)*q*jLRv[k]
                + p*q*jRRv[k];

    // finish up
    setEmissivities(lambdav,jv);
}

////////////////////////////////////////////////////////////////////
