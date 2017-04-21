/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StarburstSED.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "StringUtils.hpp"
#include "System.hpp"

//////////////////////////////////////////////////////////////////////

void StarburstSED::setupSelfBefore()
{
    StellarSED::setupSelfBefore();

    // open the resource file and skip the header
    string filename = FilePaths::resource("SED/Starburst/StarburstSED.dat");
    std::ifstream file = System::ifstream(filename);
    if (! file.is_open()) throw FATALERROR("Could not open the data file " + filename);
    find<Log>()->info("Reading SED data from file " + filename + "...");
    string line;
    for (int i=0; i<6; i++) getline(file,line);

    // get the metallicity range and the number of wavelengths
    int NZ, Nlambda;
    file >> NZ >> Nlambda;
    Array Zv(NZ);
    for (int l=0; l<NZ; l++) file >> Zv[l];
    int lL = NR::locateClip(Zv,_metallicity);
    double ZL = Zv[lL];
    double ZR = Zv[lL+1];

    // read the data from the file into local vectors
    Array lambdav(Nlambda);
    Array jv(Nlambda);
    Array logjLv(Nlambda);
    Array logjRv(Nlambda);
    double lambda, dummy;
    for (int k=0; k<Nlambda; k++)
    {
        file >> lambda;
        for (int l=0; l<lL; l++) file >> dummy;
        file >> logjLv[k] >> logjRv[k];
        for (int l=lL+2; l<NZ; l++) file >> dummy;
        lambdav[k] = lambda/1e10; // conversion from A to m
    }
    file.close();
    find<Log>()->info("File " + filename + " closed.");

    // interpolate linearly in log space
    for (int k=0; k<Nlambda; k++)
        jv[k] = pow(10., NR::interpolateLinLin(_metallicity, ZL, ZR, logjLv[k], logjRv[k]));

    // finish up
    setEmissivities(lambdav,jv);
}

////////////////////////////////////////////////////////////////////
