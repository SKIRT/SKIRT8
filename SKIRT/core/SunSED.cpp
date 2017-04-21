/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SunSED.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "System.hpp"

//////////////////////////////////////////////////////////////////////

void SunSED::setupSelfBefore()
{
    StellarSED::setupSelfBefore();

    // read the data from a resource file into local vectors lambdav[k] and jv[k]
    string filename = FilePaths::resource("SED/Sun/SunSED.dat");
    std::ifstream file = System::ifstream(filename);
    if (! file.is_open()) throw FATALERROR("Could not open the data file " + filename);
    find<Log>()->info("Reading SED data from file " + filename + "...");
    string line;
    getline(file,line);
    int Nlambda;
    file >> Nlambda;
    Array lambdav(Nlambda);
    Array jv(Nlambda);
    double lambda;
    for (int k=0; k<Nlambda; k++)
    {
        file >> lambda >> jv[k];
        lambdav[k] = lambda/1e6; // conversion from micron to m;
    }
    file.close();
    find<Log>()->info("File " + filename + " closed.");

    // finish up
    setEmissivities(lambdav,jv);
}

//////////////////////////////////////////////////////////////////////
