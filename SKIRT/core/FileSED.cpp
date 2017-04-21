/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileSED.hpp"
#include "FilePaths.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

void FileSED::setupSelfBefore()
{
    StellarSED::setupSelfBefore();

    // read the data from the specified file into local vectors lambdav[k] and jv[k]
    string filename = find<FilePaths>()->input(_filename);
    std::ifstream file = System::ifstream(filename);
    if (! file.is_open()) throw FATALERROR("Could not open the data file " + filename);
    find<Log>()->info("Reading SED data from file " + filename + "...");

    int Nlambda;
    file >> Nlambda;
    Array lambdav(Nlambda);
    Array jv(Nlambda);
    double lambda;
    for (int k=0; k<Nlambda; k++)
    {
        file >> lambda >> jv[k];
        lambdav[k] = lambda/1e6; // conversion from micron to m
    }
    file.close();
    find<Log>()->info("File " + filename + " closed.");

    // finish up
    setEmissivities(lambdav,jv);
}

//////////////////////////////////////////////////////////////////////
