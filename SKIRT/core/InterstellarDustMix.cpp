/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InterstellarDustMix.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "System.hpp"
#include <fstream>

//////////////////////////////////////////////////////////////////////

void InterstellarDustMix::setupSelfBefore()
{
    DustMix::setupSelfBefore();

    // create temporary vectors with the appropriate size
    const int Nlambda = 1064;
    Array lambdav(Nlambda);
    Array kappaabsv(Nlambda);
    Array kappascav(Nlambda);
    Array asymmparv(Nlambda);

    // read the raw data from the resource file into the temporary vectors
    string filename = FilePaths::resource("DustMix/InterstellarDustMix.dat");
    std::ifstream file = System::ifstream(filename);
    if (!file.is_open()) throw FATALERROR("Could not open the data file " + filename);
    find<Log>()->info("Reading dust mix properties from file " + filename + "...");
    string line;
    while (file.peek() == '#') getline(file,line);
    double lambda, albedo, asymmpar, Cext, Kabs, cos2;
    for (int k=Nlambda-1; k>=0; k--)
    {
        file >> lambda >> albedo >> asymmpar >> Cext >> Kabs >> cos2;

        lambda *= 1e-6;  // lambda in file is in micron, we need lambda in m
        Kabs *= 1e-1;    // Kabs in file is in cm2/g, we need Kabs in m2/kg

        lambdav[k] = lambda;
        kappaabsv[k] = Kabs;
        kappascav[k] = Kabs*albedo/(1.0-albedo);
        asymmparv[k] = asymmpar;
    }
    file.close();
    find<Log>()->info("File " + filename + " closed.");

    // add a dust population with these properties
    // since we already have kappa's, the dust mass is set to unity
    addPopulation(1., lambdav,kappaabsv,kappascav,asymmparv);
}

//////////////////////////////////////////////////////////////////////
