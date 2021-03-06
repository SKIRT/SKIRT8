/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TrustMeanDustMix.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "System.hpp"

//////////////////////////////////////////////////////////////////////

void TrustMeanDustMix::setupSelfBefore()
{
    DustMix::setupSelfBefore();

    // create temporary vectors with the appropriate size
    const int Nlambda = 1201;
    Array lambdav(Nlambda);
    Array sigmaabsv(Nlambda);
    Array sigmascav(Nlambda);
    Array asymmparv(Nlambda);

    // read the raw data from the resource file into the temporary vectors
    string filename = FilePaths::resource("DustMix/TrustMeanDustMix.dat");
    std::ifstream file = System::ifstream(filename);
    if (!file.is_open()) throw FATALERROR("Could not open the data file " + filename);
    find<Log>()->info("Reading dust mix properties from file " + filename + "...");
    string line;
    while (file.peek() == '#') getline(file,line);
    double lambda, Cabs, Csca, tauNH, albedo, asymmpar;
    for (int k=0; k<Nlambda; k++)
    {
        file >> lambda >> Cabs >> Csca >> tauNH >> albedo >> asymmpar;
        lambdav[k] = lambda * 1e-6;   // lambda in file is in micron, we need lambda in m;
        double sigmaext = tauNH * 1e-4; // tauNH in file is in cm^2 H^-1, we need sigmaext in m^2 H^-1
        sigmaabsv[k] = (1.-albedo) * sigmaext;
        sigmascav[k] = albedo * sigmaext;
        asymmparv[k] = asymmpar;
    }
    file.close();
    find<Log>()->info("File " + filename + " closed.");

    // determine the dust mass per H nucleon (cross sections in file are also per nucleon)
    const double mu = 1.434e-29;  // in kg

    // add a dust population with these properties
    addPopulation(mu, lambdav,sigmaabsv,sigmascav,asymmparv);
}

//////////////////////////////////////////////////////////////////////
