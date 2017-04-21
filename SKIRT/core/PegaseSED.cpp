/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PegaseSED.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "System.hpp"

//////////////////////////////////////////////////////////////////////

void PegaseSED::setupSelfBefore()
{
    StellarSED::setupSelfBefore();

    // determine the resource filename based on the spectral type
    string type;
    switch (_spectralType)
    {
    case SpectralType::E:  type = "E";  break;
    case SpectralType::S0: type = "S0"; break;
    case SpectralType::Sa: type = "Sa"; break;
    case SpectralType::Sb: type = "Sb"; break;
    case SpectralType::Sc: type = "Sc"; break;
    }
    string filename = FilePaths::resource("SED/Pegase/PegaseSED_" + type + ".dat");

    // Read the data from the file into local vectors lambdav[k] and jv[k]
    const int Nlambda = 1298;
    std::ifstream file = System::ifstream(filename);
    if (! file.is_open()) throw FATALERROR("Could not open the data file " + filename);
    find<Log>()->info("Reading SED data from file " + filename + "...");
    Array lambdav(Nlambda);
    Array jv(Nlambda);
    double lambda, j, dummy;
    for (int k = 0; k<Nlambda; k++)
    {
        file >> lambda >> j >> dummy;
        lambda *= 1e-6; // conversion from micron to m
        lambdav[k] = lambda;
        jv[k] = j;
    }
    file.close();
    find<Log>()->info("File " + filename + " closed.");

    // finish up
    setEmissivities(lambdav,jv);
}

//////////////////////////////////////////////////////////////////////
