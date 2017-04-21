/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TrustPolarizedMeanDustMix.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "StringUtils.hpp"
#include "System.hpp"

//////////////////////////////////////////////////////////////////////

void TrustPolarizedMeanDustMix::setupSelfBefore()
{
    DustMix::setupSelfBefore();

    // create temporary vectors with the wavelength resolution of the file
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

    // for resampling: get the simulation's wavelength grid and its length
    const Array& lambdagridv = this->simlambdav();
    int SimNlambda = lambdagridv.size();

    // create temporary vectors and tables with the appropriate size
    int Ntheta = 181;
    Array S11v(Nlambda), S12v(Nlambda), S33v(Nlambda), S34v(Nlambda);
    Array S11RSv(SimNlambda), S12RSv(SimNlambda), S33RSv(SimNlambda), S34RSv(SimNlambda);
    Table<2> S11vv(SimNlambda,Ntheta), S12vv(SimNlambda,Ntheta), S33vv(SimNlambda,Ntheta), S34vv(SimNlambda,Ntheta);

    // read the Mueller matrix components
    // hierarchy: 181 files by angles, 1201 lines by wavelengths
    string filenameX = "DustMix/TrustMDMScatMatrix/ZDA_BARE_GR_S_ESM_xxxdeg.dat";
    find<Log>()->info("Reading Mueller Matrix components from files " + filenameX + "...");
    for (int t=0; t<Ntheta; t++)
    {
        string filenameN = StringUtils::replace(filenameX, "xxx", StringUtils::toString(t,'d',0,3,'0'));
        string filepathN = FilePaths::resource(filenameN);
        std::ifstream file = System::ifstream(filepathN);
        if (!file.is_open()) throw FATALERROR("Could not open the data file " + filepathN);

        string line;
        while (file.peek() == '#') getline(file,line);
        double lambda, S11, S12, S33, S34;

        // read in tables
        for (int k=0; k<Nlambda; k++)
        {
            file >> lambda >> S11 >> S12 >> S33 >> S34;
            S11v[k] = S11;
            S12v[k] = S12;
            S33v[k] = S33;
            S34v[k] = S34;
        }

        // resample wavelengths
        S11RSv = NR::resample<NR::interpolateLogLin>(lambdagridv,lambdav,S11v);
        S12RSv = NR::resample<NR::interpolateLogLin>(lambdagridv,lambdav,S12v);
        S33RSv = NR::resample<NR::interpolateLogLin>(lambdagridv,lambdav,S33v);
        S34RSv = NR::resample<NR::interpolateLogLin>(lambdagridv,lambdav,S34v);

        for (int k=0; k<SimNlambda; k++)
        {
            S11vv(k,t) = S11RSv[k];
            S12vv(k,t) = S12RSv[k];
            S33vv(k,t) = S33RSv[k];
            S34vv(k,t) = S34RSv[k];
        }
    }

    file.close();
    find<Log>()->info("Files " + filenameX + " closed.");

    // determine the dust mass per H nucleon (cross sections in file are also per nucleon)
    const double mu = 1.434e-29;  // in kg

    // add a dust population with these properties
    addPopulation(mu, lambdav,sigmaabsv,sigmascav,asymmparv);
    addPolarization(S11vv, S12vv, S33vv, S34vv);
}

//////////////////////////////////////////////////////////////////////
