/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PinteBenchmarkDustMix.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "System.hpp"

//////////////////////////////////////////////////////////////////////

void PinteBenchmarkDustMix::setupSelfBefore()
{
    DustMix::setupSelfBefore();

    // --- basic cross sections ---

    // create temporary vectors with the appropriate size
    int Nlambda = simlambdav().size();
    Array sigmaabsv(Nlambda), sigmascav(Nlambda), asymmparv(Nlambda);

    // these values are taken from http://ipag.osug.fr/~pintec/benchmark/optical_properties.shtml
    double albedo = 0.6474786;
    double asymmpar = 0.6296066;
    double sigmaext = 4752.16 * 0.1; // convert from cm^2/g to m^2/kg
    sigmaabsv = (1.-albedo) * sigmaext;
    sigmascav = albedo * sigmaext;
    asymmparv = asymmpar;

    // add a dust population with these properties
    addPopulation(1., sigmaabsv, sigmascav, asymmparv);

    // --- Mueller matrix ---

    // create temporary vectors and tables with the appropriate size
    int Ntheta = 180;
    Table<2> S11vv(Nlambda,Ntheta), S12vv(Nlambda,Ntheta), S33vv(Nlambda,Ntheta), S34vv(Nlambda,Ntheta);

    // read the Mueller matrix components; the file has 180 angles 0.5, 1.5, ..., 179.5
    // while SKIRT expects a table from angle 0 going to and including 180; we ignore this discrepancy.
    string filename = FilePaths::resource("DustMix/Pinte2009_1mu_mueller_matrix.dat");
    find<Log>()->info("Reading Mueller Matrix components from file " + filename + "...");
    std::ifstream file = System::ifstream(filename);
    if (!file.is_open()) throw FATALERROR("Could not open the data file " + filename);

    for (int t=0; t<Ntheta; t++)
    {
        string line;
        while (file.peek() == '#') getline(file,line);
        double theta, S11, S12, S33, S34;

        // read in tables
        file >> theta >> S11 >> S12 >> S33 >> S34;

        for (int ell=0; ell<Nlambda; ell++)
        {
            S11vv(ell,t) = S11;
            S12vv(ell,t) = S12;
            S33vv(ell,t) = S33;
            S34vv(ell,t) = -S34; // input values are given in the Legacy convention!
        }
    }

    file.close();
    find<Log>()->info("File " + filename + "..." + " closed.");

    // add the polarization properties
    addPolarization(S11vv, S12vv, S33vv, S34vv);
}

//////////////////////////////////////////////////////////////////////
