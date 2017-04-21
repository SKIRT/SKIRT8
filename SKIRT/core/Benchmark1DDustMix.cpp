/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Benchmark1DDustMix.hpp"
#include "Constants.hpp"

//////////////////////////////////////////////////////////////////////

void Benchmark1DDustMix::setupSelfBefore()
{
    DustMix::setupSelfBefore();

    // get the simulation's wavelength grid
    const Array& lambdav = simlambdav();
    int Nlambda = lambdav.size();

    // create temporary vectors with the appropriate size
    Array kappaabsv(Nlambda);
    Array kappascav(Nlambda);
    Array asymmparv(Nlambda);

    // calculate the property values
    const double lambdabreak = 1e-6; // one micron
    for (int ell=0; ell<Nlambda; ell++)
    {
        double lambda = lambdav[ell];
        if (lambda <= lambdabreak)
        {
            kappaabsv[ell] = 1.0;
            kappascav[ell] = 1.0;
        }
        else
        {
            kappaabsv[ell] = lambdabreak/lambda;
            kappascav[ell] = pow(lambdabreak/lambda,4);
        }
        asymmparv[ell] = 0.0;
    }

    // add a dust population with these properties (without resampling)
    // provide a dust mass that leads to reasonable kappa values
    double Mdust = 2. / Constants::kappaV();
    addPopulation(Mdust, kappaabsv,kappascav,asymmparv);
}

//////////////////////////////////////////////////////////////////////
