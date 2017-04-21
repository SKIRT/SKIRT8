/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BlackBodySED.hpp"
#include "PlanckFunction.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void BlackBodySED::setupSelfBefore()
{
    StellarSED::setupSelfBefore();

    WavelengthGrid* lambdagrid = find<WavelengthGrid>();
    int Nlambda = lambdagrid->numWavelengths();
    PlanckFunction B(_temperature);

    Array Lv(Nlambda);
    for (int ell=0; ell<Nlambda; ell++)
    {
        // we take averages over each bin
        int N = 100;
        double loglambdamin = log10(lambdagrid->lambdamin(ell));
        double loglambdamax = log10(lambdagrid->lambdamax(ell));
        double dloglambda = (loglambdamax-loglambdamin)/N;
        double sum = 0;
        for (int i=0; i<=N; i++)
        {
            double weight = 1.0;
            if (i==0 || i==N) weight = 0.5;
            double loglambda = loglambdamin + i*dloglambda;
            double lambda = pow(10,loglambda);
            sum += weight*B(lambda)*lambda;
        }
        Lv[ell] = sum * M_LN10 * dloglambda;
    }

    // finish up
    setLuminosities(Lv);
}

//////////////////////////////////////////////////////////////////////
