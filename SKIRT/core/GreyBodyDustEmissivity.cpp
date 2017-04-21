/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GreyBodyDustEmissivity.hpp"
#include "DustMix.hpp"
#include "PlanckFunction.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

Array GreyBodyDustEmissivity::emissivity(const DustMix* mix, const Array& Jv) const
{
    // get basic information about the wavelength grid
    WavelengthGrid* lambdagrid = find<WavelengthGrid>();
    int Nlambda = lambdagrid->numWavelengths();

    // get basic information about the dust mix
    int Npop = mix->numPopulations();

    // accumulate the emissivities at the equilibrium temperature for all populations in the dust mix
    Array ev(Nlambda);
    for (int c=0; c<Npop; c++)
    {
        double T = mix->equilibrium(Jv,c);
        PlanckFunction B(T);
        for (int ell=0; ell<Nlambda; ell++)
        {
            ev[ell] += mix->sigmaabs(ell,c) * B(lambdagrid->lambda(ell));
        }
    }

    // convert emissivity from "per hydrogen atom" to "per unit mass"
    ev /= mix->mu();
    return ev;
}

////////////////////////////////////////////////////////////////////
