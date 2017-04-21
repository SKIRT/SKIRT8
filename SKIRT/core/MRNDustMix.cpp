/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MRNDustMix.hpp"
#include "DraineGraphiteGrainComposition.hpp"
#include "DraineSilicateGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // MRN grain size distributions with values taken from Weingartner & Draine (2001, ApJ, 548, 296) page 296
    const double amin = 5e-9;          // 50 Armstrong
    const double amax = 250e-9;        // 0.25 micron
    const double Cg = pow(10,-25.13) * 1e-5;  // convert value in cm^2.5 to m^2.5
    const double Cs = pow(10,-25.11) * 1e-5;  // convert value in cm^2.5 to m^2.5

    double dnda_gra(double a)
    {
        return Cg * pow(a,-3.5);
    }

    double dnda_sil(double a)
    {
        return Cs * pow(a,-3.5);
    }
}

//////////////////////////////////////////////////////////////////////

void MRNDustMix::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    addPopulations(new DraineGraphiteGrainComposition(this), amin, amax, &dnda_gra, _numGraphiteSizes);
    addPopulations(new DraineSilicateGrainComposition(this), amin, amax, &dnda_sil, _numSilicateSizes);
}

//////////////////////////////////////////////////////////////////////
