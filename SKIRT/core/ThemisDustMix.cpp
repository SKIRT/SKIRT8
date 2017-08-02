/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ThemisDustMix.hpp"
#include "AmHydrocarbonGrainComposition.hpp"
#include "EnstatiteGrainComposition.hpp"
#include "ForsteriteGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    const double aminhpl = 0.0004e-6;
    const double aminhln = 0.0005e-6;
    const double amins = 0.001e-6;
    const double amax = 4.9e-6;

    // the power-law size distribution for the hydrocarbons
    double dndahpl(double a)
    {
        const double alpha = -5.0;
        const double at = 0.01e-6;
        const double ac = 0.05e-6;
        const double C = 1.71726298266e-43;

        if (a < aminhpl || a > amax) return 0.0;
        return C * pow(a,alpha) * ( (a<=at) ? 1.0 : exp(-(a-at)/ac) );
    }

    // the lognormal size distribution for the hydrocarbons
    double dndahln(double a)
    {
        const double a0 = 0.007e-6;
        const double C = 2.05052478683e-12;

        if (a < aminhln || a > amax) return 0.0;
        double x = log(a/a0);
        return C/a * exp(-0.5*x*x);
    }

    // lognormal size distribution for the silicates (the same distribution for enstatite and forsterite)
    double dndas(double a)
    {
        const double a0 = 0.008e-6;
        const double C = 4.02595019205e-12;

        if (a < amins || a > amax) return 0.0;
        double x = log(a/a0);
        return C/a * exp(-0.5*x*x);
    }
}

//////////////////////////////////////////////////////////////////////

void ThemisDustMix::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    AmHydrocarbonGrainComposition* gchpl = new AmHydrocarbonGrainComposition(this, 1600.);
    AmHydrocarbonGrainComposition* gchln = new AmHydrocarbonGrainComposition(this, 1510.);

    ForsteriteGrainComposition* gcfor = new ForsteriteGrainComposition(this,
                                                            ForsteriteGrainComposition::GrainType::Amorphous);
    EnstatiteGrainComposition* gcens = new EnstatiteGrainComposition(this,
                                                            EnstatiteGrainComposition::GrainType::Amorphous);

    addPopulations(gchpl, aminhpl, amax, &dndahpl, _numHydrocarbonSizes);
    addPopulations(gchln, aminhln, amax, &dndahln, _numHydrocarbonSizes);
    addPopulations(gcens, amins,   amax, &dndas,   _numSilicateSizes);
    addPopulations(gcfor, amins,   amax, &dndas,   _numSilicateSizes);
}

////////////////////////////////////////////////////////////////////
