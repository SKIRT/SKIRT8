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
    const double amax = 4.9e-6;

    // parameters for the size distribution of the silicates
    const double amins = 0.001e-6;
    const double a0s = 0.008e-6;
    const double Cs = 3.74461585691147e-10;

    // parameters for the lognormal part of the size distribution for the hydrocarbons
    const double aminh1 = 0.0005e-6;
    const double a0h = 0.007e-6;
    const double Ch1 = 1.88994983508769e-10;

    // parameters for the power-law part of the size distribution for the hydrocarbons
    const double aminh2 = 0.0004e-6;
    const double alpha = -5.0;
    const double at = 0.01e-6;
    const double ac = 0.05e-6;
    const double Ch2 = 1.59726058957057e-41;
    const double aminh = min(aminh1,aminh2);

    // amorphous hydrocarbons
    double dndah(double a)
    {
        if (a > amax) return 0.0;
        double dnda = 0.0;
        if (a > aminh1) dnda += Ch1/a * exp(-0.5*pow(log(a/a0h),2));
        if (a > aminh2) dnda += Ch2 * pow(a,alpha) * ( (a<=at) ? 1.0 : exp(-(a-at)/ac) );
        return dnda;
    }

    // amorphous silicates (the same distribution for enstatite and forsterite)
    double dndas(double a)
    {
        if (a < amins || a > amax) return 0.0;
        double x = log(a/a0s);
        return Cs/a * exp(-0.5*x*x);
    }
}

//////////////////////////////////////////////////////////////////////

void ThemisDustMix::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    AmHydrocarbonGrainComposition* ahcgc = new AmHydrocarbonGrainComposition(this);
    EnstatiteGrainComposition* ensgc = new EnstatiteGrainComposition(this,
                                                            EnstatiteGrainComposition::GrainType::Amorphous);
    ForsteriteGrainComposition* forgc = new ForsteriteGrainComposition(this,
                                                            ForsteriteGrainComposition::GrainType::Amorphous);

    addPopulations(ahcgc, aminh, amax, &dndah, _numHydrocarbonSizes);
    addPopulations(ensgc, amins, amax, &dndas, _numEnstatiteSizes);
    addPopulations(forgc, amins, amax, &dndas, _numForsteriteSizes);
}

////////////////////////////////////////////////////////////////////
