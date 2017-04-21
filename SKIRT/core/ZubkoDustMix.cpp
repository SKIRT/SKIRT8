/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ZubkoDustMix.hpp"
#include "DraineGraphiteGrainComposition.hpp"
#include "DraineIonizedPAHGrainComposition.hpp"
#include "DraineNeutralPAHGrainComposition.hpp"
#include "DraineSilicateGrainComposition.hpp"
#include "ZubkoGraphiteGrainSizeDistribution.hpp"
#include "ZubkoPAHGrainSizeDistribution.hpp"
#include "ZubkoSilicateGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

void ZubkoDustMix::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    addPopulations(new DraineGraphiteGrainComposition(this), new ZubkoGraphiteGrainSizeDistribution(this),
                   _numGraphiteSizes);
    addPopulations(new DraineSilicateGrainComposition(this), new ZubkoSilicateGrainSizeDistribution(this),
                   _numSilicateSizes);
    addPopulations(new DraineNeutralPAHGrainComposition(this), new ZubkoPAHGrainSizeDistribution(this, 0.5),
                   _numPAHSizes);
    addPopulations(new DraineIonizedPAHGrainComposition(this), new ZubkoPAHGrainSizeDistribution(this, 0.5),
                   _numPAHSizes);
}

////////////////////////////////////////////////////////////////////
