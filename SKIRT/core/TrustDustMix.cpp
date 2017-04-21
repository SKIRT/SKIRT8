/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TrustDustMix.hpp"
#include "TrustGraphiteGrainComposition.hpp"
#include "TrustNeutralPAHGrainComposition.hpp"
#include "TrustSilicateGrainComposition.hpp"
#include "ZubkoGraphiteGrainSizeDistribution.hpp"
#include "ZubkoPAHGrainSizeDistribution.hpp"
#include "ZubkoSilicateGrainSizeDistribution.hpp"

//////////////////////////////////////////////////////////////////////

void TrustDustMix::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    addPopulations(new TrustGraphiteGrainComposition(this), new ZubkoGraphiteGrainSizeDistribution(this),
                   _numGraphiteSizes);
    addPopulations(new TrustSilicateGrainComposition(this), new ZubkoSilicateGrainSizeDistribution(this),
                   _numSilicateSizes);
    addPopulations(new TrustNeutralPAHGrainComposition(this), new ZubkoPAHGrainSizeDistribution(this), _numPAHSizes);
}

////////////////////////////////////////////////////////////////////
