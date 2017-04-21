/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ConfigurableDustMix.hpp"

//////////////////////////////////////////////////////////////////////

void ConfigurableDustMix::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    // add the dust populations to the dust mix
    for (DustMixPopulation* population : _populations)
    {
        population->setup();    // since we're in setupSelfBefore, our children aren't yet setup
        addPopulations(population->composition(), population->sizeDistribution(), population->numGrainSizes());
    }
}

//////////////////////////////////////////////////////////////////////
