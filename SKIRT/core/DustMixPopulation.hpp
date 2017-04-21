/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTMIXPOPULATION_HPP
#define DUSTMIXPOPULATION_HPP

#include "SimulationItem.hpp"
#include "GrainComposition.hpp"
#include "GrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** DustMixPopulation is simple class that represents a single dust population for use with the
    ConfigurableDustMix class. It combines a dust grain composition type (an instance of a
    GrainComposition subclass), a dust grain size distribution (an instance of a
    GrainSizeDistribution subclass) to define a particular dust population. In addition, it
    provides the option to split the grain size distribution into \f$N_\text{bins}\f$ bins on a
    logarithmic scale, configuring a seperate dust population for each bin. For more information
    see the MultiGrainDustMix::addpopulations() function. */
class DustMixPopulation : public SimulationItem
{
    ITEM_CONCRETE(DustMixPopulation, SimulationItem, "a dust mix population")

    PROPERTY_ITEM(composition, GrainComposition, "the dust grain composition")
        ATTRIBUTE_DEFAULT_VALUE(composition, "DraineGraphiteGrainComposition")

    PROPERTY_ITEM(sizeDistribution, GrainSizeDistribution, "the dust grain size distribution")
        ATTRIBUTE_DEFAULT_VALUE(sizeDistribution, "PowerLawGrainSizeDistribution")

    PROPERTY_INT(numGrainSizes, "the number of grain size bins")
        ATTRIBUTE_MIN_VALUE(numGrainSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numGrainSizes, "5")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
