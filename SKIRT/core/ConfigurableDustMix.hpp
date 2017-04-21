/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONFIGURABLEDUSTMIX_HPP
#define CONFIGURABLEDUSTMIX_HPP

#include "MultiGrainDustMix.hpp"
#include "DustMixPopulation.hpp"

////////////////////////////////////////////////////////////////////

/** The ConfigurableDustMix class is a subclass of the MultiGrainDustMix class and represents dust
    mixtures consisting of one or more dust populations, fully configurable through its attributes.
    Specifically, the class maintains a list DustMixPopulation instances, each of which represents
    a particular dust population with configurable grain composition and grain size distribution.
    */
class ConfigurableDustMix : public MultiGrainDustMix
{
    ITEM_CONCRETE(ConfigurableDustMix, MultiGrainDustMix, "a configurable multi-component dust mix")

    PROPERTY_ITEM_LIST(populations, DustMixPopulation, "the dust populations")
        ATTRIBUTE_DEFAULT_VALUE(populations, "DustMixPopulation")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function adds the configured dust populations to the dust mix. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
