/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef THEMISDUSTMIX_HPP
#define THEMISDUSTMIX_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The ThemisDustMix class is a subclass of the MultiGrainDustMix class and represents the THEMIS
    model for dust in the diffuse interstellar medium described by Jones et al. (2015, A&A,
    submitted) and the many references therein. In this model, there are two families of dust
    particles: amorphous hydrocarbon and amorphous silicates. For the amorphous hydrocarbon dust,
    the size distribution is a combination of a lognormal and a power-law distribution. For the
    silicates, it is assumed that 50\% of the mass is amorphous enstatite, and that the remaining
    half is amorphous forsterite, and the size distribution is a lognormal distribution (the same
    for both). The three populations in the mixture (hydrocarbon dust, enstatite and forsterite)
    can be subdivided into \f$N_{\text{ahc}}\f$, \f$N_{\text{ens}}\f$ and \f$N_{\text{for}}\f$
    subpopulations, each corresponding to a distinct grain size bin. */
class ThemisDustMix : public MultiGrainDustMix
{
    ITEM_CONCRETE(ThemisDustMix, MultiGrainDustMix, "a THEMIS multi-component dust mix")

    PROPERTY_INT(numHydrocarbonSizes, "the number of amorphous hydrocarbon grain size bins")
        ATTRIBUTE_MIN_VALUE(numHydrocarbonSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numHydrocarbonSizes, "5")

    PROPERTY_INT(numEnstatiteSizes, "the number of enstatite grain size bins")
        ATTRIBUTE_MIN_VALUE(numEnstatiteSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numEnstatiteSizes, "5")

    PROPERTY_INT(numForsteriteSizes, "the number of forsterite grain size bins")
        ATTRIBUTE_MIN_VALUE(numForsteriteSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numForsteriteSizes, "5")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function adds the requested number of dust populations based on the
        AmHydrocarbonGrainComposition, EnstatiteGrainComposition and ForsteriteGrainComposition
        grain composition classes, and on the appropriate grain size distributions. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif

