/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ZUBKODUSTMIX_HPP
#define ZUBKODUSTMIX_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The ZubkoDustMix class is a subclass of the MultiGrainDustMix class and represents a realistic
    dust mixture of bare (i.e. non-composite) graphite, silicate, neutral PAH and ionized PAH dust
    grains. The size distribution of each of these dust grain populations is finetuned in such a
    way that the global dust properties accurately reproduce the extinction, emission and abundance
    constraints on the Milky Way. The size distributions are taken from Zubko, Dwek & Arendt (2004,
    ApJS, 152, 211), model BARE_GR_S. The graphite, silicate, neutral PAH and ionized PAH
    populations can be subdivided into \f$N_{\text{gra}}\f$, \f$N_{\text{sil}}\f$, and
    \f$N_{\text{PAH}}\f$ subpopulations, each corresponding to a distinct grain size bin. */
class ZubkoDustMix : public MultiGrainDustMix
{
    ITEM_CONCRETE(ZubkoDustMix, MultiGrainDustMix, "a Zubko et al. multi-component dust mix")

    PROPERTY_INT(numGraphiteSizes, "the number of graphite grain size bins")
        ATTRIBUTE_MIN_VALUE(numGraphiteSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numGraphiteSizes, "5")

    PROPERTY_INT(numSilicateSizes, "the number of silicate grain size bins")
        ATTRIBUTE_MIN_VALUE(numSilicateSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numSilicateSizes, "5")

    PROPERTY_INT(numPAHSizes, "the number of neutral and ionized PAH size bins (each)")
        ATTRIBUTE_MIN_VALUE(numPAHSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numPAHSizes, "5")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function adds the requested number of dust populations based on the
        DraineGraphiteGrainComposition, DraineSilicateGrainComposition,
        DraineNeutralPAHGrainComposition, and DraineIonizedPAHGrainComposition grain composition
        classes. It is assumed that 50% of the PAH grains are neutral and 50% are ionized. The
        grain size distributions for the various populations are given as a complicated analytical
        formula, parameterized differently depending on the grain composition, which can be found
        in Zubko, Arendt & Dwek (2004, ApJS, 152, 211). */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
