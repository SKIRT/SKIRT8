/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WEINGARTNERDRAINEDUSTMIX_HPP
#define WEINGARTNERDRAINEDUSTMIX_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The WeingartnerDraineDustMix class is a subclass of the MultiGrainDustMix class and represents
    realistic dust mixtures consisting of populations of graphite, silicate and PAH dust grains.
    The size distributions of each of these grains are fitted in such a way that they can reproduce
    the extinction curve of the Milky Way or the LMC. For details we refer to Weingartner & Draine
    (2001, ApJ, 548, 296). The graphite, silicate, and PAH populations (for both neutral and
    ionized PAHs) can be subdivided into \f$N_{\text{gra}}\f$, \f$N_{\text{sil}}\f$, and
    \f$N_{\text{PAH}}\f$ subpopulations, each corresponding to a distinct grain size bin. */
class WeingartnerDraineDustMix : public MultiGrainDustMix
{
    /** The enumeration type indicating the environment for the Weingartner-Draine dust. */
    ENUM_DEF(Environment, MilkyWay, LMC)
    ENUM_VAL(Environment, MilkyWay, "the Milky Way")
    ENUM_VAL(Environment, LMC, "the Large Magellanic Cloud")
    ENUM_END()

    ITEM_CONCRETE(WeingartnerDraineDustMix, MultiGrainDustMix, "a Weingartner & Draine multi-component dust mix")

    PROPERTY_INT(numGraphiteSizes, "the number of graphite grain size bins")
        ATTRIBUTE_MIN_VALUE(numGraphiteSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numGraphiteSizes, "5")

    PROPERTY_INT(numSilicateSizes, "the number of silicate grain size bins")
        ATTRIBUTE_MIN_VALUE(numSilicateSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numSilicateSizes, "5")

    PROPERTY_INT(numPAHSizes, "the number of neutral and ionized PAH size bins (each)")
        ATTRIBUTE_MIN_VALUE(numPAHSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numPAHSizes, "5")

    PROPERTY_ENUM(environment, Environment, "the typical environment for the dust")
        ATTRIBUTE_DEFAULT_VALUE(environment, "MilkyWay")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function adds the requested number of dust populations based on the
        DraineGraphiteGrainComposition, DraineSilicateGrainComposition,
        DraineNeutralPAHGrainComposition, and DraineIonizedPAHGrainComposition grain composition
        classes. It is assumed that 50% of the PAH grains are neutral and 50% are ionized. The
        grain size distributions for the various populations are given as complicated analytical
        formulas with predefined constant parameters depending on the selected environment (Milky
        Way or LMC). For the graphite and silicate populations, the distribution is a power-law
        function with a curvature and an exponential cutoff. For the PAH populations, it is the sum
        of two log-normal distributions, cut off at some upper grain size. This represents the
        common distribution for neutral and ionized PAH grains; to obtain the distribution for one
        of both subpopulations, we need to divide by two (assuming an ionization fraction of 50%).
        The exact forms of the grain size distribution functions can be found in Weingartner &
        Draine (2001, ApJ, 548, 296). */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
