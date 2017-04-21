/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MRNDUSTMIX_HPP
#define MRNDUSTMIX_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The MRNDustMix class is a subclass of the MultiGrainDustMix class and represents dust mixtures
    consisting of separate populations of graphite and silicate dust grains with a grain size
    distribution according to the famous MRN distribution (Mathis, Rumpl & Nordsieck 1977, ApJ,
    217, 425). The actual size distributions are taken from Weingartner & Draine (2001, ApJ, 548,
    296), the optical properties are taken from Bruce Draine's web site. The graphite and silicate
    populations can be subdivided into \f$N_{\text{gra}}\f$ and \f$N_{\text{sil}}\f$
    subpopulations, each corresponding to a distinct grain size bin. */
class MRNDustMix : public MultiGrainDustMix
{
    ITEM_CONCRETE(MRNDustMix, MultiGrainDustMix, "an MRN multi-component dust mix")

    PROPERTY_INT(numGraphiteSizes, "the number of graphite grain size bins")
        ATTRIBUTE_MIN_VALUE(numGraphiteSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numGraphiteSizes, "5")

    PROPERTY_INT(numSilicateSizes, "the number of silicate grain size bins")
        ATTRIBUTE_MIN_VALUE(numSilicateSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numSilicateSizes, "5")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function adds the requested number of dust populations based on the
        DraineGraphiteGrainComposition and DraineSilicateGrainComposition grain composition classes
        for graphite and silicate respectively, and on grain size distributions given by \f[
        \frac{\text{d}n}{\text{d}a} = C\, a^{-3.5} \qquad \text{for}\quad a_\text{min} \leq a \leq
        a_\text{max}, \f] with \f$C=10^{-25.13}\,\text{cm}^{2.5}\f$ for graphite and
        \f$C=10^{-25.11}\,\text{cm}^{2.5}\f$ for silicate, and with
        \f$a_\text{min}=50\,\text{\AA}\f$ and \f$a_\text{max}=0.25\,\mu\text{m}\f$ for both grain
        types. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
