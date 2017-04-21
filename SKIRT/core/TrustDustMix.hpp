/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TRUSTDUSTMIX_HPP
#define TRUSTDUSTMIX_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The TrustDustMix class is a subclass of the MultiGrainDustMix class and represents a realistic
    dust mixture of bare (i.e. non-composite) graphite, silicate, and PAH dust grains according to
    the dust model used for the TRUST benchmark simulations. The underlying data is provided by
    Karel Misselt as part of a download from the TRUST web site
    (http://ipag.osug.fr/RT13/RTTRUST/opa.php) describing the BARE-GR-S model of Zubko, Dwek, and
    Arendt 2004, ApJS, 152, 211.

    The graphite, silicate, and PAH populations can be subdivided into \f$N_{\text{gra}}\f$,
    \f$N_{\text{sil}}\f$, and \f$N_{\text{PAH}}\f$ subpopulations, each corresponding to a distinct
    grain size bin. */
class TrustDustMix : public MultiGrainDustMix
{
    ITEM_CONCRETE(TrustDustMix, MultiGrainDustMix, "a multi-component dust mix from the TRUST benchmark")

    PROPERTY_INT(numGraphiteSizes, "the number of graphite grain size bins")
        ATTRIBUTE_MIN_VALUE(numGraphiteSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numGraphiteSizes, "5")

    PROPERTY_INT(numSilicateSizes, "the number of silicate grain size bins")
        ATTRIBUTE_MIN_VALUE(numSilicateSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numSilicateSizes, "5")

    PROPERTY_INT(numPAHSizes, "the number of PAH size bins")
        ATTRIBUTE_MIN_VALUE(numPAHSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numPAHSizes, "5")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function adds the requested number of dust populations based on the
        TrustGraphiteGrainComposition, TrustSilicateGrainComposition, and
        TrustNeutralPAHGrainComposition grain composition classes. The grain size distributions for
        the various populations are given as a complicated analytical formula, parameterized
        differently depending on the grain composition, which can be found in Zubko, Arendt & Dwek
        (2004, ApJS, 152, 211). */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
