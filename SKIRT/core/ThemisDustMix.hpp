/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef THEMISDUSTMIX_HPP
#define THEMISDUSTMIX_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/**
The ThemisDustMix class is a subclass of the MultiGrainDustMix class and represents the THEMIS
model for dust in the diffuse interstellar medium described by Jones et al. 2017 (A&A, 602,
A46) and the references therein. In this model, there are two families of dust particles:
amorphous hydrocarbon and amorphous silicates. For the amorphous hydrocarbon dust, the size
distribution is a combination of a lognormal and a power-law distribution. For the silicates,
it is assumed that 50\% of the mass is amorphous enstatite, and that the remaining half is
amorphous forsterite, and the size distribution is a lognormal distribution (the same for
both).

The following is an extract from the GRAIN_J17.DAT file specifying the dust model for DustEM.

\verbatim
# grain type, nsize, type keywords, Mdust/MH, rho, amin, amax, alpha/a0 [, at, ac, gamma (ED)] [, au, zeta, eta (CV)]
# cgs units
 CM20   100 plaw-ed   0.170E-02  1.600E+00  0.40E-07   4900.0E-07    -5.00E-00   10.00E-07   50.0E-07  1.000E+00
 CM20   100 logn      0.600E-03  1.510E+00  0.50E-07   4900.0E-07     7.00E-07    1.00E+00
 aOLM5  100 logn      0.255E-02  2.190E+00  1.00E-07   4900.0E-07     8.00E-07    1.00E+00
 aPyM5  100 logn      0.255E-02  2.190E+00  1.00E-07   4900.0E-07     8.00E-07    1.00E+00
\endverbatim

The four populations in the mixture (two hydrocarbon populations and two silicate populations)
can be subdivided into a number of subpopulations, each corresponding to a distinct grain size
bin.
*/
class ThemisDustMix : public MultiGrainDustMix
{
    ITEM_CONCRETE(ThemisDustMix, MultiGrainDustMix, "a THEMIS (Jones et al. 2017) multi-component dust mix")

    PROPERTY_INT(numHydrocarbonSizes, "the number of grain size bins for each of the hydrocarbon components")
        ATTRIBUTE_MIN_VALUE(numHydrocarbonSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numHydrocarbonSizes, "5")

    PROPERTY_INT(numSilicateSizes, "the number of grain size bins for each of the silicate components")
        ATTRIBUTE_MIN_VALUE(numSilicateSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numSilicateSizes, "5")

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

