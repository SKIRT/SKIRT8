/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SIMPLEOLIGODUSTMIX_HPP
#define SIMPLEOLIGODUSTMIX_HPP

#include "DustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The SimpleOligoDustMix class represents, as its name indicates, a simple dust mixture that
    can be used for oligochromatic simulations. For each wavelength in the global wavelength grid, it just
    reads in user-provided values for the extinction coefficient, the albedo and scattering asymmetry parameter.  */
class SimpleOligoDustMix : public DustMix
{
    ITEM_CONCRETE(SimpleOligoDustMix, DustMix, "a simple oligochromatic dust mix")
        ATTRIBUTE_ALLOWED_IF(SimpleOligoDustMix, "OligoMonteCarloSimulation")

    PROPERTY_DOUBLE_LIST(opacities, "the extinction coefficients, one for each wavelength")
        ATTRIBUTE_QUANTITY(opacities, "opacity")
        ATTRIBUTE_MIN_VALUE(opacities, "[0")

    PROPERTY_DOUBLE_LIST(albedos, "the scattering albedos, one for each wavelength")
        ATTRIBUTE_MIN_VALUE(albedos, "[0")
        ATTRIBUTE_MAX_VALUE(albedos, "1]")

    PROPERTY_DOUBLE_LIST(asymmetryParameters, "the asymmetry parameters, one for each wavelength")
        ATTRIBUTE_MIN_VALUE(asymmetryParameters, "[-1")
        ATTRIBUTE_MAX_VALUE(asymmetryParameters, "1]")
        ATTRIBUTE_DEFAULT_VALUE(asymmetryParameters, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

private:
    /** This function directly calculates all dust mix properties on the simulation's wavelength
        grid. It then adds a single dust population to the dust mix. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
