/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OLIGOFITSCHEME_HPP
#define OLIGOFITSCHEME_HPP

#include "FitScheme.hpp"
#include "AdjustableSkirtSimulation.hpp"
#include "Optimization.hpp"
#include "ParameterRanges.hpp"
#include "ReferenceImages.hpp"

////////////////////////////////////////////////////////////////////

/** OligoFitScheme represents a complete FitSKIRT fit scheme for oligochromatic fits/simulations.
    It contains the SKIRT simulation, parameter ranges, reference images and the optimization type.*/
class OligoFitScheme : public FitScheme
{
    ITEM_CONCRETE(OligoFitScheme, FitScheme, "an oligochromatic fit scheme")

    PROPERTY_ITEM(simulation, AdjustableSkirtSimulation, "the SKIRT simulation to be run for this fit scheme")

    PROPERTY_ITEM(parameterRanges, ParameterRanges, "the parameter ranges")

    PROPERTY_ITEM(referenceImages, ReferenceImages, "the reference images")

    PROPERTY_ITEM(optimization, Optimization, "the optimization settings")

    ITEM_END()

    //======== Construction - Setup - Run - Destruction  ===========

protected:
    /** This function actually runs the fit scheme. It assumes that setup() has been already
        performed. */
    void runSelf() override;
};

////////////////////////////////////////////////////////////////////

#endif
