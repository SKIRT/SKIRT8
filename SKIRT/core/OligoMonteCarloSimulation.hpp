/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OLIGOMONTECARLOSIMULATION_HPP
#define OLIGOMONTECARLOSIMULATION_HPP

#include "MonteCarloSimulation.hpp"
#include "OligoDustSystem.hpp"
#include "OligoWavelengthGrid.hpp"
#include "StellarSystem.hpp"

//////////////////////////////////////////////////////////////////////

/** This is a subclass of the general MonteCarloSimulation class representing an oligochromatic
    Monte Carlo simulation, i.e. operating at one or more distinct wavelengths rather than a
    discretized range of wavelengths. In such simulations there can be absorption and scattering by
    dust grains, but by definition there is no thermal dust emission. */
class OligoMonteCarloSimulation : public MonteCarloSimulation
{
    ITEM_CONCRETE(OligoMonteCarloSimulation, MonteCarloSimulation, "an oligochromatic Monte Carlo simulation")

    PROPERTY_ITEM(wavelengthGrid, OligoWavelengthGrid, "the wavelength grid")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthGrid, "OligoWavelengthGrid")

    PROPERTY_ITEM(stellarSystem, StellarSystem, "the stellar system")
        ATTRIBUTE_DEFAULT_VALUE(stellarSystem, "StellarSystem")

    PROPERTY_ITEM(dustSystem, OligoDustSystem, "the dust system")
        ATTRIBUTE_DEFAULT_VALUE(dustSystem, "OligoDustSystem")
        ATTRIBUTE_OPTIONAL(dustSystem)

    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function actually runs the simulation. For an oligochromatic simulation, this just
        includes the stellar emission phase (plus writing the results). */
    void runSelf() override;
};

////////////////////////////////////////////////////////////////////

#endif
