/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ELECTRONDUSTMIX_HPP
#define ELECTRONDUSTMIX_HPP

#include "DustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The ElectronDustMix class represents the optical properties of a population of electrons,
    including polarization. For our purposes here, we consider electrons to be a simplified version
    of dust grains. Electrons do not absorb nor emit photons. Photon scattering is wavelength
    independent and there is only one scattering cross section, the constant Thomson cross section
    \f$\sigma_\text{t}\f$. The Sxx values in the Mueller matrix depend just on \f$\cos\theta\f$ and
    are given by equation (C.7) of Wolf 2003 (Computer Physics Communications, 150, 99–115). */
class ElectronDustMix : public DustMix
{
    ITEM_CONCRETE(ElectronDustMix, DustMix, "a 'dust mix' representing a population of electrons")

    PROPERTY_BOOL(addCircularPolarization, "use a synthetic particle that adds circular polarization during scattering")
        ATTRIBUTE_DEFAULT_VALUE(addCircularPolarization, "false")
        ATTRIBUTE_SILENT(addCircularPolarization)

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function directly calculates all electron dust mix properties on the simulation's
        wavelength grid. It then adds a single dust population to the dust mix. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
