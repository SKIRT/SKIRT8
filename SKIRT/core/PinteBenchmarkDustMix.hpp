/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PINTEBENCHMARKDUSTMIX_HPP
#define PINTEBENCHMARKDUSTMIX_HPP

#include "DustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The PinteBenchmarkDustMix class represents the dust mixture used for the benchmark described by
    Pinte et al. 2009, A&A, 498, 967P. The dust mixture consists of a single population of dust
    grains with a size of 1 micron. The class supports scattering polarization assuming spherical
    grains, loading the Mueller matrix coefficients from the appropriate resource data file. */
class PinteBenchmarkDustMix : public DustMix
{
    ITEM_CONCRETE(PinteBenchmarkDustMix, DustMix, "a benchmark dust mix including polarization (Pinte et al. 2009)")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the optical properties from the appropriate resource file, and adds a
        single dust population with these properties to the dust mix. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
