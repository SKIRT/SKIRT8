/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEMISSIVITY_HPP
#define DUSTEMISSIVITY_HPP

#include "SimulationItem.hpp"
#include "Array.hpp"
class DustMix;

//////////////////////////////////////////////////////////////////////

/** The DustEmissivity class is the abstract base class for objects that calculate the wavelength
    dependent emissivity of a particular dust mix in a given radiation field. The emissivity can be
    determined from the optical properties of the dust mixture and the interstellar radiation field
    (both specified as arguments by the caller), and some additional assumptions. DustEmissivity
    subclasses implement various assumptions, and in particular, whether transient heating is taken
    into account or not. */
class DustEmissivity : public SimulationItem
{
    ITEM_ABSTRACT(DustEmissivity, SimulationItem, "a dust emissivity type")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dust emissivity \f$\varepsilon_\ell\f$ at all wavelength indices
        \f$\ell\f$ for a dust mix of the specified type residing in the specified mean radiation
        field \f$J_\ell\f$, assuming the simulation's wavelength grid. */
    virtual Array emissivity(const DustMix* mix, const Array& Jv) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
