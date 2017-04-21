/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STELLARCOMPNORMALIZATION_HPP
#define STELLARCOMPNORMALIZATION_HPP

#include "SimulationItem.hpp"
class SED;

////////////////////////////////////////////////////////////////////

/** StellarCompNormalization is a class that describes the way in which the normalization of a
    stellar component should be accomplished. Normalization just means setting the luminosity of
    the stellar component (the spatial and relative spectral distributions are determined
    separately). It is an abstract class; the concrete subclasses define distinct ways of
    normalization. A stellar component can, for example, be normalized by specifying the total
    bolometric luminosity or the absolute magnitude in a given filter band. */
class StellarCompNormalization : public SimulationItem
{
    ITEM_ABSTRACT(StellarCompNormalization, SimulationItem, "the normalization for a stellar component")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This virtual function returns the bolometric luminosity of a (virtual) stellar
        component that would have a given %SED. */
    virtual double totalLuminosity(SED* sed) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
