/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTCOMPNORMALIZATION_HPP
#define DUSTCOMPNORMALIZATION_HPP

#include "SimulationItem.hpp"
class Geometry;
class DustMix;

//////////////////////////////////////////////////////////////////////

/** DustCompNormalization is an abstract class that describes the way in which the normalization of
    a dust component should be accomplished. */
class DustCompNormalization : public SimulationItem
{
    ITEM_ABSTRACT(DustCompNormalization, SimulationItem, "the normalization for a dust component")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This pure virtual function returns the appropriate normalization factor for the specified
        geometry and dust mixture. */
    virtual double normalizationFactor(const Geometry* geom, const DustMix* mix) const = 0;
};

//////////////////////////////////////////////////////////////////////

#endif
