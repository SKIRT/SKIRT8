/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MESHDUSTCOMPONENT_HPP
#define MESHDUSTCOMPONENT_HPP

#include "SimulationItem.hpp"
#include "DustMix.hpp"

//////////////////////////////////////////////////////////////////////

/** A MeshDustComponent instance represents a dust component in a dust distribution of class
    AdaptiveMeshDustDistribution or VoronoiDustDistribution. Its properties specify the column(s)
    in the mesh data file defining the dust density distribution for the component, and the dust
    mix to be used for the component.

    The property \em densityIndex specifies the index \f$g_d\;(0\le g_d \le N_{fields}-1)\f$ of the
    field that should be interpreted as a density distribution \f$D\f$ over the domain. If the
    property \em multiplierIndex is nonnegative (the default value is -1), it specifies the index
    \f$g_m\;(0\le g_m \le N_{fields}-1)\f$ for the field that will serve as a multiplication factor
    for the density field. Finally, the density is always multiplied by the fraction \f$f\f$
    specified by the property \em densityFraction (with a default value of 1). In other words the
    density field value for each cell is determined by \f$D=F_{g_d}\times F_{g_m}\times f\f$. */
class MeshDustComponent : public SimulationItem
{
    ITEM_CONCRETE(MeshDustComponent, SimulationItem, "a mesh dust component")

    PROPERTY_INT(densityIndex, "the index of the column defining the density distribution for the dust component")
        ATTRIBUTE_MIN_VALUE(densityIndex, "0")
        ATTRIBUTE_MAX_VALUE(densityIndex, "99")
        ATTRIBUTE_DEFAULT_VALUE(densityIndex, "0")

    PROPERTY_INT(multiplierIndex, "the index of the column defining an extra multiplication factor, or -1")
        ATTRIBUTE_MIN_VALUE(multiplierIndex, "-1")
        ATTRIBUTE_MAX_VALUE(multiplierIndex, "99")
        ATTRIBUTE_DEFAULT_VALUE(multiplierIndex, "-1")

    PROPERTY_DOUBLE(densityFraction, "the fraction of the density actually locked up in dust grains")
        ATTRIBUTE_MIN_VALUE(densityFraction, "]0")
        ATTRIBUTE_MAX_VALUE(densityFraction, "1]")
        ATTRIBUTE_DEFAULT_VALUE(densityFraction, "1")

    PROPERTY_ITEM(mix, DustMix, "the dust mixture for the dust component")
        ATTRIBUTE_DEFAULT_VALUE(mix, "InterstellarDustMix")

    ITEM_END()
};

//////////////////////////////////////////////////////////////////////

#endif
