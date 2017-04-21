/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARAMETERRANGES_HPP
#define PARAMETERRANGES_HPP

#include "ParameterRange.hpp"

////////////////////////////////////////////////////////////////////

/** The ParameterRanges class represents a complete set of parameter ranges.
    Objects of this class are essentially lists of pointers to ParameterRange
    objects. */
class ParameterRanges: public SimulationItem
{
    ITEM_CONCRETE(ParameterRanges, SimulationItem, "a list of parameter ranges")

    PROPERTY_ITEM_LIST(ranges, ParameterRange, "the parameter ranges")
        ATTRIBUTE_DEFAULT_VALUE(ranges, "ParameterRange")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
