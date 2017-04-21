/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParameterRange.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void ParameterRange::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    if (_maxValue < _minValue)
        throw FATALERROR("Maximum value of parameter range cannot be smaller than minimum value");
}

////////////////////////////////////////////////////////////////////

string ParameterRange::quantityString() const
{
    switch (_quantityType)
    {
    case PhysicalQuantity::dimless: return string();
    case PhysicalQuantity::length: return "length";
    case PhysicalQuantity::distance: return "distance";
    case PhysicalQuantity::mass: return "mass";
    case PhysicalQuantity::posangle: return "posangle";
    }
    return string();
}

////////////////////////////////////////////////////////////////////
