/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BoxDustGrid.hpp"
#include "FatalError.hpp"
#include "Log.hpp"

//////////////////////////////////////////////////////////////////////

void BoxDustGrid::setupSelfBefore()
{
    DustGrid::setupSelfBefore();

    // copy the configured values into our inherited Box
    _xmin = _minX; _ymin = _minY; _zmin = _minZ;
    _xmax = _maxX; _ymax = _maxY; _zmax = _maxZ;

    // verify the values
    if (_xmax <= _xmin) throw FATALERROR("The extent of the box should be positive in the X direction");
    if (_ymax <= _ymin) throw FATALERROR("The extent of the box should be positive in the Y direction");
    if (_zmax <= _zmin) throw FATALERROR("The extent of the box should be positive in the Z direction");
}

//////////////////////////////////////////////////////////////////////

int BoxDustGrid::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

Box BoxDustGrid::boundingBox() const
{
    // use our inherited Box
    return extent();
}

//////////////////////////////////////////////////////////////////////
