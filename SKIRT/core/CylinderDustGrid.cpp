/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CylinderDustGrid.hpp"
#include "FatalError.hpp"

//////////////////////////////////////////////////////////////////////

void CylinderDustGrid::setupSelfBefore()
{
    DustGrid::setupSelfBefore();
    if (_maxZ <= _minZ) throw FATALERROR("The extent of the cylinder should be positive in the Z direction");
}

//////////////////////////////////////////////////////////////////////

Box CylinderDustGrid::boundingBox() const
{
    return Box(-_maxRadius,-_maxRadius,_minZ, _maxRadius,_maxRadius,_maxZ);
}

//////////////////////////////////////////////////////////////////////