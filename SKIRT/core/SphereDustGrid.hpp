/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHEREDUSTGRID_HPP
#define SPHEREDUSTGRID_HPP

#include "DustGrid.hpp"

////////////////////////////////////////////////////////////////////

/** The SphereDustGrid class is an abstract subclass of the general DustGrid class, and represents
    any dust grid defined within a spherical configuration space, centered on the origin of the
    system. */
class SphereDustGrid : public DustGrid
{
    ITEM_ABSTRACT(SphereDustGrid, DustGrid, "a dust grid bounded by a sphere")

    PROPERTY_DOUBLE(maxRadius, "the outer radius of the grid")
        ATTRIBUTE_QUANTITY(maxRadius, "length")
        ATTRIBUTE_MIN_VALUE(maxRadius, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the bounding box that encloses the dust grid. */
    Box boundingBox() const override;
};

////////////////////////////////////////////////////////////////////

#endif
