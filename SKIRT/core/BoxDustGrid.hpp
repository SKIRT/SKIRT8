/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOXDUSTGRID_HPP
#define BOXDUSTGRID_HPP

#include "DustGrid.hpp"

////////////////////////////////////////////////////////////////////

/** The BoxDustGrid class is an abstract subclass of the general DustGrid class, and represents
    any dust grid defined within cuboidal configuration space in which the faces are aligned with
    the planes of the coordinate system (a box). The class also inherits from the Box class. */
class BoxDustGrid : public DustGrid, public Box
{
    ITEM_ABSTRACT(BoxDustGrid, DustGrid, "a dust grid bounded by a box")

    PROPERTY_DOUBLE(minX, "the start point of the box in the X direction")
        ATTRIBUTE_QUANTITY(minX, "length")

    PROPERTY_DOUBLE(maxX, "the end point of the box in the X direction")
        ATTRIBUTE_QUANTITY(maxX, "length")

    PROPERTY_DOUBLE(minY, "the start point of the box in the Y direction")
        ATTRIBUTE_QUANTITY(minY, "length")

    PROPERTY_DOUBLE(maxY, "the end point of the box in the Y direction")
        ATTRIBUTE_QUANTITY(maxY, "length")

    PROPERTY_DOUBLE(minZ, "the start point of the box in the Z direction")
        ATTRIBUTE_QUANTITY(minZ, "length")

    PROPERTY_DOUBLE(maxZ, "the end point of the box in the Z direction")
        ATTRIBUTE_QUANTITY(maxZ, "length")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the box has a positive volume. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is 3 for all subclasses of this class. */
    int dimension() const override;

    /** This function returns the bounding box that encloses the dust grid. */
    Box boundingBox() const override;
};

////////////////////////////////////////////////////////////////////

#endif
