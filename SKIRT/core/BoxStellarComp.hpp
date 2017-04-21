/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOXSTELLARCOMP_HPP
#define BOXSTELLARCOMP_HPP

#include "StellarComp.hpp"
#include "Box.hpp"

//////////////////////////////////////////////////////////////////////

/** The BoxStellarComp class is an abstract subclass of the StellarComp class, and represents any
    stellar component defined within a cuboidal volume with faces that are aligned with the planes
    of the coordinate system (a box). The class also inherits from the Box class. */
class BoxStellarComp : public StellarComp, public Box
{
    ITEM_ABSTRACT(BoxStellarComp, StellarComp, "a stellar component bounded by a box")

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
    /** This function returns the dimension of the stellar component, which for this class is always 3
        since there are no symmetries in the geometry. */
    int dimension() const override;
};

//////////////////////////////////////////////////////////////////////

#endif
