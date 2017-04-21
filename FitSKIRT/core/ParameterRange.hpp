/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARAMETERRANGE_HPP
#define PARAMETERRANGE_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The ParameterRange class represents a numeric parameter range, including a label,
     the type of physical quantity, and a minimum and a maximum value. */
class ParameterRange : public SimulationItem
{
    /** The enumeration type indicating the type of physical quantity represented by this parameter
        range. Implementation note: the enumeration identifiers must match quantity strings defined
        in the SkirtUnitDef class (with the exception of the identifier for dimensionless
        quantities). */
    ENUM_DEF(PhysicalQuantity, dimless, length, distance, mass, posangle)
    ENUM_VAL(PhysicalQuantity, dimless, "dimensionless")
    ENUM_VAL(PhysicalQuantity, length, "length")
    ENUM_VAL(PhysicalQuantity, distance, "distance")
    ENUM_VAL(PhysicalQuantity, mass, "mass")
    ENUM_VAL(PhysicalQuantity, posangle, "position angle")
    ENUM_END()

    ITEM_CONCRETE(ParameterRange, SimulationItem, "a parameter range")

    PROPERTY_STRING(label, "the label identifying this parameter range")

    PROPERTY_ENUM(quantityType, PhysicalQuantity, "the type of physical quantity represented by this parameter range")
        ATTRIBUTE_DEFAULT_VALUE(quantityType, "length")

    PROPERTY_DOUBLE(minValue, "the minimum value in this parameter range")
        ATTRIBUTE_QUANTITY(minimumValue, "@quantityType")

    PROPERTY_DOUBLE(maxValue, "the maximum value in this parameter range")
        ATTRIBUTE_QUANTITY(maximumValue, "@quantityType")

    ITEM_END()

    //============ Construction - Setup - Destruction  =============

protected:
    /** This function verifies that the maximum value is larger than the minimum value. */
    void setupSelfBefore() override;

    //====================== Other functions =======================

public:
    /** Returns a string value indicating the physical quantity represented by this parameter
        range. This is one of the quantity strings defined in the SkirtUnitDef class, or the empty
        string for dimensionless quantities. */
    string quantityString() const;
};

////////////////////////////////////////////////////////////////////

#endif
