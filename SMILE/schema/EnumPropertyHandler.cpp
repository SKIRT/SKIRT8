/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EnumPropertyHandler.hpp"
#include "Item.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

bool EnumPropertyHandler::isTrueInCondition() const
{
    return property()->trueIf() == value();
}

////////////////////////////////////////////////////////////////////

bool EnumPropertyHandler::hasDefaultValue() const
{
    return isValid(property()->defaultValue());
}

////////////////////////////////////////////////////////////////////

void EnumPropertyHandler::acceptVisitor(PropertyHandlerVisitor* visitor)
{
    visitor->visitPropertyHandler(this);
}

////////////////////////////////////////////////////////////////////

vector<string> EnumPropertyHandler::values() const
{
    return property()->enumNames();
}

////////////////////////////////////////////////////////////////////

vector<string> EnumPropertyHandler::titlesForValues() const
{
    return property()->enumTitles();
}

////////////////////////////////////////////////////////////////////

string EnumPropertyHandler::defaultValue() const
{
    return property()->defaultValue();
}

////////////////////////////////////////////////////////////////////

string EnumPropertyHandler::value() const
{
    return target()->getEnumProperty(property());
}

////////////////////////////////////////////////////////////////////

string EnumPropertyHandler::titleForValue() const
{
    return property()->enumTitle(value());
}

////////////////////////////////////////////////////////////////////

bool EnumPropertyHandler::isValid(string value) const
{
    return StringUtils::contains(property()->enumNames(), value);
}

////////////////////////////////////////////////////////////////////

void EnumPropertyHandler::setValue(string value)
{
    target()->setEnumProperty(property(), value);
    setChanged();
}

////////////////////////////////////////////////////////////////////
