/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PERSPECTIVEINSTRUMENT_HPP
#define PERSPECTIVEINSTRUMENT_HPP

#include "Instrument.hpp"
#include "HomogeneousTransform.hpp"
#include "ParallelDataCube.hpp"

////////////////////////////////////////////////////////////////////

/** The PerspectiveInstrument class implements a full perspective view of the simulated model, with
    arbitrary placement of the viewport outside or inside the model. For each wavelength the
    instrument maintains the total luminosity per pixel for all photon packages arriving from the
    front; light emitted behind the viewport is ignored. The instrument does \em not maintain the
    luminosity fractions caused by various phenomena (such as scattered versus direct light), nor
    does it keep track of the integrated luminosity across the viewport.

    The perspective instrument is intended mostly for making movies. Each movie frame is generated
    by a separate perspective instrument with the appropriate parameters.
*/
class PerspectiveInstrument : public Instrument
{
    ITEM_CONCRETE(PerspectiveInstrument, Instrument, "a perspective instrument (mostly for making movies)")

    PROPERTY_INT(numPixelsX, "the number of viewport pixels in the horizontal direction")
        ATTRIBUTE_MIN_VALUE(numPixelsX, "25")
        ATTRIBUTE_MAX_VALUE(numPixelsX, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsX, "250")

    PROPERTY_INT(numPixelsY, "the number of viewport pixels in the vertical direction")
        ATTRIBUTE_MIN_VALUE(numPixelsY, "25")
        ATTRIBUTE_MAX_VALUE(numPixelsY, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsY, "250")

    PROPERTY_DOUBLE(width, "the width of the viewport")
        ATTRIBUTE_QUANTITY(width, "length")
        ATTRIBUTE_MIN_VALUE(width, "]0")

    PROPERTY_DOUBLE(viewX, "the position of the viewport origin, x component")
        ATTRIBUTE_QUANTITY(viewX, "length")

    PROPERTY_DOUBLE(viewY, "the position of the viewport origin, y component")
        ATTRIBUTE_QUANTITY(viewY, "length")

    PROPERTY_DOUBLE(viewZ, "the position of the viewport origin, z component")
        ATTRIBUTE_QUANTITY(viewZ, "length")

    PROPERTY_DOUBLE(crossX, "the position of the crosshair, x component")
        ATTRIBUTE_QUANTITY(crossX, "length")

    PROPERTY_DOUBLE(crossY, "the position of the crosshair, y component")
        ATTRIBUTE_QUANTITY(crossY, "length")

    PROPERTY_DOUBLE(crossZ, "the position of the crosshair, z component")
        ATTRIBUTE_QUANTITY(crossZ, "length")

    PROPERTY_DOUBLE(upX, "the upwards direction, x component")
        ATTRIBUTE_QUANTITY(upX, "length")

    PROPERTY_DOUBLE(upY, "the upwards direction, y component")
        ATTRIBUTE_QUANTITY(upY, "length")

    PROPERTY_DOUBLE(upZ, "the upwards direction, z component")
        ATTRIBUTE_QUANTITY(upZ, "length")

    PROPERTY_DOUBLE(focal, "the distance from the eye to the viewport origin")
        ATTRIBUTE_QUANTITY(focal, "length")
        ATTRIBUTE_MIN_VALUE(focal, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that all attribute values have been appropriately set and performs
        setup for the instrument. Its most important task is to determine the appropriate
        perspective transformation. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** Returns the direction towards the eye from the given photon package launching position. */
    Direction bfkobs(const Position& bfr) const override;

    /** Returns the direction along the positive x-axis of the instrument frame, expressed in model
        coordinates. */
    Direction bfkx() const override;

    /** Returns the direction along the positive y-axis of the instrument frame, expressed in model
        coordinates. */
    Direction bfky() const override;

protected:
    /** This function simulates the detection of a photon package by the instrument. */
    void detect(PhotonPackage* pp) override;

    /** This function calibrates and outputs the instrument data. */
    void write() override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const int& _Nx{_numPixelsX};
    const int& _Ny{_numPixelsY};
    const double& _Sx{_width};
    const double& _Vx{_viewX};
    const double& _Vy{_viewY};
    const double& _Vz{_viewZ};
    const double& _Cx{_crossX};
    const double& _Cy{_crossY};
    const double& _Cz{_crossZ};
    const double& _Ux{_upX};
    const double& _Uy{_upY};
    const double& _Uz{_upZ};
    const double& _Fe{_focal};

    // data members derived from the published attributes during setup
    double _s{0.};          // width and height of a pixel
    double _Ex{0.}, _Ey{0.}, _Ez{0.};  // eye position
    Direction _bfkx;        // unit vector along the viewport's x-axis
    Direction _bfky;        // unit vector along the viewport's y-axis
    HomogeneousTransform _transform;   // transform from world to pixel coordinates

    // data cube
    ParallelDataCube _ftotv;
};

////////////////////////////////////////////////////////////////////

#endif
