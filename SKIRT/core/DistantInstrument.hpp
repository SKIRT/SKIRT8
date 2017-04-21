/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DISTANTINSTRUMENT_HPP
#define DISTANTINSTRUMENT_HPP

#include "Instrument.hpp"

////////////////////////////////////////////////////////////////////

/** DistantInstrument is an abstract class for representing instruments positioned at a
    sufficiently large distance from the system, so that the observable sky is a plane
    perpendicular to the line of sight rather than a part of a sphere. Consequently
    parallel projection is used and the distance is only important for flux calibration.

    The direction \f${\boldsymbol{k}}_{\text{obs}} = (\theta,\varphi)\f$ towards the instrument is
    specified through the inclination angle \f$\theta\f$ and the azimuth angle \f$\varphi\f$.
    Finally, the instrument can be rotated about the line of sight by specifying its position angle
    \f$\omega\f$. The table below lists some typical values for these angles, in degrees.

    <TABLE>
    <TR><TD><I>View</I></TD>     <TD>\f$\theta\f$</TD> <TD>\f$\varphi\f$</TD> <TD>\f$\omega\f$</TD> </TR>
    <TR><TD>XY-plane</TD>        <TD>0</TD>   <TD>0</TD>   <TD>90</TD>  </TR>
    <TR><TD>XZ-plane</TD>        <TD>90</TD>  <TD>-90</TD> <TD>0</TD>   </TR>
    <TR><TD>YZ-plane</TD>        <TD>90</TD>  <TD>0</TD>   <TD>0</TD>   </TR>
    <TR><TD>first octant</TD>    <TD>45</TD>  <TD>45</TD>  <TD>0</TD>   </TR>
    </TABLE>
    */
class DistantInstrument : public Instrument
{
    ITEM_ABSTRACT(DistantInstrument, Instrument, "a distant instrument")

    PROPERTY_DOUBLE(distance, "the distance to the system")
        ATTRIBUTE_QUANTITY(distance, "distance")
        ATTRIBUTE_MIN_VALUE(distance, "]0")

    PROPERTY_DOUBLE(inclination, "the inclination angle θ of the detector")
        ATTRIBUTE_QUANTITY(inclination, "posangle")
        ATTRIBUTE_MIN_VALUE(inclination, "0 deg")
        ATTRIBUTE_MAX_VALUE(inclination, "180 deg")
        ATTRIBUTE_DEFAULT_VALUE(inclination, "0 deg")

    PROPERTY_DOUBLE(azimuth, "the azimuth angle φ of the detector")
        ATTRIBUTE_QUANTITY(azimuth, "posangle")
        ATTRIBUTE_MIN_VALUE(azimuth, "-360 deg")
        ATTRIBUTE_MAX_VALUE(azimuth, "360 deg")
        ATTRIBUTE_DEFAULT_VALUE(azimuth, "0 deg")

    PROPERTY_DOUBLE(positionAngle, "the position angle ω of the detector")
        ATTRIBUTE_QUANTITY(positionAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(positionAngle, "-360 deg")
        ATTRIBUTE_MAX_VALUE(positionAngle, "360 deg")
        ATTRIBUTE_DEFAULT_VALUE(positionAngle, "0 deg")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates some frequently used values for the instrument. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** Returns the direction towards the observer, as seen from the origin of the coordinate
        system. The provided photon package's launching position is not used; it is considered to
        be very close to the coordinate origin from the observer's standpoint, since the distance
        is sufficiently large. */
    Direction bfkobs(const Position& bfr) const override;

    /** Returns the direction along the positive x-axis of the instrument frame, expressed in model
        coordinates. The function applies the inverse instrument transformation to the pixel
        frame's x-axis. */
    Direction bfkx() const override;

    /** Returns the direction along the positive y-axis of the instrument frame, expressed in model
        coordinates. The function applies the inverse instrument transformation to the pixel
        frame's y-axis. */
    Direction bfky() const override;

protected:
    /** This convenience function calibrates one or more integrated luminosity data vectors
        gathered by a DistantInstrument subclass, and outputs them as columns in a single SED text
        file. The incoming data is organized as a list of data arrays and a second list of
        corresponding human-readable names. Each array in the first list is a 1D vector containing
        a luminosity value per wavelength in the simulation's wavelength grid. The strings in
        the second list are used to identify the columns in the output file. The two lists must
        have the same size, and should not be empty. The generated SED text file is named for the
        instrument as in <tt>prefix_instrument_SED.dat</tt>. The first column lists the wavelength;
        subsequent columns provide the corresponding data point from one of the specified arrays.
        The calibration performed by this function takes care of the
        conversion from bolometric luminosity units to flux density units. Typical units for the
        quantities in the SED file are are \f$\text{W}\,\text{m}^{-2}\f$. The calibration is
        performed in-place in the arrays, so the incoming data is overwritten. */
    void calibrateAndWriteSEDs(const vector<Array*>& Farrays, const vector<string>& Fnames);

    //======================== Data Members ========================

    // data members derived from the published attributes during setup
protected:
    double _cosphi{0.}, _sinphi{0.};
    double _costheta{0.}, _sintheta{0.};
    double _cospa{0.}, _sinpa{0.};
private:
    Direction _bfkobs;
    Direction _bfkx;
    Direction _bfky;
};

////////////////////////////////////////////////////////////////////

#endif
