/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LOGWAVELENGTHGRID_HPP
#define LOGWAVELENGTHGRID_HPP

#include "PanWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** LogWavelengthGrid is a subclass of the PanWavelengthGrid class representing logarithmically
    distributed wavelength grids. */
class LogWavelengthGrid : public PanWavelengthGrid
{
    ITEM_CONCRETE(LogWavelengthGrid, PanWavelengthGrid, "a logarithmic wavelength grid")

    PROPERTY_DOUBLE(minWavelength, "the shortest wavelength")
        ATTRIBUTE_QUANTITY(minWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(minWavelength, "1 A")
        ATTRIBUTE_MAX_VALUE(minWavelength, "1 m")

    PROPERTY_DOUBLE(maxWavelength, "the longest wavelength")
        ATTRIBUTE_QUANTITY(maxWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(maxWavelength, "1 A")
        ATTRIBUTE_MAX_VALUE(maxWavelength, "1 m")

    PROPERTY_INT(numWavelengths, "the number of wavelength grid points")
        ATTRIBUTE_MIN_VALUE(numWavelengths, "3")
        ATTRIBUTE_DEFAULT_VALUE(numWavelengths, "25")

    ITEM_END()

    /** \fn numWavelengths
        Implementation note: this function hides the WavelengthGrid::numWavelengths() function with
        the same name; this is not a problem because once the wavelength grid has been constructed,
        both functions return the same value. */

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs the vector with wavelengths. The \f$N\f$ wavelength grid points
        are distributed logarithmically between \f$\lambda_{\text{min}}\f$ and
        \f$\lambda_{\text{max}}\f$. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
