/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OLIGOWAVELENGTHGRID_HPP
#define OLIGOWAVELENGTHGRID_HPP

#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** OligoWavelengthGrid is a subclass of the general WavelengthGrid class representing one or more
    distinct wavelengths rather than a discretized wavelength range. It is intended for use with
    oligochromatic simulations, which don't calculate the dust temperature by integrating over a
    wavelength range. */
class OligoWavelengthGrid : public WavelengthGrid
{
    ITEM_CONCRETE(OligoWavelengthGrid, WavelengthGrid, "a list of one or more distinct wavelengths")

    PROPERTY_DOUBLE_LIST(wavelengths, "the wavelengths")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets the wavelength grid to the central wavelengths widths specified by the
        \em wavelengths property, and sets artificial bin widths. Since an oligochromatic
        simulation consists of multiple monochromatic simulations for distinct wavelengths, the
        bins are taken to be very small and independent from the other wavelengths. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function's implementation for this class always returns false, since an
        OligoWavelengthGrid contains individual distinct wavelengths for use by oligochromatic
        simulations. */
    bool isSampledRange() const override;
};

//////////////////////////////////////////////////////////////////////

#endif
