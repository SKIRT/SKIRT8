/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PANWAVELENGTHGRID_HPP
#define PANWAVELENGTHGRID_HPP

#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** PanWavelengthGrid is an abstract subclass of the general WavelengthGrid class representing
    wavelength grids for use in panchromatic simulations. It calculates appropriate wavelength bin
    widths based on the wavelength vector setup by a subclass in its setupSelfAfter() function. */
class PanWavelengthGrid : public WavelengthGrid
{
    ITEM_ABSTRACT(PanWavelengthGrid, WavelengthGrid, "a wavelength grid for use in a panchromatic simulation")

    PROPERTY_BOOL(writeWavelengths, "output a data file listing the wavelength grid points")
        ATTRIBUTE_DEFAULT_VALUE(writeWavelengths, "false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets the wavelength grid to the central wavelengths widths given as an
        argument, and calculates and sets the corresponding bin widths. It should be called from
        the setupSelfBefore() function in each PanWavelengthGrid subclass. Because the wavelength
        grid in a panchromatic simulation represents a sampled wavelength range, there must be at
        least three wavelengths in the specified array. If this requirement is not met, the
        function throws a fatal error.

        This function calculates the wavelength bin widths \f$\Delta_\ell\f$ for each of the
        wavelength points in the wavelength grid. The border of each bin is determined as the
        geometric mean of the wavelengths bordering it. In concreto, the width of the \f$\ell\f$'th
        bin is \f[ \Delta_\ell = \sqrt{\lambda_\ell\lambda_{\ell+1}} -
        \sqrt{\lambda_{\ell-1}\lambda_\ell}. \f] For the first bin this becomes \f[ \Delta_0 =
        \sqrt{\lambda_0\lambda_1} - \lambda_0 = \sqrt{\lambda_{\text{min}}\lambda_1} -
        \lambda_{\text{min}} \f] while for the last bin we obtain \f[ \Delta_{N-1} = \lambda_{N-1}
        - \sqrt{\lambda_{N-2}\lambda_{N-1}} = \lambda_{\text{max}} -
        \sqrt{\lambda_{N-2}\lambda_{\text{max}}}. \f]

        If the \em writeWavelengths flag is turned on, this function also outputs a data file
        called <tt>prefix_wavelengths.dat</tt> listing the wavelength grid points and bin widths.
        Apart from the header comment line, the text file contains a line for each of the \f$N\f$
        wavelength grid points. Each line has two columns: the first column specifies the
        wavelength \f$\lambda_\ell\f$ and the second columns specifies the corresponding bin width
        \f$\Delta_\ell\f$. */
    void setWavelengths(const Array& lambdav);

    //======================== Other Functions =======================

public:
    /** This function's implementation for this class always returns true since a PanWavelengthGrid
        represents a sampled wavelength range, as required for panchromatic simulations. */
    bool isSampledRange() const override;
};

//////////////////////////////////////////////////////////////////////

#endif
