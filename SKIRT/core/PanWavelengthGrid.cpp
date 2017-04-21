/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PanWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void PanWavelengthGrid::setWavelengths(const Array& lambdav)
{
    // a wavelength range should be sampled with at least 3 bins for integration algorithms not to crash
    int Nlambda = lambdav.size();
    if (Nlambda < 3) throw FATALERROR("There must be at least three bins in a panchromatic wavelength grid");

    // determine wavelength bin widths
    Array dlambdav(Nlambda);
    setWavelengthBins(lambdav, dlambdav);   // set the grid with temporary widths so that lambdamin/max work
    for (int ell=0; ell<Nlambda; ell++)
    {
        dlambdav[ell] = lambdamax(ell)-lambdamin(ell);
    }

    // if requested, write a data file with the wavelengths and bin widths
    if (_writeWavelengths)
    {
        Units* units = find<Units>();

        // Create a text file
        TextOutFile file(this, "wavelengths", "wavelengths");

        // Write the header
        file.addColumn("lambda (" + units->uwavelength() + ")", 'e', 8);
        file.addColumn("delta lambda (" + units->uwavelength() + ")", 'e', 8);

        // Write the body
        for (int ell=0; ell<Nlambda; ell++)
        {
            file.writeRow(vector<double>{ units->owavelength(lambdav[ell]), units->owavelength(dlambdav[ell]) });
        }
    }

    // set the wavelength grid
    setWavelengthBins(lambdav, dlambdav);
}

//////////////////////////////////////////////////////////////////////

bool PanWavelengthGrid::isSampledRange() const
{
    return true;
}

//////////////////////////////////////////////////////////////////////
