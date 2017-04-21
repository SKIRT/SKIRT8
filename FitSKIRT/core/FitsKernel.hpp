/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FITSKERNEL_HPP
#define FITSKERNEL_HPP

#include "ConvolutionKernel.hpp"

////////////////////////////////////////////////////////////////////

/** The FitsKernel class represents a convolution kernel that is imported from a FITS file. */
class FitsKernel : public ConvolutionKernel
{
    ITEM_CONCRETE(FitsKernel, ConvolutionKernel, "a convolution kernel read from a FITS file")

    PROPERTY_STRING(filename, "the name of the input image file")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function makes sure that the FITS file is imported. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
