/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GAUSSIANKERNEL_HPP
#define GAUSSIANKERNEL_HPP

#include "ConvolutionKernel.hpp"

////////////////////////////////////////////////////////////////////

/** The GaussianKernel class represents a convolution kernel that is described by a 2D (symmetric)
    Gaussian function. */
class GaussianKernel : public ConvolutionKernel
{
    ITEM_CONCRETE(GaussianKernel, ConvolutionKernel, "a Gaussian convolution kernel")

    PROPERTY_DOUBLE(fwhm, "the full width at half max in pixels")
        ATTRIBUTE_MIN_VALUE(fwhm, "]0")
        ATTRIBUTE_MAX_VALUE(fwhm, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(fwhm, "3")

    PROPERTY_INT(dimension, "the dimension in pixels of the kernel frame")
        ATTRIBUTE_MIN_VALUE(dimension, "1")
        ATTRIBUTE_MAX_VALUE(dimension, "1000")
        ATTRIBUTE_DEFAULT_VALUE(dimension, "6")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function creates the image frame that describes the Gaussian kernel. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
