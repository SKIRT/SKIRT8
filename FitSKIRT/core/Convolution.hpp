/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONVOLUTION_HPP
#define CONVOLUTION_HPP

#include "Basics.hpp"
class Image;
class ConvolutionKernel;

////////////////////////////////////////////////////////////////////

/** This static class offers functions to convolve a given image with a given convolution kernel.
    */
class Convolution final
{
private:
    /** This function convolves a given image with a given convolution kernel using the Fast
        Fourier Transform (FFT) method. */
    static void fft(Image& image, const ConvolutionKernel& kernel);

    /** This function convolves a given image with a given convolution kernel using nested loops (a
        loop over the image pixels inside a loop over the kernel pixels). */
    static void loop(Image& image, const ConvolutionKernel& kernel);

public:
    /** This function convolves a given image with a given convolution kernel using one of the
        private methods defined in this class, depending on availability. */
    static void convolve(Image& Image, const ConvolutionKernel& kernel);
};

////////////////////////////////////////////////////////////////////

#endif
