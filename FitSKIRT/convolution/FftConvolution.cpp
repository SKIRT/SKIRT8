/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FftConvolution.hpp"
#include "WorkSpace.hpp"

////////////////////////////////////////////////////////////////////

FftConvolution::FftConvolution(int input_xsize, int input_ysize, int kernel_xsize, int kernel_ysize)
{
#ifdef BUILD_WITH_FFT
    // Create a new workspace
    _ws = new WorkSpace();

    // Initialize the workspace
    _ws->initialize(ConvolutionMode::LINEAR_SAME, input_xsize, input_ysize, kernel_xsize, kernel_ysize);
#else
    _ws = nullptr;
    (void)input_xsize; (void)input_ysize; (void)kernel_xsize; (void)kernel_ysize;
#endif
}

////////////////////////////////////////////////////////////////////

FftConvolution::~FftConvolution()
{
#ifdef BUILD_WITH_FFT
    // Clear the workspace
    _ws->clear();

    // Delete the workspace
    delete _ws;
#endif
}

////////////////////////////////////////////////////////////////////

void FftConvolution::perform(const Array& input, const Array& kernel, Array& output)
{
#ifdef BUILD_WITH_FFT
    // Do the convolution
    _ws->convolve(input, kernel, output);
#else
    (void)input; (void)kernel; (void)output;
#endif
}

////////////////////////////////////////////////////////////////////

bool FftConvolution::enabled()
{
#ifdef BUILD_WITH_FFT
    return true;
#else
    return false;
#endif
}

////////////////////////////////////////////////////////////////////
