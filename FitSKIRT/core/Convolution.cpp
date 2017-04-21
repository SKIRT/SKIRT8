/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Convolution.hpp"
#include "ConvolutionKernel.hpp"
#include "FatalError.hpp"
#include "FftConvolution.hpp"
#include "Image.hpp"

////////////////////////////////////////////////////////////////////

void Convolution::fft(Image& image, const ConvolutionKernel& kernel)
{
    // Test whether FFT convolution is possible
    if (!FftConvolution::enabled()) throw FATALERROR("The library required for FFT is not present");

    // Initialize an output array
    Array output(image.size());

    // Create an FftConvolution object
    FftConvolution fftc(image.sizeX(), image.sizeY(), kernel.sizeX(), kernel.sizeY());

    // Perform the convolution
    fftc.perform(image.data(), kernel.data(), output);

    // Move the output array to the image
    image.moveData(std::move(output));
}

////////////////////////////////////////////////////////////////////

void Convolution::loop(Image& image, const ConvolutionKernel& kernel)
{
    // Initialize a convolved image
    Array output(image.size());

    // Do the convolution by looping over all input frame pixels and all kernel pixels
    for (int l = 0; l < image.size(); l++)
    {
        // Determine the x and y coordinate of the input image
        int yi = l/image.sizeX();
        int xi = l-image.sizeX()*yi;

        // Loop over the kernel pixels in the x direction
        for (int xk = 0; xk < kernel.sizeX(); xk++)
        {
            // Loop over the kernel pixels in the y direction
            for (int yk = 0; yk < kernel.sizeY(); yk++)
            {
                int startX = xi - (kernel.sizeX()-1)/2;
                int startY = yi - (kernel.sizeY()-1)/2;

                // Exclude pixels out of bound
                if (startX + xk >= 0 && startY + yk >= 0 && startX + xk < image.sizeX() && startY + yk < image.sizeY())
                {
                    output[(startX + xk)+image.sizeX()*(startY + yk)] += image[xi+image.sizeX()*(yi)] * kernel(xk,yk);
                }
            }
        }
    }

    // Move the output array to the image
    image.moveData(std::move(output));
}

////////////////////////////////////////////////////////////////////

void Convolution::convolve(Image& image, const ConvolutionKernel& kernel)
{
    // Use the Fast Fourier Transform method for the convolution if the kernel is sufficiently large
    if (kernel.size() > 200 && FftConvolution::enabled()) fft(image, kernel);

    // Else, use the standard nested loop that iterates over all image and kernel pixels
    else loop(image, kernel);
}

////////////////////////////////////////////////////////////////////
