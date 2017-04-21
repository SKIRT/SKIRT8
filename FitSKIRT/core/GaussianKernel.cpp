/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GaussianKernel.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void GaussianKernel::setupSelfBefore()
{
    ConvolutionKernel::setupSelfBefore();

    // From the FWHM, calculate the standard deviation
    double sigma = _fwhm/2.3548;

    // Resize the kernel image
    resize(_dimension, _dimension);

    // Set the kernel data
    for (int yk=1; yk<_dimension+1; yk++)
    {
        for (int xk=1; xk< _dimension+1; xk++)
        {
            int Xi = (_dimension+1)/2 - xk;
            int Yi = (_dimension+1)/2 - yk;

            (*this)(xk-1,yk-1) = exp(-(Xi*Xi+Yi*Yi)/(2*sigma*sigma));
        }
    }
}

////////////////////////////////////////////////////////////////////
