/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONVOLUTIONKERNEL_HPP
#define CONVOLUTIONKERNEL_HPP

#include "SimulationItem.hpp"
#include "Image.hpp"

////////////////////////////////////////////////////////////////////

/** The ConvolutionKernel class is used to describe a general convolution kernel. Subclasses of
    this class represent specific types of convolution kernels. */
class ConvolutionKernel : public SimulationItem, public Image
{
    ITEM_ABSTRACT(ConvolutionKernel, SimulationItem, "a convolution kernel")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** The purpose of this function, which is re-implemented in each derived class, is to initialize
        the data that describes the kernel. */
    void setupSelfBefore() override;

    /** This function, only implemented in the base class, makes sure that the kernel is properly
        normalized. */
    void setupSelfAfter() override;
};

////////////////////////////////////////////////////////////////////

#endif
