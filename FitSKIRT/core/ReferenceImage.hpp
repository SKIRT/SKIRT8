/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef REFERENCEIMAGE_HPP
#define REFERENCEIMAGE_HPP

#include "SimulationItem.hpp"
#include "Image.hpp"
#include "ConvolutionKernel.hpp"

////////////////////////////////////////////////////////////////////

/** The ReferenceImage class contains all information for a certain reference image.
    It contains the path to the file and information about the PSF.
    The \f$\chi^2\f$ value for a certain input frame can be calculated. */
class ReferenceImage : public SimulationItem, public Image
{
    ITEM_CONCRETE(ReferenceImage, SimulationItem, "a reference image")

    PROPERTY_STRING(filename, "the name (or relative path) of this reference image")

    PROPERTY_ITEM(kernel, ConvolutionKernel, "the convolution kernel")
        ATTRIBUTE_DEFAULT_VALUE(kernel, "GaussianKernel")

    PROPERTY_DOUBLE_LIST(minLuminosities, "the minimum luminosity factors, one per stellar component")
        ATTRIBUTE_MIN_VALUE(minLuminosities, "]0")

    PROPERTY_DOUBLE_LIST(maxLuminosities, "the maximum luminosity factors, one per stellar component")
        ATTRIBUTE_MIN_VALUE(maxLuminosities, "]0")

    ITEM_END()

    //============ Construction - Setup - Destruction  =============

protected:
    /** This function reads the reference image file with the given name into memory. */
    void setupSelfBefore() override;

    //====================== Other functions =======================

public:
    /** This function finds the luminosities that correspond to an optimal match of the weighted
        sum of the given input frames to the reference image. There must be one input frame per
        luminosity component in the simulation. The function returns the \f$\chi^2\f$ value for the
        optimal match, and sets the \em luminosities vector to the list of optimally matching
        luminosities (one per input frame). Furthermore, the \em inputFrames are altered in place:
        each frame is convolved with the reference image kernel and adjusted to contain the same
        masks as the reference image. */
    double optimizeLuminosities(vector<Image>& inputFrames, vector<double>& luminosities) const;

    /** Given the adjusted input frames and the optimal luminosities returned by the
        optimizeLuminosities() function, this function calculates the the optimally matching total
        image, i.e. a weighted sum of the input frames, and the residual image, i.e. the relative
        difference between the total image and the reference image. */
    void getTotalAndResidual(const vector<Image>& inputFrames, const vector<double>& luminosities,
                             Image& totalImage, Image& residualImage) const;
};

////////////////////////////////////////////////////////////////////

#endif
