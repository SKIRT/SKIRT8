/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef REFERENCEIMAGES_HPP
#define REFERENCEIMAGES_HPP

#include "ReferenceImage.hpp"

////////////////////////////////////////////////////////////////////

/** The ReferenceImages class represents a complete set of reference images. Objects of this class
    are essentially lists of pointers to ReferenceImage objects. */
class ReferenceImages: public SimulationItem
{
    ITEM_CONCRETE(ReferenceImages, SimulationItem, "a list of reference images")

    PROPERTY_ITEM_LIST(images, ReferenceImage, "the reference images")
        ATTRIBUTE_DEFAULT_VALUE(images, "ReferenceImage")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the number of reference images matches the number of instrument
        frames in the simulation, and that each reference image matches the size of the
        corresponding instrument frame. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function finds the luminosities that correspond to an optimal match of the weighted
        sum of the given input frames to the reference image for each wavelength. The input frame
        table must contain a frame per luminosity component (inner index) for each wavelength
        (outer index). The function returns the sum of the \f$\chi^2\f$ values for the optimal
        match, sets the \em chis vector to the \f$\chi^2\f$ values for each wavelength, and sets
        the \em luminosities table to optimally matching luminosities per luminosity component
        (inner index) and per wavelength (outer index). Furthermore, the \em inputFrames are
        altered in place: each frame is convolved with the reference image kernel and adjusted to
        contain the same masks as the reference image. */
    double optimizeLuminosities(vector<vector<Image>>& inputFrames,
                                vector<vector<double>>& luminosities, vector<double>& chis);

    /** Given the adjusted input frames and the optimal luminosities returned by the
        optimizeLuminosities() function, this function calculates the the optimally matching total
        images, i.e. a weighted sum of the input frames, and the residual images, i.e. the relative
        difference between the total image and the reference image. */
    void getTotalAndResidual(const vector<vector<Image>>& inputFrames, const vector<vector<double>>& luminosities,
                             vector<Image>& totalImages, vector<Image>& residualImages) const;
};

////////////////////////////////////////////////////////////////////

#endif
