/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ReferenceImage.hpp"
#include "AdjustableSkirtSimulation.hpp"
#include "Convolution.hpp"
#include "FatalError.hpp"
#include "LumFit1.hpp"
#include "LumFit2.hpp"
#include "LumFitN.hpp"

////////////////////////////////////////////////////////////////////

void ReferenceImage::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // verify the number of luminosity boundaries vs. the number of stellar components
    size_t ncomp = find<AdjustableSkirtSimulation>()->numComponents();
    if (_minLuminosities.size() != ncomp || _maxLuminosities.size() != ncomp)
        throw FATALERROR("Number of luminosity boundaries differs from number of components " + std::to_string(ncomp));

    // verify that maximum boundary is not smaller than minimum boundary
    for (size_t k=0; k < ncomp; ++k)
    {
        if (_maxLuminosities[k] < _minLuminosities[k])
            throw FATALERROR("Maximum luminosity is smaller than minimum for component " + std::to_string(k));
    }

    // Import the reference image
    load(this, _filename);
}

////////////////////////////////////////////////////////////////////

double ReferenceImage::optimizeLuminosities(vector<Image>& inputFrames, vector<double>& luminosities) const
{
    // verify the number of input frames
    size_t ncomp = _minLuminosities.size();
    if (inputFrames.size() != ncomp)
        throw FATALERROR("Number of input frames differs from number of components " + std::to_string(ncomp));

    // convolve the input frames
    for (size_t k = 0; k < ncomp; k++)
    {
        Convolution::convolve(inputFrames[k], *_kernel);
    }

    // perform optimization algorithm depending on number of components
    luminosities.clear();
    double chi_value = 0.;
    const Image& refImage = *this;
    if (ncomp == 1)
    {
        LumFit1 lumfit;
        lumfit.setMinLum(_minLuminosities[0]);
        lumfit.setMaxLum(_maxLuminosities[0]);
        double lum;
        lumfit.optimize(refImage, inputFrames[0], lum, chi_value);
        luminosities.push_back(lum);
    }
    else if (ncomp == 2)
    {
        LumFit2 lumfit;
        lumfit.setMinLumA(_minLuminosities[0]);
        lumfit.setMaxLumA(_maxLuminosities[0]);
        lumfit.setMinLumB(_minLuminosities[1]);
        lumfit.setMaxLumB(_maxLuminosities[1]);
        double lumA, lumB;
        lumfit.optimize(refImage, inputFrames[0], inputFrames[1], lumA, lumB, chi_value);
        luminosities.push_back(lumA);
        luminosities.push_back(lumB);
    }
    else
    {
        LumFitN lumfit;
        lumfit.setMinLuminosities(_minLuminosities);
        lumfit.setMaxLuminosities(_maxLuminosities);
        lumfit.optimize(refImage, inputFrames, luminosities, chi_value);
    }
    return chi_value;
}

//////////////////////////////////////////////////////////////////////

void ReferenceImage::getTotalAndResidual(const vector<Image>& inputFrames, const vector<double>& luminosities,
                                         Image& totalImage, Image& residualImage) const
{
    size_t ncomp = inputFrames.size();
    totalImage = inputFrames[0]*luminosities[0];
    for (size_t k = 1; k < ncomp; k++) totalImage += inputFrames[k]*luminosities[k];

    const Image& refImage = *this;
    residualImage = ((refImage-totalImage)/refImage).abs();
}

////////////////////////////////////////////////////////////////////
