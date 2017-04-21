/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ReferenceImages.hpp"
#include "AdjustableSkirtSimulation.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "Units.hpp"

//////////////////////////////////////////////////////////////////////

void ReferenceImages::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // verify the number of reference images
    auto simulation = find<AdjustableSkirtSimulation>();
    size_t numImages = _images.size();
    if (simulation->numWavelengths() != numImages)
        throw FATALERROR("Number of wavelengths in the simulation does not match number of reference images");

    // verify the size of each reference image
    for (size_t ell=0; ell<numImages; ++ell)
    {
        if (simulation->pixelsX(ell) != _images[ell]->sizeX() ||
            simulation->pixelsY(ell) != _images[ell]->sizeY())
        {
            throw FATALERROR("Instrument frame and reference image have different dimensions "
                             "for wavelength index ell = " + std::to_string(ell));
        }
    }
}

//////////////////////////////////////////////////////////////////////

double ReferenceImages::optimizeLuminosities(vector<vector<Image>>& inputFrames,
                                             vector<vector<double>>& luminosities, vector<double>& chis)
{
    if (inputFrames.size() != _images.size())
        throw FATALERROR("Number of input images does not match the number of reference images");

    luminosities.clear();
    chis.clear();
    double chi2_sum = 0;
    for (size_t ell=0; ell < _images.size() ; ++ell)
    {
        vector<double> lumin;
        double chi = _images[ell]->optimizeLuminosities(inputFrames[ell], lumin);
        luminosities.push_back(lumin);
        chis.push_back(chi);
        chi2_sum += chi;
    }
    return chi2_sum;
}

//////////////////////////////////////////////////////////////////////

void ReferenceImages::getTotalAndResidual(const vector<vector<Image>>& inputFrames,
                                          const vector<vector<double>>& luminosities,
                                          vector<Image>& totalImages, vector<Image>& residualImages) const
{
    size_t numImages = _images.size();
    totalImages.resize(numImages);
    residualImages.resize(numImages);

    for (size_t ell=0; ell < numImages ; ++ell)
    {
        _images[ell]->getTotalAndResidual(inputFrames[ell], luminosities[ell], totalImages[ell], residualImages[ell]);
    }
}

//////////////////////////////////////////////////////////////////////
