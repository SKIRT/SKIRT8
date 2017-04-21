/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "WavelengthGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "StaggeredAssigner.hpp"

////////////////////////////////////////////////////////////////////

void WavelengthGrid::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // the subclass should have added wavelengths
    int Nlambda = _lambdav.size();
    if (!Nlambda) throw FATALERROR("Wavelength grid should have been initialized by subclass");

    // in data parallelization mode, create an assigner
    PeerToPeerCommunicator* comm = find<PeerToPeerCommunicator>();
    if (comm->dataParallel())
    {
        if (comm->size() > Nlambda) throw FATALERROR("In data parallelization mode, the number of wavelengths "
                                                     "must be at least the number of processes.");
        _assigner = new StaggeredAssigner(Nlambda, this);
    }
}

////////////////////////////////////////////////////////////////////

void WavelengthGrid::setWavelengthBins(const Array& lambdav, const Array& dlambdav)
{
    // verify length of incoming arrays
    int Nlambda = lambdav.size();
    if (!Nlambda) throw FATALERROR("There must be at least one wavelength in the grid");
    if (lambdav.size() != dlambdav.size()) throw FATALERROR("Number of wavelengths and bin widths must match");

    // wavelengths should be positive and sorted in ascending order
    if (lambdav[0] <= 0.0) throw FATALERROR("All wavelengths should be positive");
    for (int ell = 1; ell < Nlambda; ell++)
    {
        if (lambdav[ell] <= lambdav[ell-1]) throw FATALERROR("Wavelengths should be sorted in ascending order");
    }

    // remember the wavelength grid
    _lambdav = lambdav;
    _dlambdav = dlambdav;
}

////////////////////////////////////////////////////////////////////

const ProcessAssigner* WavelengthGrid::assigner() const
{
    return _assigner;
}

////////////////////////////////////////////////////////////////////

int WavelengthGrid::numWavelengths() const
{
    return _lambdav.size();
}

////////////////////////////////////////////////////////////////////

double WavelengthGrid::lambda(int ell) const
{
    return _lambdav[ell];
}

////////////////////////////////////////////////////////////////////

double WavelengthGrid::dlambda(int ell) const
{
    return _dlambdav[ell];
}

////////////////////////////////////////////////////////////////////

double WavelengthGrid::lambdamin(int ell) const
{
    return (ell==0) ? _lambdav[0] : sqrt(_lambdav[ell-1]*_lambdav[ell]);
}

////////////////////////////////////////////////////////////////////

double WavelengthGrid::lambdamax(int ell) const
{
    int Nlambda = _lambdav.size();
    return (ell==Nlambda-1) ? _lambdav[Nlambda-1] : sqrt(_lambdav[ell]*_lambdav[ell+1]);
}

////////////////////////////////////////////////////////////////////

int WavelengthGrid::nearest(double lambda) const
{
    int ell = NR::locateFail(_lambdav,lambda);
    if (ell<0) return -1;
    double lambdac = sqrt(_lambdav[ell]*_lambdav[ell+1]);
    if (lambda<lambdac) return ell;
    else return ell+1;
}

////////////////////////////////////////////////////////////////////

const Array& WavelengthGrid::lambdav() const
{
    return _lambdav;
}

////////////////////////////////////////////////////////////////////

const Array& WavelengthGrid::dlambdav() const
{
    return _dlambdav;
}

////////////////////////////////////////////////////////////////////
