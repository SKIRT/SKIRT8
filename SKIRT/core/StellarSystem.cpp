////////////////////////////////////////////////////////////////////

#include "StellarSystem.hpp"
#include "NR.hpp"
#include "PhotonPackage.hpp"
#include "Random.hpp"
#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

void StellarSystem::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // cache the random generator
    _random = find<Random>();
}

//////////////////////////////////////////////////////////////////////

void StellarSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    int Nlambda = find<WavelengthGrid>()->numWavelengths();
    int Ncomp = _components.size();

    // Fill the vector _Lv with the luminosity for every wavelength bin
    _Lv.resize(Nlambda);
    for (int ell=0; ell<Nlambda; ell++)
        for (int h=0; h<Ncomp; h++)
            _Lv[ell] += _components[h]->luminosity(ell);

    // Fill the vectors _Xvv with the normalized cumulative luminosities (per wavelength bin)
    _Xvv.resize(Nlambda,0);
    for (int ell=0; ell<Nlambda; ell++)
        NR::cdf(_Xvv[ell], Ncomp, [this,ell](int h){return _components[h]->luminosity(ell);} );
}

//////////////////////////////////////////////////////////////////////

double StellarSystem::luminosity(int ell) const
{
    return _Lv[ell];
}

//////////////////////////////////////////////////////////////////////

int StellarSystem::dimension() const
{
    int result = 1;
    for (StellarComp* sc : _components) result = max(result, sc->dimension());
    return result;
}

//////////////////////////////////////////////////////////////////////

int StellarSystem::numComponents() const
{
    return _components.size();
}

//////////////////////////////////////////////////////////////////////

void StellarSystem::launch(PhotonPackage* pp, int ell, double L) const
{
    // if there is only one component, simply launch from it
    int N = numComponents();
    if (N == 1)
    {
        _components[0]->launch(pp,ell,L);
        pp->setStellarOrigin(0);
    }
    // otherwise select a component using the appriopriate biased distribution
    else
    {
        int h = 0;
        double X = _random->uniform();
        if (X<_emissionBias)
        {
            // select component from uniform distribution
            // rescale the deviate from [0,_emissionBias[ to [0,N[
            h = max(0,min(N-1,static_cast<int>(N*X/_emissionBias)));
        }
        else
        {
            // select component based on luminosity distribution
            // rescale the deviate from [_emissionBias,1[ to [0,1[
            h = NR::locateClip(_Xvv[ell],(X-_emissionBias)/(1.0-_emissionBias));
        }
        StellarComp* sc = _components[h];

        // launch a photon package from the selected component only if it has a nonzero luminosity for this wavelength
        double Lh = sc->luminosity(ell);
        if (Lh > 0)
        {
            double Lmean = _Lv[ell]/N; // the mean luminosity emitted from each stellar component
            double weight = 1.0/(1.0-_emissionBias+_emissionBias*Lmean/Lh);
            sc->launch(pp,ell,L*weight);
        }
        else
        {
            pp->launch(0., ell, Position(), Direction());
        }
        pp->setStellarOrigin(h);
    }
}

//////////////////////////////////////////////////////////////////////
