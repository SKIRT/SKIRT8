/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CompDustDistribution.hpp"
#include "Geometry.hpp"
#include "NR.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void CompDustDistribution::setupSelfAfter()
{
    DustDistribution::setupSelfAfter();

    // construct a vector with the normalized cumulative masses
    NR::cdf(_cumrhov, _components.size(), [this](int i){return _components[i]->mass();} );
}

////////////////////////////////////////////////////////////////////

int CompDustDistribution::dimension() const
{
    int result = 1;
    for (auto comp : _components) result = max(result, comp->dimension());
    return result;
}

//////////////////////////////////////////////////////////////////////

int CompDustDistribution::numComponents() const
{
    return _components.size();
}

//////////////////////////////////////////////////////////////////////

DustMix* CompDustDistribution::mix(int h) const
{
    return _components[h]->mix();
}

//////////////////////////////////////////////////////////////////////

double CompDustDistribution::density(int h, Position bfr) const
{
    return _components[h]->density(bfr);
}

//////////////////////////////////////////////////////////////////////

double CompDustDistribution::density(Position bfr) const
{
    double sum = 0.0;
    for (auto comp : _components) sum += comp->density(bfr);
    return sum;
}

//////////////////////////////////////////////////////////////////////

Position CompDustDistribution::generatePosition() const
{
    Random* random = find<Random>();
    int h = NR::locateClip(_cumrhov, random->uniform());
    return _components[h]->geometry()->generatePosition();
}

//////////////////////////////////////////////////////////////////////

double CompDustDistribution::mass(int h) const
{
    return _components[h]->mass();
}

//////////////////////////////////////////////////////////////////////

double CompDustDistribution::mass() const
{
    double sum = 0.0;
    for (auto comp : _components) sum += comp->mass();
    return sum;
}

//////////////////////////////////////////////////////////////////////

double CompDustDistribution::SigmaX() const
{
    double sum = 0.0;
    for (auto comp : _components) sum += comp->SigmaX();
    return sum;
}

//////////////////////////////////////////////////////////////////////

double CompDustDistribution::SigmaY() const
{
    double sum = 0.0;
    for (auto comp : _components) sum += comp->SigmaY();
    return sum;
}

//////////////////////////////////////////////////////////////////////

double CompDustDistribution::SigmaZ() const
{
    double sum = 0.0;
    for (auto comp : _components) sum += comp->SigmaZ();
    return sum;
}

//////////////////////////////////////////////////////////////////////

vector<SimulationItem*> CompDustDistribution::interfaceCandidates(const std::type_info& interfaceTypeInfo)
{
    vector<SimulationItem*> candidates = DustDistribution::interfaceCandidates(interfaceTypeInfo);
    if (_components.size() == 1) candidates.push_back(_components[0]->geometry());
    return candidates;
}

//////////////////////////////////////////////////////////////////////
