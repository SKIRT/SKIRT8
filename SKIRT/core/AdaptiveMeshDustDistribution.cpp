/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshDustDistribution.hpp"
#include "AdaptiveMesh.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "MeshDustComponent.hpp"
#include "NR.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

AdaptiveMeshDustDistribution::~AdaptiveMeshDustDistribution()
{
    delete _mesh;
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshDustDistribution::setupSelfAfter()
{
    BoxDustDistribution::setupSelfAfter();

    // make a list of the field indices needed by any of our components
    vector<int> fieldIndices;
    for (auto comp : _components)
    {
        fieldIndices.push_back(comp->densityIndex());
        fieldIndices.push_back(comp->multiplierIndex());
    }

    // import the adaptive mesh
    _mesh = new AdaptiveMesh(_adaptiveMeshFile, fieldIndices, extent(), find<Log>());
    find<Log>()->info("Adaptive mesh data was successfully imported: " + std::to_string(_mesh->numCells()) + " cells.");

    // add a density field for each of our components, so that the mesh holds the total density
    for (auto comp : _components)
    {
        _mesh->addDensityDistribution(comp->densityIndex(), comp->multiplierIndex(), comp->densityFraction());
    }

    // construct a vector with the normalized cumulative masses
    NR::cdf(_cumrhov, _mesh->numCells(), [this](int i){return _mesh->density(i)*_mesh->volume(i);} );
}

//////////////////////////////////////////////////////////////////////

int AdaptiveMeshDustDistribution::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

int AdaptiveMeshDustDistribution::numComponents() const
{
    return _components.size();
}

//////////////////////////////////////////////////////////////////////

DustMix* AdaptiveMeshDustDistribution::mix(int h) const
{
    return _components[h]->mix();
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshDustDistribution::density(int h, Position bfr) const
{
    return _densityUnits * _mesh->density(h, bfr);
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshDustDistribution::density(Position bfr) const
{
    return _densityUnits * _mesh->density(bfr);
}

//////////////////////////////////////////////////////////////////////

Position AdaptiveMeshDustDistribution::generatePosition() const
{
    Random* random = find<Random>();
    int m = NR::locateClip(_cumrhov, random->uniform());
    return _mesh->randomPosition(random, m);
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshDustDistribution::mass() const
{
    return _densityUnits * _mesh->integratedDensity();
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshDustDistribution::mass(int h) const
{
    return _densityUnits * _mesh->integratedDensity(h);
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshDustDistribution::SigmaX() const
{
    return _densityUnits * _mesh->SigmaX();
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshDustDistribution::SigmaY() const
{
    return _densityUnits * _mesh->SigmaY();
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshDustDistribution::SigmaZ() const
{
    return _densityUnits * _mesh->SigmaZ();
}

//////////////////////////////////////////////////////////////////////

AdaptiveMesh* AdaptiveMeshDustDistribution::mesh() const
{
    return _mesh;
}

//////////////////////////////////////////////////////////////////////
