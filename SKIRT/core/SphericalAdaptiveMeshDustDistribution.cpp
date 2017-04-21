/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SphericalAdaptiveMeshDustDistribution.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SphericalAdaptiveMesh.hpp"

////////////////////////////////////////////////////////////////////

SphericalAdaptiveMeshDustDistribution::~SphericalAdaptiveMeshDustDistribution()
{
    delete _mesh;
}

//////////////////////////////////////////////////////////////////////

void SphericalAdaptiveMeshDustDistribution::setupSelfBefore()
{
    DustDistribution::setupSelfBefore();

    // verify property values
    if (_maxRadius <= _minRadius) throw FATALERROR("Domain size should be positive");
}

//////////////////////////////////////////////////////////////////////

void SphericalAdaptiveMeshDustDistribution::setupSelfAfter()
{
    DustDistribution::setupSelfAfter();

    // make a list of the field indices needed by any of our components
    vector<int> fieldIndices;
    for (auto comp : _components)
    {
        fieldIndices.push_back(comp->densityIndex());
        fieldIndices.push_back(comp->multiplierIndex());
    }

    // import the adaptive mesh
    _mesh = new SphericalAdaptiveMesh(_adaptiveMeshFile, fieldIndices, _minRadius, _maxRadius);
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

int SphericalAdaptiveMeshDustDistribution::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

int SphericalAdaptiveMeshDustDistribution::numComponents() const
{
    return _components.size();
}

//////////////////////////////////////////////////////////////////////

DustMix* SphericalAdaptiveMeshDustDistribution::mix(int h) const
{
    return _components[h]->mix();
}

//////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMeshDustDistribution::density(int h, Position bfr) const
{
    return _densityUnits * _mesh->density(h, bfr);
}

//////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMeshDustDistribution::density(Position bfr) const
{
    return _densityUnits * _mesh->density(bfr);
}

//////////////////////////////////////////////////////////////////////

Position SphericalAdaptiveMeshDustDistribution::generatePosition() const
{
    Random* random = find<Random>();
    int m = NR::locateClip(_cumrhov, random->uniform());
    return _mesh->randomPosition(random, m);
}

//////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMeshDustDistribution::mass(int h) const
{
    return _densityUnits * _mesh->integratedDensity(h);
}

//////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMeshDustDistribution::mass() const
{
    return _densityUnits * _mesh->integratedDensity();
}

//////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMeshDustDistribution::SigmaX() const
{
    return _densityUnits * _mesh->SigmaX();
}

//////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMeshDustDistribution::SigmaY() const
{
    return _densityUnits * _mesh->SigmaY();
}

//////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMeshDustDistribution::SigmaZ() const
{
    return _densityUnits * _mesh->SigmaZ();
}

//////////////////////////////////////////////////////////////////////
