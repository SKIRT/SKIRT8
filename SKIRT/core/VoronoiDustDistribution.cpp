/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoronoiDustDistribution.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "VoronoiMesh.hpp"

////////////////////////////////////////////////////////////////////

VoronoiDustDistribution::~VoronoiDustDistribution()
{
    delete _mesh;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustDistribution::setupSelfAfter()
{
    BoxDustDistribution::setupSelfAfter();

    // make a list of the field indices needed by any of our components
    vector<int> fieldIndices;
    for (auto comp : _components)
    {
        fieldIndices.push_back(comp->densityIndex());
        fieldIndices.push_back(comp->multiplierIndex());
    }

    // import the Voronoi mesh
    _mesh = new VoronoiMesh(_voronoiMeshFile, fieldIndices, extent());
    find<Log>()->info("Voronoi mesh data was successfully imported: " + std::to_string(_mesh->numCells()) + " cells.");

    // add a density field for each of our components, so that the mesh holds the total density
    for (auto comp : _components)
    {
        _mesh->addDensityDistribution(comp->densityIndex(), comp->multiplierIndex(), comp->densityFraction());
    }

    // construct a vector with the normalized cumulative masses
    NR::cdf(_cumrhov, _mesh->numCells(), [this](int i){return _mesh->density(i)*_mesh->volume(i);} );
}

//////////////////////////////////////////////////////////////////////

int VoronoiDustDistribution::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

int VoronoiDustDistribution::numComponents() const
{
    return _components.size();
}

//////////////////////////////////////////////////////////////////////

DustMix* VoronoiDustDistribution::mix(int h) const
{
    return _components[h]->mix();
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustDistribution::density(int h, Position bfr) const
{
    return _densityUnits * _mesh->density(h, bfr);
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustDistribution::density(Position bfr) const
{
    return _densityUnits * _mesh->density(bfr);
}

//////////////////////////////////////////////////////////////////////

Position VoronoiDustDistribution::generatePosition() const
{
    Random* random = find<Random>();
    int m = NR::locateClip(_cumrhov, random->uniform());
    return _mesh->randomPosition(random, m);
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustDistribution::mass(int h) const
{
    return _densityUnits * _mesh->integratedDensity(h);
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustDistribution::mass() const
{
    return _densityUnits * _mesh->integratedDensity();
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustDistribution::SigmaX() const
{
    return _densityUnits * _mesh->SigmaX();
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustDistribution::SigmaY() const
{
    return _densityUnits * _mesh->SigmaY();
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustDistribution::SigmaZ() const
{
    return _densityUnits * _mesh->SigmaZ();
}

//////////////////////////////////////////////////////////////////////

VoronoiMesh* VoronoiDustDistribution::mesh() const
{
    return _mesh;
}

//////////////////////////////////////////////////////////////////////

int VoronoiDustDistribution::numParticles() const
{
    return _mesh->numCells();
}

//////////////////////////////////////////////////////////////////////

Vec VoronoiDustDistribution::particleCenter(int index) const
{
    return _mesh->particlePosition(index);
}

//////////////////////////////////////////////////////////////////////
