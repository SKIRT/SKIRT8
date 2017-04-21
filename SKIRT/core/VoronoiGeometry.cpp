/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoronoiGeometry.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "VoronoiMesh.hpp"
#include "VoronoiMeshFile.hpp"

////////////////////////////////////////////////////////////////////

VoronoiGeometry::~VoronoiGeometry()
{
    delete _mesh;
}

////////////////////////////////////////////////////////////////////

void VoronoiGeometry::setupSelfBefore()
{
    BoxGeometry::setupSelfBefore();

    // import the Voronoi mesh
    _mesh = new VoronoiMesh(_voronoiMeshFile, vector<int>({_densityIndex, _multiplierIndex}), extent());
    _mesh->addDensityDistribution(_densityIndex, _multiplierIndex);
    find<Log>()->info("Voronoi mesh data was successfully imported: " + std::to_string(_mesh->numCells()) + " cells.");

    // construct a vector with the normalized cumulative masses
    NR::cdf(_cumrhov, _mesh->numCells(), [this](int i){return _mesh->density(i)*_mesh->volume(i);} );
}

//////////////////////////////////////////////////////////////////////

double VoronoiGeometry::density(Position bfr) const
{
    return _mesh->density(bfr) / _mesh->integratedDensity();
}

//////////////////////////////////////////////////////////////////////

Position VoronoiGeometry::generatePosition() const
{
    int m = NR::locateClip(_cumrhov, random()->uniform());
    return _mesh->randomPosition(random(), m);
}

//////////////////////////////////////////////////////////////////////

double VoronoiGeometry::SigmaX() const
{
    return _mesh->SigmaX() / _mesh->integratedDensity();
}

//////////////////////////////////////////////////////////////////////

double VoronoiGeometry::SigmaY() const
{
    return _mesh->SigmaY() / _mesh->integratedDensity();
}

//////////////////////////////////////////////////////////////////////

double VoronoiGeometry::SigmaZ() const
{
    return _mesh->SigmaZ() / _mesh->integratedDensity();
}

///////////////////////////////////////////////////////////////////////////////

VoronoiMesh* VoronoiGeometry::mesh() const
{
    return _mesh;
}

//////////////////////////////////////////////////////////////////////

int VoronoiGeometry::numParticles() const
{
    return _mesh->numCells();
}

//////////////////////////////////////////////////////////////////////

Vec VoronoiGeometry::particleCenter(int index) const
{
    return _mesh->particlePosition(index);
}

//////////////////////////////////////////////////////////////////////
