/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshGeometry.hpp"
#include "AdaptiveMesh.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

AdaptiveMeshGeometry::~AdaptiveMeshGeometry()
{
    delete _mesh;
}

////////////////////////////////////////////////////////////////////

void AdaptiveMeshGeometry::setupSelfBefore()
{
    BoxGeometry::setupSelfBefore();

    // import the adaptive mesh
    _mesh = new AdaptiveMesh(_adaptiveMeshFile, vector<int>({_densityIndex,_multiplierIndex}), extent(), find<Log>());
    _mesh->addDensityDistribution(_densityIndex, _multiplierIndex);
    find<Log>()->info("Adaptive mesh data was successfully imported: " + std::to_string(_mesh->numCells()) + " cells.");

    // construct a vector with the normalized cumulative masses
    NR::cdf(_cumrhov, _mesh->numCells(), [this](int i){return _mesh->density(i)*_mesh->volume(i);} );
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshGeometry::density(Position bfr) const
{
    return _mesh->density(bfr) / _mesh->integratedDensity();
}

//////////////////////////////////////////////////////////////////////

Position AdaptiveMeshGeometry::generatePosition() const
{
    int m = NR::locateClip(_cumrhov, random()->uniform());
    return _mesh->randomPosition(random(), m);
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshGeometry::SigmaX() const
{
    return _mesh->SigmaX() / _mesh->integratedDensity();
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshGeometry::SigmaY() const
{
    return _mesh->SigmaY() / _mesh->integratedDensity();
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshGeometry::SigmaZ() const
{
    return _mesh->SigmaZ() / _mesh->integratedDensity();
}

//////////////////////////////////////////////////////////////////////

AdaptiveMesh* AdaptiveMeshGeometry::mesh() const
{
    return _mesh;
}

///////////////////////////////////////////////////////////////////////////////
