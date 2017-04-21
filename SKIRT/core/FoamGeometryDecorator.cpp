/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FoamGeometryDecorator.hpp"
#include "FatalError.hpp"
#include "Foam.hpp"
#include "Log.hpp"
#include "Random.hpp"

///////////////////////////////////////////////////////////////////////////////

FoamGeometryDecorator::~FoamGeometryDecorator()
{
    delete _foam;
}

//////////////////////////////////////////////////////////////////////

void FoamGeometryDecorator::setupSelfAfter()
{
    BoxGeometry::setupSelfAfter();
    _jacobian = volume();
    _foam = Foam::createFoam(find<Log>(), random(), this, 3, _numCells);
}

///////////////////////////////////////////////////////////////////////////////

double FoamGeometryDecorator::density(Position bfr) const
{
    return _geometry->density(bfr);
}

///////////////////////////////////////////////////////////////////////////////

double FoamGeometryDecorator::SigmaX() const
{
    return _geometry->SigmaX();
}

///////////////////////////////////////////////////////////////////////////////

double FoamGeometryDecorator::SigmaY() const
{
    return _geometry->SigmaY();
}

///////////////////////////////////////////////////////////////////////////////

double FoamGeometryDecorator::SigmaZ() const
{
    return _geometry->SigmaZ();
}

///////////////////////////////////////////////////////////////////////////////

Position FoamGeometryDecorator::generatePosition() const
{
    double par[3];
    _foam->MCgenerate(par);
    return Position(fracPos(par[0], par[1], par[2]));
}

/////////////////////////////////////////////////////////////////////

double FoamGeometryDecorator::foamDensity(int ndim, double* par) const
{
    if (ndim != 3) throw FATALERROR("Incorrect dimension (ndim = " + std::to_string(ndim) + ")");
    return _geometry->density(Position(fracPos(par[0], par[1], par[2]))) * _jacobian;
}

///////////////////////////////////////////////////////////////////////////////
