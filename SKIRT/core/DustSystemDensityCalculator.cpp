/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustSystemDensityCalculator.hpp"
#include "DustDistribution.hpp"
#include "DustGrid.hpp"
#include "DustSystem.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

DustSystemDensityCalculator::DustSystemDensityCalculator(const DustSystem *ds, int numBodies, int numSamplesPerBody)
    : _ds(ds), _dd(ds->dustDistribution()), _grid(ds->dustGrid()), _random(ds->find<Random>()),
      _numBodies(numBodies), _numSamplesPerBody(numSamplesPerBody),
      _extent(_grid->boundingBox()),
      _drhov(numBodies), _drho2v(numBodies),
      _meanDelta(0), _stddevDelta(0),
      _consolidated(false)
{
}

//////////////////////////////////////////////////////////////////////

void DustSystemDensityCalculator::body(size_t n)
{
    int k = _numSamplesPerBody;
    while (k--)
    {
        Position pos = _random->position(_extent);
        double rhot = _dd->density(pos);
        double rhog = _ds->density(_grid->whichCell(pos));
        double drho = fabs(rhog-rhot);
        _drhov[n]  += drho;
        _drho2v[n] += drho*drho;
    }
    _drhov[n]  /= _numSamplesPerBody;
    _drho2v[n] /= _numSamplesPerBody;
}

//////////////////////////////////////////////////////////////////////

void DustSystemDensityCalculator::consolidate()
{
    _meanDelta = _drhov.sum()/_numBodies;
    _stddevDelta = sqrt(_drho2v.sum()/_numBodies - _meanDelta*_meanDelta);
    _consolidated = true;
}

//////////////////////////////////////////////////////////////////////

double DustSystemDensityCalculator::meanDelta()
{
    if (!_consolidated) consolidate();
    return _meanDelta;
}

//////////////////////////////////////////////////////////////////////

double DustSystemDensityCalculator::stdDevDelta()
{
    if (!_consolidated) consolidate();
    return _stddevDelta;
}

//////////////////////////////////////////////////////////////////////
