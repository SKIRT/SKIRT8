/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TwoPhaseDustGrid.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void TwoPhaseDustGrid::setupSelfAfter()
{
    CartesianDustGrid::setupSelfAfter();
    Random* random = find<Random>();

    // construction of the weight matrix
    _weightv.resize(numCells());
    for (int m=0; m<numCells(); m++)
    {
        double X = random->uniform();
        _weightv[m] = (X<_fillingFactor) ?
                          _contrast/(_contrast*_fillingFactor+1.0-_fillingFactor) :
                          1.0/(_contrast*_fillingFactor+1.0-_fillingFactor);
    }
}

//////////////////////////////////////////////////////////////////////

double TwoPhaseDustGrid::weight(int m) const
{
    return (m==-1) ? 0.0 : _weightv[m];
}

//////////////////////////////////////////////////////////////////////
