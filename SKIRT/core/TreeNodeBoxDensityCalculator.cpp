/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TreeNodeBoxDensityCalculator.hpp"
#include "Constants.hpp"
#include "DustMassInBoxInterface.hpp"
#include "FatalError.hpp"
#include "TreeNode.hpp"

//////////////////////////////////////////////////////////////////////

TreeNodeBoxDensityCalculator::TreeNodeBoxDensityCalculator(DustMassInBoxInterface* dmib, TreeNode* node)
    : _dmib(dmib), _extent(node->extent()), _mass(-1.)
{
}

//////////////////////////////////////////////////////////////////////

double TreeNodeBoxDensityCalculator::volume() const
{
    return _extent.volume();
}

//////////////////////////////////////////////////////////////////////

double TreeNodeBoxDensityCalculator::mass() const
{
    if (_mass < 0)
    {
        const_cast<TreeNodeBoxDensityCalculator*>(this)->_mass = _dmib->massInBox(_extent);
    }
    return _mass;
}

//////////////////////////////////////////////////////////////////////

Vec TreeNodeBoxDensityCalculator::barycenter() const
{
    throw FATALERROR("Calculation is not supported");
}

//////////////////////////////////////////////////////////////////////

double TreeNodeBoxDensityCalculator::densityDispersion() const
{
    throw FATALERROR("Calculation is not supported");
}

//////////////////////////////////////////////////////////////////////

double TreeNodeBoxDensityCalculator::opticalDepth() const
{
    return Constants::kappaV() * mass() / pow(volume(),2./3.);
}

//////////////////////////////////////////////////////////////////////
