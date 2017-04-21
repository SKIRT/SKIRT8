/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BinTreeDustGrid.hpp"
#include "BaryBinTreeNode.hpp"
#include "FatalError.hpp"

//////////////////////////////////////////////////////////////////////

void BinTreeDustGrid::setupSelfBefore()
{
    TreeDustGrid::setupSelfBefore();

    if (searchMethod() == TreeDustGrid::SearchMethod::Bookkeeping)
        throw FATALERROR("Bookkeeping method is not compatible with binary tree");
}

//////////////////////////////////////////////////////////////////////

TreeNode* BinTreeDustGrid::createRoot(const Box& extent)
{
    switch (_directionMethod)
    {
    case DirectionMethod::Barycenter:
        return new BaryBinTreeNode(0, 0, extent);
    case DirectionMethod::Alternating:
        return new BinTreeNode(0, 0, extent);
    }
    return nullptr;
}

//////////////////////////////////////////////////////////////////////

bool BinTreeDustGrid::canUseDmibForSubdivide()
{
    return _directionMethod != DirectionMethod::Barycenter;
}

//////////////////////////////////////////////////////////////////////
