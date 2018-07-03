/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BaryBinTreeNode.hpp"
#include "TreeNodeDensityCalculator.hpp"

//////////////////////////////////////////////////////////////////////

// constants for splitting direction (not an enum because we use modulo arithmetic)
#define XDIR 0
#define YDIR 1
#define ZDIR 2

//////////////////////////////////////////////////////////////////////

BaryBinTreeNode::BaryBinTreeNode(TreeNode* father, int id, const Box& extent)
    : BinTreeNode(father, id, extent)
{
}

//////////////////////////////////////////////////////////////////////

TreeNode* BaryBinTreeNode::createNode(TreeNode* father, int id, const Box& extent)
{
    return new BaryBinTreeNode(father, id, extent);
}

//////////////////////////////////////////////////////////////////////

void BaryBinTreeNode::createChildren(int id, const TreeNodeDensityCalculator* calc)
{
    Vec b = calc->barycenter();

    // fractional distance to nearest wall, in each direction
    double dx = min(b.x()-xmin(), xmax()-b.x())/(xmax()-xmin());
    double dy = min(b.y()-ymin(), ymax()-b.y())/(ymax()-ymin());
    double dz = min(b.z()-zmin(), zmax()-b.z())/(zmax()-zmin());

    // select the direction with the smallest relative distance
    int dir;
    if (dx < dy)
    {
        if (dx < dz) dir = XDIR;
        else dir = ZDIR;
    }
    else
    {
        if (dy < dz) dir = YDIR;
        else dir = ZDIR;
    }

    // split along that direction
    createChildrenSplitDir(id, dir);
}

//////////////////////////////////////////////////////////////////////
