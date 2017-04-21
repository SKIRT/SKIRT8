/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BINTREEDUSTGRID_HPP
#define BINTREEDUSTGRID_HPP

#include "TreeDustGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** BinTreeDustGrid is a subclass of the TreeDustGrid class that
    implements an binary tree dust grid (2 children per node), which is in fact
    a 3-dimensional k-d tree. */
class BinTreeDustGrid : public TreeDustGrid
{
    /** The enumeration type indicating the method to be used for determining the orientation for
        each node subdivision. The Alternating method (the default) alternates repeatedly between
        x, y, and z directions in a consistent fashion. The Barycenter method chooses a subdividing
        plane parallel to the cell wall that is nearest the barycenter of the cell. */
    ENUM_DEF(DirectionMethod, Alternating, Barycenter)
    ENUM_VAL(DirectionMethod, Alternating, "alternating between x, y, and z directions")
    ENUM_VAL(DirectionMethod, Barycenter, "parallel to the cell wall nearest the barycenter")
    ENUM_END()

    ITEM_CONCRETE(BinTreeDustGrid, TreeDustGrid, "a k-d tree (binary tree) dust grid")

    PROPERTY_ENUM(directionMethod, DirectionMethod, "the method determining subdivision orientation at each level")
        ATTRIBUTE_DEFAULT_VALUE(directionMethod, "Alternating")

    ITEM_END()

protected:
    /** This function verifies that the search method has not been set to Bookkeeping,
        since that method is not compatible with a binary tree node. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

protected:
    /** This function creates a root node of type BinTreeNode using a node identifier of zero
        and the specified spatial extent, and returns a pointer to it. The caller must take
        ownership of the newly created object. */
    TreeNode* createRoot(const Box& extent) override;

    /** This virtual function returns true if this type of TreeDustGrid instance allows the use of
        the DustMassInBoxInterface for determining subdivisions, and false if not. For the
        BinTreeDustGrid class, the function returns true when using alternating subdivision, and
        false when using barycentric subdivision, because in that case the implementation of the
        BaryBinTreeNode class needs the barycenter of the node. */
    bool canUseDmibForSubdivide() override;
};

//////////////////////////////////////////////////////////////////////

#endif
