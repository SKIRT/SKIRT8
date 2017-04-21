/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OCTTREEDUSTGRID_HPP
#define OCTTREEDUSTGRID_HPP

#include "TreeDustGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** OctTreeDustGrid is a subclass of the TreeDustGrid class that implements an octtree dust grid (8
    children per node). The class supports geometric and barycentric subdivision of cells. If the
    \em barycentric flag is set to false (the default), cells are subdivided in their geometric
    center. If the flag is set to true, cells are subdivided in their center of mass (barycenter).
    */
class OctTreeDustGrid : public TreeDustGrid
{
    ITEM_CONCRETE(OctTreeDustGrid, TreeDustGrid, "an octtree dust grid")

    PROPERTY_BOOL(useBarycentric, "use barycentric subdivision")
        ATTRIBUTE_DEFAULT_VALUE(useBarycentric, "false")

    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function creates a root node of type OctTreeNode using a node identifier of zero
        and the specified spatial extent, and returns a pointer to it. The caller must take
        ownership of the newly created object. */
    TreeNode* createRoot(const Box& extent) override;

    /** This virtual function returns true if this type of TreeDustGrid instance allows the use of
        the DustMassInBoxInterface for determining subdivisions, and false if not. For the
        OctTreeDustGrid class, the function returns true when using geometric subdivision, and
        false when using barycentric subdivision, because in that case the implementation of the
        BaryOctTreeNode class needs the barycenter of the node. */
    bool canUseDmibForSubdivide() override;
};

//////////////////////////////////////////////////////////////////////

#endif
