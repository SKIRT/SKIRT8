/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OctTreeNode.hpp"

//////////////////////////////////////////////////////////////////////

// macros for easily accessing a particular child
#define CHILD_0 _children[0]
#define CHILD_1 _children[1]
#define CHILD_2 _children[2]
#define CHILD_3 _children[3]
#define CHILD_4 _children[4]
#define CHILD_5 _children[5]
#define CHILD_6 _children[6]
#define CHILD_7 _children[7]

//////////////////////////////////////////////////////////////////////

OctTreeNode::OctTreeNode(TreeNode* father, int id, const Box& extent)
    : TreeNode(father, id, extent)
{
}

//////////////////////////////////////////////////////////////////////

TreeNode* OctTreeNode::createNode(TreeNode* father, int id, const Box& extent)
{
    return new OctTreeNode(father, id, extent);
}

//////////////////////////////////////////////////////////////////////

void OctTreeNode::createChildrenSplitPoint(int id, Vec r)
{
    _children.resize(8);
    CHILD_0 = createNode(this, id++, Box(_xmin, _ymin, _zmin, r.x(), r.y(), r.z()));
    CHILD_1 = createNode(this, id++, Box(r.x(), _ymin, _zmin, _xmax, r.y(), r.z()));
    CHILD_2 = createNode(this, id++, Box(_xmin, r.y(), _zmin, r.x(), _ymax, r.z()));
    CHILD_3 = createNode(this, id++, Box(r.x(), r.y(), _zmin, _xmax, _ymax, r.z()));
    CHILD_4 = createNode(this, id++, Box(_xmin, _ymin, r.z(), r.x(), r.y(), _zmax));
    CHILD_5 = createNode(this, id++, Box(r.x(), _ymin, r.z(), _xmax, r.y(), _zmax));
    CHILD_6 = createNode(this, id++, Box(_xmin, r.y(), r.z(), r.x(), _ymax, _zmax));
    CHILD_7 = createNode(this, id++, Box(r.x(), r.y(), r.z(), _xmax, _ymax, _zmax));
}

//////////////////////////////////////////////////////////////////////

void OctTreeNode::createChildren(int id)
{
    createChildrenSplitPoint(id, center());
}

//////////////////////////////////////////////////////////////////////

void OctTreeNode::createChildren(int id, const TreeNodeDensityCalculator* /*calc*/)
{
    createChildren(id);
}

//////////////////////////////////////////////////////////////////////

void OctTreeNode::addNeighbors()
{
    // if we don't have any children, we can't add neighbors
    if (_children.empty()) return;

    // ensure that all involved nodes have a neighbor list for each of the walls
    ensureNeighborLists();
    CHILD_0->ensureNeighborLists();
    CHILD_1->ensureNeighborLists();
    CHILD_2->ensureNeighborLists();
    CHILD_3->ensureNeighborLists();
    CHILD_4->ensureNeighborLists();
    CHILD_5->ensureNeighborLists();
    CHILD_6->ensureNeighborLists();
    CHILD_7->ensureNeighborLists();

    // Internal neighbors: each of the 8 new children has 3 neighbors among its siblings
    makeNeighbors(FRONT, CHILD_0, CHILD_1);
    makeNeighbors(RIGHT, CHILD_0, CHILD_2);
    makeNeighbors(TOP  , CHILD_0, CHILD_4);
    makeNeighbors(RIGHT, CHILD_1, CHILD_3);
    makeNeighbors(TOP  , CHILD_1, CHILD_5);
    makeNeighbors(FRONT, CHILD_2, CHILD_3);
    makeNeighbors(TOP  , CHILD_2, CHILD_6);
    makeNeighbors(TOP  , CHILD_3, CHILD_7);
    makeNeighbors(FRONT, CHILD_4, CHILD_5);
    makeNeighbors(RIGHT, CHILD_4, CHILD_6);
    makeNeighbors(RIGHT, CHILD_5, CHILD_7);
    makeNeighbors(FRONT, CHILD_6, CHILD_7);

    // The point where this node is split into its children
    double xc = CHILD_0->xmax();
    double yc = CHILD_0->ymax();
    double zc = CHILD_0->zmax();

    // The BACK neighbors of this node
    {
        const vector<TreeNode*>& neighbors = _neighbors[BACK];
        for (unsigned int i=0; i<neighbors.size(); i++)
        {
            TreeNode* neighbor = neighbors[i];
            neighbor->deleteNeighbor(FRONT, this);
            if (neighbor->ymin() <= yc  &&  neighbor->zmin() <= zc) makeNeighbors(FRONT, neighbor, CHILD_0);
            if (neighbor->ymax() >= yc  &&  neighbor->zmin() <= zc) makeNeighbors(FRONT, neighbor, CHILD_2);
            if (neighbor->ymin() <= yc  &&  neighbor->zmax() >= zc) makeNeighbors(FRONT, neighbor, CHILD_4);
            if (neighbor->ymax() >= yc  &&  neighbor->zmax() >= zc) makeNeighbors(FRONT, neighbor, CHILD_6);
        }
    }
    // The FRONT neighbors of this node
    {
        const vector<TreeNode*>& neighbors = _neighbors[FRONT];
        for (unsigned int i=0; i<neighbors.size(); i++)
        {
            TreeNode* neighbor = neighbors[i];
            neighbor->deleteNeighbor(BACK, this);
            if (neighbor->ymin() <= yc  &&  neighbor->zmin() <= zc) makeNeighbors(BACK, neighbor, CHILD_1);
            if (neighbor->ymax() >= yc  &&  neighbor->zmin() <= zc) makeNeighbors(BACK, neighbor, CHILD_3);
            if (neighbor->ymin() <= yc  &&  neighbor->zmax() >= zc) makeNeighbors(BACK, neighbor, CHILD_5);
            if (neighbor->ymax() >= yc  &&  neighbor->zmax() >= zc) makeNeighbors(BACK, neighbor, CHILD_7);
        }
    }
    // The LEFT neighbors of this node
    {
        const vector<TreeNode*>& neighbors = _neighbors[LEFT];
        for (unsigned int i=0; i<neighbors.size(); i++)
        {
            TreeNode* neighbor = neighbors[i];
            neighbor->deleteNeighbor(RIGHT, this);
            if (neighbor->xmin() <= xc  &&  neighbor->zmin() <= zc) makeNeighbors(RIGHT, neighbor, CHILD_0);
            if (neighbor->xmax() >= xc  &&  neighbor->zmin() <= zc) makeNeighbors(RIGHT, neighbor, CHILD_1);
            if (neighbor->xmin() <= xc  &&  neighbor->zmax() >= zc) makeNeighbors(RIGHT, neighbor, CHILD_4);
            if (neighbor->xmax() >= xc  &&  neighbor->zmax() >= zc) makeNeighbors(RIGHT, neighbor, CHILD_5);
        }
    }
    // The RIGHT neighbors of this node
    {
        const vector<TreeNode*>& neighbors = _neighbors[RIGHT];
        for (unsigned int i=0; i<neighbors.size(); i++)
        {
            TreeNode* neighbor = neighbors[i];
            neighbor->deleteNeighbor(LEFT, this);
            if (neighbor->xmin() <= xc  &&  neighbor->zmin() <= zc) makeNeighbors(LEFT, neighbor, CHILD_2);
            if (neighbor->xmax() >= xc  &&  neighbor->zmin() <= zc) makeNeighbors(LEFT, neighbor, CHILD_3);
            if (neighbor->xmin() <= xc  &&  neighbor->zmax() >= zc) makeNeighbors(LEFT, neighbor, CHILD_6);
            if (neighbor->xmax() >= xc  &&  neighbor->zmax() >= zc) makeNeighbors(LEFT, neighbor, CHILD_7);
        }
    }
    // The BOTTOM neighbors of this node
    {
        const vector<TreeNode*>& neighbors = _neighbors[BOTTOM];
        for (unsigned int i=0; i<neighbors.size(); i++)
        {
            TreeNode* neighbor = neighbors[i];
            neighbor->deleteNeighbor(TOP, this);
            if (neighbor->xmin() <= xc  &&  neighbor->ymin() <= yc) makeNeighbors(TOP, neighbor, CHILD_0);
            if (neighbor->xmax() >= xc  &&  neighbor->ymin() <= yc) makeNeighbors(TOP, neighbor, CHILD_1);
            if (neighbor->xmin() <= xc  &&  neighbor->ymax() >= yc) makeNeighbors(TOP, neighbor, CHILD_2);
            if (neighbor->xmax() >= xc  &&  neighbor->ymax() >= yc) makeNeighbors(TOP, neighbor, CHILD_3);
        }
    }
    // The TOP neighbors of this node
    {
        const vector<TreeNode*>& neighbors = _neighbors[TOP];
        for (unsigned int i=0; i<neighbors.size(); i++)
        {
            TreeNode* neighbor = neighbors[i];
            neighbor->deleteNeighbor(BOTTOM, this);
            if (neighbor->xmin() <= xc  &&  neighbor->ymin() <= yc) makeNeighbors(BOTTOM, neighbor, CHILD_4);
            if (neighbor->xmax() >= xc  &&  neighbor->ymin() <= yc) makeNeighbors(BOTTOM, neighbor, CHILD_5);
            if (neighbor->xmin() <= xc  &&  neighbor->ymax() >= yc) makeNeighbors(BOTTOM, neighbor, CHILD_6);
            if (neighbor->xmax() >= xc  &&  neighbor->ymax() >= yc) makeNeighbors(BOTTOM, neighbor, CHILD_7);
        }
    }
}

//////////////////////////////////////////////////////////////////////

TreeNode* OctTreeNode::child(Vec r) const
{
    Vec rc = CHILD_0->rmax();
    int l = (r.x()<rc.x() ? 0 : 1) + (r.y()<rc.y() ? 0 : 2) + (r.z()<rc.z() ? 0 : 4);
    return _children[l];
}

//////////////////////////////////////////////////////////////////////

