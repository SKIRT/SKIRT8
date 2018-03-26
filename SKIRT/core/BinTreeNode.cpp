/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BinTreeNode.hpp"
#include "FatalError.hpp"

//////////////////////////////////////////////////////////////////////

// macros for easily accessing a particular child
#define CHILD_0 _children[0]
#define CHILD_1 _children[1]

//////////////////////////////////////////////////////////////////////

// constants for splitting direction (not an enum because we use modulo arithmetic)
#define XDIR 0
#define YDIR 1
#define ZDIR 2

//////////////////////////////////////////////////////////////////////

BinTreeNode::BinTreeNode(TreeNode* father, int id, const Box& extent)
    : TreeNode(father, id, extent), _dir(level() % 3)
{
}

//////////////////////////////////////////////////////////////////////

TreeNode* BinTreeNode::createNode(TreeNode* father, int id, const Box& extent)
{
    return new BinTreeNode(father, id, extent);
}

//////////////////////////////////////////////////////////////////////

void BinTreeNode::createChildrenSplitDir(int id, int dir)
{
    _children.resize(2);
    _dir = dir;
    switch (dir)
    {
    case XDIR:
        {
            double xc = 0.5*(_xmin+_xmax);
            CHILD_0 = createNode(this, id++, Box(_xmin, _ymin, _zmin,    xc, _ymax, _zmax));
            CHILD_1 = createNode(this, id++, Box(   xc, _ymin, _zmin, _xmax, _ymax, _zmax));
        }
        break;
    case YDIR:
        {
            double yc = 0.5*(_ymin+_ymax);
            CHILD_0 = createNode(this, id++, Box(_xmin, _ymin, _zmin, _xmax,    yc, _zmax));
            CHILD_1 = createNode(this, id++, Box(_xmin,    yc, _zmin, _xmax, _ymax, _zmax));
        }
        break;
    case ZDIR:
        {
            double zc = 0.5*(_zmin+_zmax);
            CHILD_0 = createNode(this, id++, Box(_xmin, _ymin, _zmin, _xmax, _ymax,    zc));
            CHILD_1 = createNode(this, id++, Box(_xmin, _ymin,    zc, _xmax, _ymax, _zmax));
        }
        break;
    default:
        throw FATALERROR("Incorrect value for subdivision direction");
    }
}

//////////////////////////////////////////////////////////////////////

void BinTreeNode::createChildren(int id)
{
    createChildrenSplitDir(id, level() % 3);
}

//////////////////////////////////////////////////////////////////////

void BinTreeNode::createChildren(int id, const TreeNodeDensityCalculator* /*calc*/)
{
    createChildren(id);
}

//////////////////////////////////////////////////////////////////////

void BinTreeNode::addNeighbors()
{
    // if we don't have any children, we can't add neighbors
    if (_children.empty()) return;

    // ensure that all involved nodes have a neighbor list for each of the walls
    ensureNeighborLists();
    CHILD_0->ensureNeighborLists();
    CHILD_1->ensureNeighborLists();

    switch (_dir)
    {
    case XDIR:
        {
            double xc = CHILD_0->xmax();

            // Internal neighbors
            makeNeighbors(FRONT, CHILD_0, CHILD_1);

            // The BACK neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[BACK];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(FRONT, this);
                    makeNeighbors(FRONT, neighbor, CHILD_0);
                }
            }
            // The FRONT neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[FRONT];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(BACK, this);
                    makeNeighbors(BACK, neighbor, CHILD_1);
                }
            }
            // The LEFT neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[LEFT];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(RIGHT, this);
                    if (neighbor->xmin() <= xc) makeNeighbors(RIGHT, neighbor, CHILD_0);
                    if (neighbor->xmax() >= xc) makeNeighbors(RIGHT, neighbor, CHILD_1);
                }
            }
            // The RIGHT neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[RIGHT];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(LEFT, this);
                    if (neighbor->xmin() <= xc) makeNeighbors(LEFT, neighbor, CHILD_0);
                    if (neighbor->xmax() >= xc) makeNeighbors(LEFT, neighbor, CHILD_1);
                }
            }
            // The BOTTOM neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[BOTTOM];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(TOP, this);
                    if (neighbor->xmin() <= xc) makeNeighbors(TOP, neighbor, CHILD_0);
                    if (neighbor->xmax() >= xc) makeNeighbors(TOP, neighbor, CHILD_1);
                }
            }
            // The TOP neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[TOP];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(BOTTOM, this);
                    if (neighbor->xmin() <= xc) makeNeighbors(BOTTOM, neighbor, CHILD_0);
                    if (neighbor->xmax() >= xc) makeNeighbors(BOTTOM, neighbor, CHILD_1);
                }
            }
        }
        break;
    case YDIR:
        {
            double yc = CHILD_0->ymax();

            // Internal neighbors
            makeNeighbors(RIGHT, CHILD_0, CHILD_1);

            // The BACK neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[BACK];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(FRONT, this);
                    if (neighbor->ymin() <= yc) makeNeighbors(FRONT, neighbor, CHILD_0);
                    if (neighbor->ymax() >= yc) makeNeighbors(FRONT, neighbor, CHILD_1);
                }
            }
            // The FRONT neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[FRONT];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(BACK, this);
                    if (neighbor->ymin() <= yc) makeNeighbors(BACK, neighbor, CHILD_0);
                    if (neighbor->ymax() >= yc) makeNeighbors(BACK, neighbor, CHILD_1);
                }
            }
            // The LEFT neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[LEFT];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(RIGHT, this);
                    makeNeighbors(RIGHT, neighbor, CHILD_0);
                }
            }
            // The RIGHT neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[RIGHT];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(LEFT, this);
                    makeNeighbors(LEFT, neighbor, CHILD_1);
                }
            }
            // The BOTTOM neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[BOTTOM];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(TOP, this);
                    if (neighbor->ymin() <= yc) makeNeighbors(TOP, neighbor, CHILD_0);
                    if (neighbor->ymax() >= yc) makeNeighbors(TOP, neighbor, CHILD_1);
                }
            }
            // The TOP neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[TOP];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(BOTTOM, this);
                    if (neighbor->ymin() <= yc) makeNeighbors(BOTTOM, neighbor, CHILD_0);
                    if (neighbor->ymax() >= yc) makeNeighbors(BOTTOM, neighbor, CHILD_1);
                }
            }
        }
        break;
    case ZDIR:
        {
            double zc = CHILD_0->zmax();

            // Internal neighbors
            makeNeighbors(TOP, CHILD_0, CHILD_1);

            // The BACK neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[BACK];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(FRONT, this);
                    if (neighbor->zmin() <= zc) makeNeighbors(FRONT, neighbor, CHILD_0);
                    if (neighbor->zmax() >= zc) makeNeighbors(FRONT, neighbor, CHILD_1);
                }
            }
            // The FRONT neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[FRONT];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(BACK, this);
                    if (neighbor->zmin() <= zc) makeNeighbors(BACK, neighbor, CHILD_0);
                    if (neighbor->zmax() >= zc) makeNeighbors(BACK, neighbor, CHILD_1);
                }
            }
            // The LEFT neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[LEFT];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(RIGHT, this);
                    if (neighbor->zmin() <= zc) makeNeighbors(RIGHT, neighbor, CHILD_0);
                    if (neighbor->zmax() >= zc) makeNeighbors(RIGHT, neighbor, CHILD_1);
                }
            }
            // The RIGHT neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[RIGHT];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(LEFT, this);
                    if (neighbor->zmin() <= zc) makeNeighbors(LEFT, neighbor, CHILD_0);
                    if (neighbor->zmax() >= zc) makeNeighbors(LEFT, neighbor, CHILD_1);
                }
            }
            // The BOTTOM neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[BOTTOM];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(TOP, this);
                    makeNeighbors(TOP, neighbor, CHILD_0);
                }
            }
            // The TOP neighbors of this node
            {
                const vector<TreeNode*>& neighbors = _neighbors[TOP];
                for (unsigned int i=0; i<neighbors.size(); i++)
                {
                    TreeNode* neighbor = neighbors[i];
                    neighbor->deleteNeighbor(BOTTOM, this);
                    makeNeighbors(BOTTOM, neighbor, CHILD_1);
                }
            }
        }
        break;
    default:
        throw FATALERROR("Incorrect value for subdivision direction");
    }
}

//////////////////////////////////////////////////////////////////////

TreeNode* BinTreeNode::child(Vec r) const
{
    switch (_dir)
    {
    case XDIR:  return  r.x() < CHILD_0->xmax()  ?  CHILD_0 : CHILD_1;
    case YDIR:  return  r.y() < CHILD_0->ymax()  ?  CHILD_0 : CHILD_1;
    case ZDIR:  return  r.z() < CHILD_0->zmax()  ?  CHILD_0 : CHILD_1;
    }
    throw FATALERROR("Incorrect value for subdivision direction");
}

//////////////////////////////////////////////////////////////////////
