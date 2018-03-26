/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileTreeDustGrid.hpp"
#include "BinTreeNode.hpp"
#include "DustDistribution.hpp"
#include "DustGridPath.hpp"
#include "DustGridPlotFile.hpp"
#include "DustMassInBoxInterface.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "OctTreeNode.hpp"
#include "Random.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

FileTreeDustGrid::~FileTreeDustGrid()
{
    for (int l=0; l<_Nnodes; l++)
        delete _tree[l];
}

////////////////////////////////////////////////////////////////////

void FileTreeDustGrid::setupSelfBefore()
{
    DustGrid::setupSelfBefore();

    // open the file representing the tree
    TextInFile infile(this, _filename, "dust grid tree");

    // extract the conversion factor for the coordinates
    string line = infile.readHeaderLine("coordinate");
    auto i1 = line.find('(');
    auto i2 = line.find(')');
    if (i1==string::npos || i2==string::npos || i2<=i1+1)
        throw FATALERROR("Dust grid tree file does not specify coordinate units");
    string unit = line.substr(i1+1, i2-i1-1);
    double factor = find<Units>()->in("length", unit, 1.);

    // initialize temporary storage
    enum { Unknown, BinTree, OctTree } treetype = Unknown;
    vector<vector<int>> children;
    Array row;

    // read the tree nodes from the input file
    while (infile.readRow(row, 17))
    {
        // extract node ID and dust cell index
        int id = static_cast<int>(row[0]);
        int m = static_cast<int>(row[1]);

        // extract IDs of the father and child nodes
        int father = static_cast<int>(row[8]);
        int child0 = static_cast<int>(row[9]);
        int child1 = static_cast<int>(row[10]);
        int child2 = static_cast<int>(row[11]);
        int child3 = static_cast<int>(row[12]);
        int child4 = static_cast<int>(row[13]);
        int child5 = static_cast<int>(row[14]);
        int child6 = static_cast<int>(row[15]);
        int child7 = static_cast<int>(row[16]);

        // extract node extent (note: order differs)
        Box extent(row[2]*factor, row[4]*factor, row[6]*factor, row[3]*factor, row[5]*factor, row[7]*factor);

        // figure out the tree type from the root node (the first node in the list)
        if (treetype == Unknown)
        {
            if (child0 == -1) throw FATALERROR("Root node in dust grid tree file has no children");
            treetype = (child2 == -1) ? BinTree : OctTree;
        }

        // remember children IDs
        if (child0 == -1)
        {
            children.emplace_back();  // if there are no children, add empty list
        }
        else
        {
            if (treetype == BinTree) children.emplace_back(std::initializer_list<int>{child0, child1});
            else children.emplace_back(std::initializer_list<int>{child0, child1, child2, child3,
                                                                  child4, child5, child6, child7});
        }

        // retrieve a pointer to the father node, which has already been created
        TreeNode* fatherNode = father >= 0 ? _tree[father] : nullptr;

        // create and add a new node of the appropriate type
        if (treetype == BinTree) _tree.push_back(new BinTreeNode(fatherNode, id, extent));
        else _tree.push_back(new OctTreeNode(fatherNode, id, extent));

        // add the cell index for the node
        _cellnumberv.push_back(m);
    }

    // remember the number of nodes and the extent of the complete domain
    _Nnodes = _tree.size();
    static_cast<Box>(*this) = root()->extent();
    _eps = 1e-12 * extent().widths().norm();

    // add children to the nodes
    for (int l=0; l<_Nnodes; l++)
    {
        auto node = _tree[l];
        for (auto child : children[l]) node->addChild(_tree[child]);
    }

    // construct a vector that contains the node IDs of all leaf nodes (i.e. actual dust cells)
    for (int l=0; l<_Nnodes; l++)
    {
        if (_tree[l]->isChildless()) _idv.push_back(l);
    }
    int Ncells = _idv.size();

    // log the number of nodes and cells
    Log* log = find<Log>();
    log->info("Tree import has finished.");
    log->info("  Total number of nodes: " + std::to_string(_Nnodes));
    log->info("  Total number of leaves: " + std::to_string(Ncells));

    // log the structure of the tree
    vector<int> countv;
    int maxLevel = 0;
    for (int m=0; m<Ncells; m++)
    {
        auto node = _tree[_idv[m]];
        int level = node->level();
        if (level > maxLevel)
        {
            maxLevel = level;
            countv.resize(maxLevel+1);
        }
        countv[level]++;
    }
    log->info("  Number of leaf cells of each level:");
    for (int level=0; level<=maxLevel; level++)
        log->info("    Level " + std::to_string(level) + ": " + std::to_string(countv[level]) + " cells");

    // determine the number of levels to be included in 3D grid output (if such output is requested)
    if (writeGrid())
    {
        int cumulativeCells = 0;
        for (_highestWriteLevel=0; _highestWriteLevel<=maxLevel; _highestWriteLevel++)
        {
            cumulativeCells += countv[_highestWriteLevel];
            if (cumulativeCells > 1500) break;          // experimental number
        }
        if (_highestWriteLevel<maxLevel)
            log->info("Will be outputting 3D grid data up to level " + std::to_string(_highestWriteLevel) +
                      ", i.e. " + std::to_string(cumulativeCells) + " cells.");
    }

    // add neighbors to the tree structure (if required for the search method)
    if (_searchMethod == SearchMethod::Neighbor)
    {
        log->info("Adding neighbors to the tree nodes...");
        for (int l=0; l<_Nnodes; l++) _tree[l]->addNeighbors();
        for (int l=0; l<_Nnodes; l++) _tree[l]->sortNeighbors();
    }

    // Verify that bookkeeping method is not used with binary tree
    if (searchMethod() == SearchMethod::Bookkeeping && treetype == BinTree)
        throw FATALERROR("Bookkeeping method is not compatible with binary tree");

    // Cache pointers used elsewhere
    _random = find<Random>();
    _dmib = find<DustDistribution>()->interface<DustMassInBoxInterface>();
}

////////////////////////////////////////////////////////////////////

int FileTreeDustGrid::dimension() const
{
    return 3;
}

////////////////////////////////////////////////////////////////////

Box FileTreeDustGrid::boundingBox() const
{
    return extent();
}

////////////////////////////////////////////////////////////////////

double FileTreeDustGrid::volume(int m) const
{
    if (m<0 || m>numCells())
        throw FATALERROR("Invalid cell number: " + std::to_string(m));
    TreeNode* node = getNode(m);
    return node->xwidth() * node->ywidth() * node->zwidth();
}

////////////////////////////////////////////////////////////////////

int FileTreeDustGrid::numCells() const
{
    return _idv.size();
}

////////////////////////////////////////////////////////////////////

int FileTreeDustGrid::whichCell(Position bfr) const
{
    const TreeNode* node = root()->whichNode(bfr);
    return node ? cellNumber(node) : -1;
}

////////////////////////////////////////////////////////////////////

Position FileTreeDustGrid::centralPositionInCell(int m) const
{
    return Position(getNode(m)->extent().center());
}

////////////////////////////////////////////////////////////////////

Position FileTreeDustGrid::randomPositionInCell(int m) const
{
    return _random->position(getNode(m)->extent());
}

////////////////////////////////////////////////////////////////////

void FileTreeDustGrid::path(DustGridPath* path) const
{
    // Initialize the path
    path->clear();

    // If the photon package starts outside the dust grid, move it into the first grid cell that it will pass
    Position bfr = path->moveInside(extent(), _eps);

    // Get the node containing the current location;
    // if the position is not inside the grid, return an empty path
    const TreeNode* node = root()->whichNode(bfr);
    if (!node) return path->clear();

    // Start the loop over nodes/path segments until we leave the grid.
    // Use a different code segment depending on the search method.
    double x,y,z;
    bfr.cartesian(x,y,z);
    double kx,ky,kz;
    path->direction().cartesian(kx,ky,kz);

    // ----------- Top-down -----------

    if (_searchMethod == SearchMethod::TopDown)
    {
        while (node)
        {
            double xnext = (kx<0.0) ? node->xmin() : node->xmax();
            double ynext = (ky<0.0) ? node->ymin() : node->ymax();
            double znext = (kz<0.0) ? node->zmin() : node->zmax();
            double dsx = (fabs(kx)>1e-15) ? (xnext-x)/kx : DBL_MAX;
            double dsy = (fabs(ky)>1e-15) ? (ynext-y)/ky : DBL_MAX;
            double dsz = (fabs(kz)>1e-15) ? (znext-z)/kz : DBL_MAX;

            double ds;
            if (dsx<=dsy && dsx<=dsz) ds = dsx;
            else if (dsy<=dsx && dsy<=dsz) ds = dsy;
            else ds = dsz;
            path->addSegment(cellNumber(node), ds);
            x += (ds+_eps)*kx;
            y += (ds+_eps)*ky;
            z += (ds+_eps)*kz;

            // always search from the root node down
            const TreeNode* oldnode = node;
            node = root()->whichNode(Vec(x,y,z));

            // if we're stuck in the same node...
            if (node==oldnode)
            {
                // try to escape by advancing the position to the next representable coordinates
                find<Log>()->warning("Photon package seems stuck in dust cell "
                                     + std::to_string(node->id()) + " -- escaping");
                x = nextafter(x, (kx<0.0) ? -DBL_MAX : DBL_MAX);
                y = nextafter(y, (ky<0.0) ? -DBL_MAX : DBL_MAX);
                z = nextafter(z, (kz<0.0) ? -DBL_MAX : DBL_MAX);
                node = root()->whichNode(Vec(x,y,z));

                // if that didn't work, terminate the path
                if (node==oldnode)
                {
                    find<Log>()->warning("Photon package is stuck in dust cell "
                                         + std::to_string(node->id()) + " -- terminating this path");
                    break;
                }
            }
        }
    }

    // ----------- Neighbor -----------

    else if (_searchMethod == SearchMethod::Neighbor)
    {
        while (node)
        {
            double xnext = (kx<0.0) ? node->xmin() : node->xmax();
            double ynext = (ky<0.0) ? node->ymin() : node->ymax();
            double znext = (kz<0.0) ? node->zmin() : node->zmax();
            double dsx = (fabs(kx)>1e-15) ? (xnext-x)/kx : DBL_MAX;
            double dsy = (fabs(ky)>1e-15) ? (ynext-y)/ky : DBL_MAX;
            double dsz = (fabs(kz)>1e-15) ? (znext-z)/kz : DBL_MAX;

            double ds;
            TreeNode::Wall wall;
            if (dsx<=dsy && dsx<=dsz)
            {
                ds = dsx;
                wall = (kx<0.0) ? TreeNode::BACK : TreeNode::FRONT;
            }
            else if (dsy<=dsx && dsy<=dsz)
            {
                ds = dsy;
                wall = (ky<0.0) ? TreeNode::LEFT : TreeNode::RIGHT;
            }
            else
            {
                ds = dsz;
                wall = (kz<0.0) ? TreeNode::BOTTOM : TreeNode::TOP;
            }
            path->addSegment(cellNumber(node), ds);
            x += (ds+_eps)*kx;
            y += (ds+_eps)*ky;
            z += (ds+_eps)*kz;

            // attempt to find the new node among the neighbors of the current node;
            // this should not fail unless the new location is outside the grid,
            // however on rare occasions it fails due to rounding errors (e.g. in a corner),
            // thus we use top-down search as a fall-back
            const TreeNode* oldnode = node;
            node = node->whichNode(wall, Vec(x,y,z));
            if (!node) node = root()->whichNode(Vec(x,y,z));

            // if we're stuck in the same node...
            if (node==oldnode)
            {
                // try to escape by advancing the position to the next representable coordinates
                find<Log>()->warning("Photon package seems stuck in dust cell "
                                     + std::to_string(node->id()) + " -- escaping");
                x = nextafter(x, (kx<0.0) ? -DBL_MAX : DBL_MAX);
                y = nextafter(y, (ky<0.0) ? -DBL_MAX : DBL_MAX);
                z = nextafter(z, (kz<0.0) ? -DBL_MAX : DBL_MAX);
                node = root()->whichNode(Vec(x,y,z));

                // if that didn't work, terminate the path
                if (node==oldnode)
                {
                    find<Log>()->warning("Photon package is stuck in dust cell "
                                         + std::to_string(node->id()) + " -- terminating this path");
                    break;
                }
            }
        }
    }

    // ----------- Bookkeeping -----------

    // !! This code section relies on the fact that an octtree node is used !!

    else if (_searchMethod == SearchMethod::Bookkeeping)
    {
        int l = node->id();  // node index in the tree
        while (true)
        {
            double xnext = (kx<0.0) ? _tree[l]->xmin() : _tree[l]->xmax();
            double ynext = (ky<0.0) ? _tree[l]->ymin() : _tree[l]->ymax();
            double znext = (kz<0.0) ? _tree[l]->zmin() : _tree[l]->zmax();
            double dsx = (fabs(kx)>1e-15) ? (xnext-x)/kx : DBL_MAX;
            double dsy = (fabs(ky)>1e-15) ? (ynext-y)/ky : DBL_MAX;
            double dsz = (fabs(kz)>1e-15) ? (znext-z)/kz : DBL_MAX;

            // First option: the x-wall is hit first. After moving
            // towards the boundary, we have to find the next cell. First we
            // check whether the node is on the right or left side of his
            // father node. If the movement is towards positive x (i.e. if
            // kx>0) we move up in the tree until we find a node on the left
            // side. The next cell will then be the corresponding right node
            // (if it is a leaf) or one of its children. If we have to move
            // up until we hit the root node, this means our path has ended.

            if (dsx<=dsy && dsx<=dsz)
            {
                path->addSegment(_cellnumberv[l], dsx);
                x = xnext;
                y += ky*dsx;
                z += kz*dsx;
                while (true)
                {
                    int oct = ((l - 1) % 8) + 1;
                    bool place = (kx<0.0) ? (oct % 2 == 1) : (oct % 2 == 0);
                    if (!place) break;
                    l = _tree[l]->father()->id();
                    if (l == 0) return;
                }
                l += (kx<0.0) ? -1 : 1;
                while (_cellnumberv[l] == -1)
                {
                    double yM = _tree[l]->child(0)->ymax();
                    double zM = _tree[l]->child(0)->zmax();
                    if (kx<0.0)
                    {
                        if (y<=yM)
                            l = (z<=zM) ? _tree[l]->child(1)->id() : _tree[l]->child(5)->id();
                        else
                            l = (z<=zM) ? _tree[l]->child(3)->id() : _tree[l]->child(7)->id();
                    }
                    else
                    {
                        if (y<=yM)
                            l = (z<=zM) ? _tree[l]->child(0)->id() : _tree[l]->child(4)->id();
                        else
                            l = (z<=zM) ? _tree[l]->child(2)->id() : _tree[l]->child(6)->id();
                    }
                }
            }

            // Repeat the same exercise, but now the y-wall is hit first...

            else if (dsy<dsx && dsy<=dsz)
            {
                path->addSegment(_cellnumberv[l], dsy);
                x += kx*dsy;
                y  = ynext;
                z += kz*dsy;
                while (true)
                {
                    bool place = (ky<0.0) ? ((l-1) % 4 < 2) : ((l-1) % 4 > 1);
                    if (!place) break;
                    l = _tree[l]->father()->id();
                    if (l == 0) return;
                }
                l += (ky<0.0) ? -2 : 2;
                while (_cellnumberv[l] == -1)
                {
                    double xM = _tree[l]->child(0)->xmax();
                    double zM = _tree[l]->child(0)->zmax();
                    if (ky<0.0)
                    {
                        if (x<=xM)
                            l = (z<=zM) ? _tree[l]->child(2)->id() : _tree[l]->child(6)->id();
                        else
                            l = (z<=zM) ? _tree[l]->child(3)->id() : _tree[l]->child(7)->id();
                    }
                    else
                    {
                        if (x<=xM)
                            l = (z<=zM) ? _tree[l]->child(0)->id() : _tree[l]->child(4)->id();
                        else
                            l = (z<=zM) ? _tree[l]->child(1)->id() : _tree[l]->child(5)->id();
                    }
                }
            }

            // Finally, repeat the same exercise, but now the z-wall is hit first...

            else if (dsz< dsx && dsz< dsy)
            {
                path->addSegment(_cellnumberv[l], dsz);
                x += kx*dsz;
                y += ky*dsz;
                z  = znext;
                while (true)
                {
                    int oct = ((l-1) % 8) + 1;
                    bool place = (kz<0.0) ? (oct < 5) : (oct > 4);
                    if (!place) break;
                    l = _tree[l]->father()->id();
                    if (l == 0) return;
                }
                l += (kz<0.0) ? -4 : 4;
                while (_cellnumberv[l] == -1)
                {
                    double xM = _tree[l]->child(0)->xmax();
                    double yM = _tree[l]->child(0)->ymax();
                    if (kz<0.0)
                    {
                        if (x<=xM)
                            l = (y<=yM) ? _tree[l]->child(4)->id() : _tree[l]->child(6)->id();
                        else
                            l = (y<=yM) ? _tree[l]->child(5)->id() : _tree[l]->child(7)->id();
                    }
                    else
                    {
                        if (x<=xM)
                            l = (y<=yM) ? _tree[l]->child(0)->id() : _tree[l]->child(2)->id();
                        else
                            l = (y<=yM) ? _tree[l]->child(1)->id() : _tree[l]->child(3)->id();
                    }
                }
            }
        }
    }

    // ------------------------------
}

////////////////////////////////////////////////////////////////////

vector<SimulationItem*> FileTreeDustGrid::interfaceCandidates(const std::type_info& interfaceTypeInfo)
{
    if (interfaceTypeInfo == typeid(DustGridDensityInterface) && !_dmib)
        return vector<SimulationItem*>();
    return DustGrid::interfaceCandidates(interfaceTypeInfo);
}

////////////////////////////////////////////////////////////////////

double FileTreeDustGrid::density(int h, int m) const
{
    TreeNode* node = getNode(m);
    return _dmib->massInBox(h, node->extent()) / node->volume();
}

////////////////////////////////////////////////////////////////////

void FileTreeDustGrid::write_xy(DustGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(_xmin, _ymin, _xmax, _ymax);
    int Ncells = numCells();
    for (int m=0; m<Ncells; m++)
    {
        TreeNode* node = getNode(m);
        if (fabs(node->zmin()) < 1e-8*extent().zwidth())
        {
            outfile->writeRectangle(node->xmin(), node->ymin(), node->xmax(), node->ymax());
        }
    }
}

////////////////////////////////////////////////////////////////////

void FileTreeDustGrid::write_xz(DustGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(_xmin, _zmin, _xmax, _zmax);
    int Ncells = numCells();
    for (int m=0; m<Ncells; m++)
    {
        TreeNode* node = getNode(m);
        if (fabs(node->ymin()) < 1e-8*extent().ywidth())
        {
            outfile->writeRectangle(node->xmin(), node->zmin(), node->xmax(), node->zmax());
        }
    }
}

////////////////////////////////////////////////////////////////////

void FileTreeDustGrid::write_yz(DustGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(_ymin, _zmin, _ymax, _zmax);
    int Ncells = numCells();
    for (int m=0; m<Ncells; m++)
    {
        TreeNode* node = getNode(m);
        if (fabs(node->xmin()) < 1e-8*extent().xwidth())
        {
            outfile->writeRectangle(node->ymin(), node->zmin(), node->ymax(), node->zmax());
        }
    }
}

////////////////////////////////////////////////////////////////////

void FileTreeDustGrid::write_xyz(DustGridPlotFile* outfile) const
{
    // Output all leaf cells up to a certain level
    int Ncells = numCells();
    for (int m=0; m<Ncells; m++)
    {
        TreeNode* node = getNode(m);
        if (node->level() <= _highestWriteLevel)
            outfile->writeCube(node->xmin(), node->ymin(), node->zmin(), node->xmax(), node->ymax(), node->zmax());
    }
}

////////////////////////////////////////////////////////////////////

TreeNode* FileTreeDustGrid::root() const
{
    return _tree[0];
}

////////////////////////////////////////////////////////////////////

TreeNode* FileTreeDustGrid::getNode(int m) const
{
    return _tree[_idv[m]];
}

////////////////////////////////////////////////////////////////////

int FileTreeDustGrid::cellNumber(const TreeNode* node) const
{
    return _cellnumberv[node->id()];
}

////////////////////////////////////////////////////////////////////
