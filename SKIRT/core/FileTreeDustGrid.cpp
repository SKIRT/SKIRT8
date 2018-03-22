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
    Log* log = find<Log>();

    // open the file representing the tree
    TextInFile infile(this, _filename, "dust grid tree");

    // extract the conversion factor for the coordinates
    string line = infile.readHeaderLine("coordinate");
    auto i1 = line.find('(');
    auto i2 = line.find(')');
    if (i1==string::npos || i2==string::npos || i2<=i1+1)
        throw FATALERROR("Dust grid tree file does not specify coordinate units");
    string unit = line.substr(i1+1, i2-i1-1);

    log->warning(unit);
    //double factor = 1.;
    int treeType = 0;
    // Load the tree
    /*

    // Initialize vectors
    std::vector<int> fathers;
    std::vector<std::vector<int>> children;

    // Inform the user
    log->info("Importing nodes...");

    // load the SPH gas particles
    TextInFile infile(this, _filename, "dust grid tree");
    double id = -1; // node ID
    double m = -1; // dust cell index
    double xmin = 0;  // minimum x coordinate of the node
    double xmax = 0;  // maximum x coordinate of the node
    double ymin = 0;  // minimum y coordinate of the node
    double ymax = 0;  // maximum y coordinate of the node
    double zmin = 0;  // minimum z coordinate of the node
    double zmax = 0;  // maximum z coordinate of the node
    double father = 0;  // ID of the father node
    double child0 = 0; // file.addColumn("ID of child node 0", 'd');
    double child1 = 0; // file.addColumn("ID of child node 1", 'd');
    double child2 = 0; // file.addColumn("ID of child node 2", 'd');
    double child3 = 0; // file.addColumn("ID of child node 3", 'd');
    double child4 = 0; // file.addColumn("ID of child node 4", 'd');
    double child5 = 0; // file.addColumn("ID of child node 5", 'd');
    double child6 = 0; // file.addColumn("ID of child node 6", 'd');
    double child7 = 0;// file.addColumn("ID of child node 7", 'd');

    int counter = 0;
    auto lines = infile.readHeader();
    cout << "file header:" << endl;
    for (size_t index = 0; index<lines.size(); index++)
    {
        cout << lines[index] << endl;
    }

    while (infile.readRow(0, id, m, xmin, xmax, ymin, ymax, zmin, zmax, father, child0, child1, child2, child3, child4, child5, child6, child7))
    {
        //counter++;
        //cout << counter << endl;

        // Cast
        int idi = (int)round(id);
        int mi = (int)round(m);
        int fatheri = (int)round(father);
        int child0i = (int)round(child0);
        int child1i = (int)round(child1);
        int child2i = (int)round(child2);
        int child3i = (int)round(child3);
        int child4i = (int)round(child4);
        int child5i = (int)round(child5);
        int child6i = (int)round(child6);
        int child7i = (int)round(child7);

        // CHECK ID
        if (idi != counter) throw FATALERROR("Tree file is invalid: id index " + QString::number(idi) + " does not match the position of that node in the list");

        // Increment counter
        counter++;

        // Remember father ID
        fathers.push_back(fatheri);

        // Remember father and children IDs
        if (child2i == -1)
        {
            // if child0i == -1, there are no children: add empty list
            if (child0i == -1)
            {
                // we cannot say anything about the node type
                children.push_back(std::vector<int>());
            }
            else
            {
                _treetype = "bintree";
                auto childreni = std::vector<int>();
                childreni.push_back(child0i);
                childreni.push_back(child1i);

                // Add children IDs
                children.push_back(childreni);
            }
        }
        else
        {
            _treetype = "octtree";
            auto childreni = std::vector<int>();
            childreni.push_back(child0i);
            childreni.push_back(child1i);
            childreni.push_back(child2i);
            childreni.push_back(child3i);
            childreni.push_back(child4i);
            childreni.push_back(child5i);
            childreni.push_back(child6i);
            childreni.push_back(child7i);

            // Add children IDs
            children.push_back(childreni);
        }

        // Create extent
        Box extent = Box(xmin*pc, ymin*pc, zmin*pc, xmax*pc, ymax*pc, zmax*pc);

        // father must be pointer!

        // Create new oct tree node
        TreeNode* node = 0;
        //if (child2i == -1)
        if (_treetype == "bintree")
        {
            node = new BinTreeNode(0, id, extent);
        }
        else if (_treetype == "octtree")
        {
            node = new OctTreeNode(0, id, extent);
        }
        else throw FATALERROR("Tree type not set: " + _treetype);

        // Add the new node
        //_tree.insert(_tree.end(), node->children().begin(), node->children().end());

        // Get cell number
        //int m = _cellnumberv[l];
        //_cellnumberv[idi] = m;
        _cellnumberv.push_back(mi);

        // Add the new node
        _tree.push_back(node);
    }

    // Inform the user
    log->info("Setting fathers and children...");

    // Set fathers and children
    for (size_t l=0; l<fathers.size(); l++)
    {
        //cout << l << endl;
        //cout << "  father: " << fathers[l] << endl;

        // Set father
        if (fathers[l] != -1)
        {
            // Get the father pointer
            auto father = _tree[fathers[l]];

            // Set the father pointer to the tree node
            _tree[l]->setFather(father);
            // NOW ALSO AUTOMATICALLY SETS LEVEL: NO, ITERATE OVER THE NODES AGAIN TO SET THE LEVELS CORRECTLY

            // The constructor sets the level of the new node to be one higher than the level of the father.
            // If the pointer to the father is null, the level of the new cell is zero.

            // Set the level
            //int level = father->level() + 1;
            //_tree[l]->setLevel(level);

            //cout << level << endl;

            // Update maximum level
            //if (level > _maxlevel) _maxlevel = level;
        }
        else // no father
        {
            // Level must be zero, but this is guarenteed by the constructor of TreeNode when no father is passed
        }

        //cout << "here" << endl;
        //cout << "  number of children: " << children[l].size() << endl;

        // Loop over the children
        for (size_t k=0; k<children[l].size(); k++)
        {
            auto child = _tree[children[l][k]];
            _tree[l]->addChild(child);
        }
    }

    // Set children
    //for (size_t j=0; j<children.size(); j++)
    //{
    //    // Loop over the children
    //    for (size_t k=0; k<children[j].size(); k++)
    //    {
    //        auto child = _tree[children[j][k]];
    //        _tree[l]->addChild(child);
    //    }
    //}

    // Set the number of nodes and the number of cells
    _Nnodes = _tree.size();

    // Inform the user
    log->info("Setting tree node levels...");

    // Set levels
    for (int l=0; l<_Nnodes; l++)
    {
        int level = _tree[l]->nAncestors();
        _tree[l]->setLevel(level);

        // Update maximum level
        if (level > _maxlevel) _maxlevel = level;
    }

    // Only for binary trees
    if (_treetype == "bintree")
    {
        // Inform the user
        log->info("Setting the direction methods...");

        for (int l=0; l<_Nnodes; l++)
        {
            // Determine the direction
            int dir = _tree[l]->level() % 3;

            // CAST
            BinTreeNode* bintreenode = dynamic_cast<BinTreeNode*>(_tree[l]);

            // Set the direction
            bintreenode->setDir(dir);
        }
    }

    // Set the extent
    //_extent = extent;
    _extent = root()->extent();

    // Set EPS here

    // Construction of a vector _idv that contains the node IDs of all
    // leaves. This is the actual dust cell vector (only the leaves will
    // eventually become valid dust cells). We also create a vector
    // _cellnumberv with the cell numbers of all the nodes (i.e. the
    // rank m of the node in the vector _idv if the node is a leaf, and
    // -1 if the node is not a leaf).

    //int m = 0;
    //_cellnumberv.resize(_Nnodes,-1);
    for (int l=0; l<_Nnodes; l++)
    {
        if (_tree[l]->ynchildless())
        {
            _idv.push_back(l);
            //_cellnumberv[l] = m;
            //m++;
        }
    }
    //int Ncells = _idv.size();
    int Ncells = _idv.size();

    // Log the number of cells

    log->info("Construction of the tree finished.");
    log->info("  Total number of nodes: " + QString::number(_Nnodes));
    log->info("  Total number of leaves: " + QString::number(Ncells));

    vector<int> countv(_maxlevel+1);
    for (int m=0; m<Ncells; m++)
    {
        TreeNode* node = _tree[_idv[m]];
        int level = node->level();
        countv[level]++;
    }
    log->info("  Number of leaf cells of each level:");
    for (int level=0; level<=_maxlevel; level++)
        log->info("    Level " + QString::number(level) + ": " + QString::number(countv[level]) + " cells");

    // Determine the number of levels to be included in 3D grid output (if such output is requested)

    if (writeGrid())
    {
        int cumulativeCells = 0;
        for (_highestWriteLevel=0; _highestWriteLevel<=_maxlevel; _highestWriteLevel++)
        {
            cumulativeCells += countv[_highestWriteLevel];
            if (cumulativeCells > 1500) break;          // experimental number
        }
        if (_highestWriteLevel<_maxlevel)
            log->info("Will be outputting 3D grid data up to level " + QString::number(_highestWriteLevel) +
                      ", i.e. " + QString::number(cumulativeCells) + " cells.");
    }

    // Add neighbors to the tree structure (but only if required for the search method)

    if (_search == Neighbor)
    {
        log->info("Adding neighbors to the tree nodes...");
        for (int l=0; l<_Nnodes; l++) _tree[l]->addneighbors();
        for (int l=0; l<_Nnodes; l++) _tree[l]->sortneighbors();
    }

    // Determine the number of levels to be included in 3D grid output (if such output is requested)

    int maxLevel=9;
    vector<int> countv(maxLevel+1);

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

*/

    // Cache various pieces of information used elsewhere
    _Nnodes = _tree.size();
    static_cast<Box>(*this) = root()->extent();
    _eps = 1e-12 * extent().widths().norm();
    _random = find<Random>();
    _dmib = find<DustDistribution>()->interface<DustMassInBoxInterface>();

    // Add neighbors to the tree structure (but only if required for the search method)
    if (_searchMethod == SearchMethod::Neighbor)
    {
        log->info("Adding neighbors to the tree nodes...");
        for (int l=0; l<_Nnodes; l++) _tree[l]->addNeighbors();
        for (int l=0; l<_Nnodes; l++) _tree[l]->sortNeighbors();
    }

    // Verify that bookkeeping method is not used with binary tree
    if (searchMethod() == SearchMethod::Bookkeeping && treeType == 0)
        throw FATALERROR("Bookkeeping method is not compatible with binary tree");
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
