/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TreeDustGrid.hpp"
#include "DustDistribution.hpp"
#include "DustGridPath.hpp"
#include "DustGridPlotFile.hpp"
#include "DustMassInBoxInterface.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "Random.hpp"
#include "TreeNode.hpp"
#include "TreeNodeBoxDensityCalculator.hpp"
#include "TreeNodeSampleDensityCalculator.hpp"

//////////////////////////////////////////////////////////////////////

TreeDustGrid::~TreeDustGrid()
{
    for (int l=0; l<_Nnodes; l++)
        delete _tree[l];
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::setupSelfBefore()
{
    BoxDustGrid::setupSelfBefore();
    Log* log = find<Log>();

    // Validate attribute values
    if (_maxLevel <= _minLevel) throw FATALERROR("Maximum tree level should be larger than minimum tree level");

    // Cache some often used values
    // A Parallel instance is created with a limited amount of threads (4) for performance reasons
    // or with a single thread in case of multiprocessing, to ensure consistency.
    size_t nthreads = find<PeerToPeerCommunicator>()->isMultiProc() ? 1 : 4;
    _random = find<Random>();
    _parallel = find<ParallelFactory>()->parallel(nthreads);
    _dd = find<DustDistribution>();
    _dmib = _dd->interface<DustMassInBoxInterface>();
    _useDmibForSubdivide = _dmib && !_maxDensityDispersion && canUseDmibForSubdivide();
    _totalmass = _dd->mass();
    _eps = 1e-12 * extent().widths().norm();

    // Create the root node

    _tree.push_back(createRoot(extent()));

    // Recursively subdivide the root node until all nodes satisfy the
    // necessary criteria. When finished, set the number _Nnodes.

    int currentlevel = -1;
    unsigned int l = 0;
    while (true)
    {
        TreeNode* node = _tree[l];
        int level = node->level();
        if (level>currentlevel)
        {
            log->info("Starting subdivision of level " + std::to_string(level) + "...");
            currentlevel = level;
        }
        if (l%50000 == 0)
            log->info("Subdividing node number " + std::to_string(l) + "...");
        if (node->isChildless())
            subdivide(node);
        l++;
        if (l>=_tree.size()) break;
    }
    _Nnodes = _tree.size();

    // Construction of a vector _idv that contains the node IDs of all
    // leaves. This is the actual dust cell vector (only the leaves will
    // eventually become valid dust cells). We also create a vector
    // _cellnumberv with the cell numbers of all the nodes (i.e. the
    // rank m of the node in the vector _idv if the node is a leaf, and
    // -1 if the node is not a leaf).

    int m = 0;
    _cellnumberv.resize(_Nnodes,-1);
    for (int l=0; l<_Nnodes; l++)
    {
        if (_tree[l]->isChildless())
        {
            _idv.push_back(l);
            _cellnumberv[l] = m;
            m++;
        }
    }
    int Ncells = _idv.size();

    // Log the number of cells

    log->info("Construction of the tree finished.");
    log->info("  Total number of nodes: " + std::to_string(_Nnodes));
    log->info("  Total number of leaves: " + std::to_string(Ncells));
    vector<int> countv(_maxLevel+1);
    for (int m=0; m<Ncells; m++)
    {
        TreeNode* node = _tree[_idv[m]];
        int level = node->level();
        countv[level]++;
    }
    log->info("  Number of leaf cells of each level:");
    for (int level=0; level<=_maxLevel; level++)
        log->info("    Level " + std::to_string(level) + ": " + std::to_string(countv[level]) + " cells");

    // Determine the number of levels to be included in 3D grid output (if such output is requested)

    if (writeGrid())
    {
        int cumulativeCells = 0;
        for (_highestWriteLevel=0; _highestWriteLevel<=_maxLevel; _highestWriteLevel++)
        {
            cumulativeCells += countv[_highestWriteLevel];
            if (cumulativeCells > 1500) break;          // experimental number
        }
        if (_highestWriteLevel<_maxLevel)
            log->info("Will be outputting 3D grid data up to level " + std::to_string(_highestWriteLevel) +
                      ", i.e. " + std::to_string(cumulativeCells) + " cells.");
    }

    // Add neighbors to the tree structure (but only if required for the search method)

    if (_searchMethod == SearchMethod::Neighbor)
    {
        log->info("Adding neighbors to the tree nodes...");
        for (int l=0; l<_Nnodes; l++) _tree[l]->addNeighbors();
        for (int l=0; l<_Nnodes; l++) _tree[l]->sortNeighbors();
    }
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::subdivide(TreeNode* node)
{
    // If level is below or at minlevel, there is always subdivision, and the subdivision is "regular"
    int level = node->level();
    if (level <= _minLevel)
    {
        node->createChildren(_tree.size());
        _tree.insert(_tree.end(), node->children().begin(), node->children().end());
    }

    // if level is below maxlevel, there may be subdivision depending on various stopping criteria
    else if (level < _maxLevel)
    {
        // construct an appropriate density calculator to estimate properties for stopping criteria and division
        TreeNodeDensityCalculator* calc;
        if (_useDmibForSubdivide)
        {
            // use the DustMassInBox interface
            calc = new TreeNodeBoxDensityCalculator(_dmib, node);
        }
        else
        {
            // sample the density in the cell
            TreeNodeSampleDensityCalculator* sampleCalc =
                    new TreeNodeSampleDensityCalculator(_random, _numSamples, _dd, node);
            _parallel->call(sampleCalc, _numSamples);
            calc = sampleCalc;
        }

        // if no stopping criteria are enabled, we keep subdividing indefinitely
        bool needDivision = (_maxOpticalDepth == 0 && _maxMassFraction == 0 && _maxDensityDispersion == 0);

        // otherwise, we subdivide if at least one stopping criterion is not satisfied

        // check mass fraction
        if (!needDivision && _maxMassFraction > 0)
        {
            double massfraction = calc->mass() / _totalmass;
            if (massfraction >= _maxMassFraction) needDivision = true;
        }

        // check optical depth
        if (!needDivision && _maxOpticalDepth > 0)
        {
            double opticaldepth = calc->opticalDepth();
            if (opticaldepth >= _maxOpticalDepth) needDivision = true;
        }

        // check density dispersion fraction
        if (!needDivision && _maxDensityDispersion > 0)
        {
            double densdispfraction = calc->densityDispersion();
            if (densdispfraction >= _maxDensityDispersion) needDivision = true;
        }

        if (needDivision)
        {
            // there is subdivision, possibly using calculated properties such as barycenter
            node->createChildren(_tree.size(), calc);
            _tree.insert(_tree.end(), node->children().begin(), node->children().end());
        }

        delete calc;
    }
}

////////////////////////////////////////////////////////////////////

double TreeDustGrid::volume(int m) const
{
    if (m<0 || m>numCells())
        throw FATALERROR("Invalid cell number: " + std::to_string(m));
    TreeNode* node = getNode(m);
    return node->xwidth() * node->ywidth() * node->zwidth();
}

//////////////////////////////////////////////////////////////////////

int TreeDustGrid::numCells() const
{
    return _idv.size();
}

//////////////////////////////////////////////////////////////////////

int TreeDustGrid::whichCell(Position bfr) const
{
    const TreeNode* node = root()->whichNode(bfr);
    return node ? cellNumber(node) : -1;
}

//////////////////////////////////////////////////////////////////////

Position TreeDustGrid::centralPositionInCell(int m) const
{
    return Position(getNode(m)->extent().center());
}

//////////////////////////////////////////////////////////////////////

Position TreeDustGrid::randomPositionInCell(int m) const
{
    return _random->position(getNode(m)->extent());
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::path(DustGridPath* path) const
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

//////////////////////////////////////////////////////////////////////

vector<SimulationItem*> TreeDustGrid::interfaceCandidates(const std::type_info& interfaceTypeInfo)
{
    if (interfaceTypeInfo == typeid(DustGridDensityInterface) && !_dmib)
        return vector<SimulationItem*>();
    return BoxDustGrid::interfaceCandidates(interfaceTypeInfo);
}

//////////////////////////////////////////////////////////////////////

double TreeDustGrid::density(int h, int m) const
{
    TreeNode* node = getNode(m);
    return _dmib->massInBox(h, node->extent()) / node->volume();
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::write_xy(DustGridPlotFile* outfile) const
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

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::write_xz(DustGridPlotFile* outfile) const
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

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::write_yz(DustGridPlotFile* outfile) const
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

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::write_xyz(DustGridPlotFile* outfile) const
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

//////////////////////////////////////////////////////////////////////

TreeNode* TreeDustGrid::root() const
{
    return _tree[0];
}

//////////////////////////////////////////////////////////////////////

TreeNode* TreeDustGrid::getNode(int m) const
{
    return _tree[_idv[m]];
}

//////////////////////////////////////////////////////////////////////

int TreeDustGrid::cellNumber(const TreeNode* node) const
{
    return _cellnumberv[node->id()];
}

//////////////////////////////////////////////////////////////////////
