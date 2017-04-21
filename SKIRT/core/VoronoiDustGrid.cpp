/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoronoiDustGrid.hpp"
#include "DustDistribution.hpp"
#include "DustGridPlotFile.hpp"
#include "DustParticleInterface.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "VoronoiMesh.hpp"
#include "VoronoiMeshFile.hpp"
#include "VoronoiMeshInterface.hpp"
#include "container.hh"

//////////////////////////////////////////////////////////////////////

VoronoiDustGrid::~VoronoiDustGrid()
{
    if (_meshOwned) delete _mesh;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setupSelfBefore()
{
    BoxDustGrid::setupSelfBefore();

    // Cache the random number generator
    _random = find<Random>();

    Log* log = find<Log>();

    // Determine an appropriate set of particles and construct the Voronoi mesh
    switch (_distribution)
    {
    case Distribution::Uniform:
        {
            vector<Vec> rv(_numParticles);
            for (int m=0; m<_numParticles; m++)
            {
                rv[m] = _random->position(extent());
            }
            log->info("Computing Voronoi tesselation for " + std::to_string(_numParticles)
                      + " uniformly distributed random particles...");
            _mesh = new VoronoiMesh(rv, extent(), log);
            break;
        }
    case Distribution::CentralPeak:
        {
            const int a = 1000;                     // steepness of the peak; the central 1/a portion is NOT covered
            const double rscale = extent().rmax().norm();
            vector<Vec> rv(_numParticles);
            for (int m=1; m<_numParticles; m++)     // skip first particle so that it remains (0,0,0)
            {
                while (true)
                {
                    double r = rscale * pow(1./a, _random->uniform());   // random distribution according to 1/x
                    Direction k = _random->direction();
                    Position p = Position(r,k);
                    if (extent().contains(p))       // discard any points outside of the domain
                    {
                        rv[m] = p;
                        break;
                    }
                }
            }
            log->info("Computing Voronoi tesselation for " + std::to_string(_numParticles)
                      + " random particles distributed in a central peak...");
            _mesh = new VoronoiMesh(rv, extent(), log);
            break;
        }
    case Distribution::DustDensity:
        {
            DustDistribution* dd = find<DustDistribution>();
            vector<Vec> rv(_numParticles);
            for (int m=0; m<_numParticles; m++)
            {
                while (true)
                {
                    Position p = dd->generatePosition();
                    if (extent().contains(p))       // discard any points outside of the domain
                    {
                        rv[m] = p;
                        break;
                    }
                }
            }
            log->info("Computing Voronoi tesselation for " + std::to_string(_numParticles)
                      + " random particles distributed according to dust density...");
            _mesh = new VoronoiMesh(rv, extent(), log);
            break;
        }
    case Distribution::DustTesselation:
        {
            VoronoiMeshInterface* vmi = find<DustDistribution>()->interface<VoronoiMeshInterface>();
            if (!vmi) throw FATALERROR("Can't retrieve Voronoi mesh from this dust distribution");
            _mesh = vmi->mesh();
            _meshOwned = false;
            log->info("Using Voronoi tesselation from dust distribution with " + std::to_string(_mesh->numCells())
                      + " particles...");
            break;
        }
    case Distribution::SPHParticles:
        {
            DustParticleInterface* dpi = find<DustDistribution>()->interface<DustParticleInterface>();
            if (!dpi) throw FATALERROR("Can't retrieve particle locations from this dust distribution");
            log->info("Computing Voronoi tesselation for " + std::to_string(dpi->numParticles())
                      + " dust distribution particles...");
            _mesh = new VoronoiMesh(dpi, extent(), log);
            break;
        }
    case Distribution::File:
        {
            log->info("Computing Voronoi tesselation for particles loaded from file "
                      + _voronoiMeshFile->filename() + "...");
            _mesh = new VoronoiMesh(_voronoiMeshFile, vector<int>(), extent(), log);
            break;
        }
    }

    int Ncells = _mesh->numCells();

    // Log statistics on the cell neighbors
    double avgNeighbors;
    int minNeighbors, maxNeighbors;
    _mesh->neighborStatistics(avgNeighbors, minNeighbors, maxNeighbors);
    log->info("Computed Voronoi tesselation with " + std::to_string(Ncells) + " cells:");
    log->info("  Average number of neighbors per cell: " + StringUtils::toString(avgNeighbors,'f',1));
    log->info("  Minimum number of neighbors per cell: " + std::to_string(minNeighbors));
    log->info("  Maximum number of neighbors per cell: " + std::to_string(maxNeighbors));

    // Log statistics on the block lists
    int nblocks = _mesh->numBlocks();
    double avgRefsPerBlock;
    int minRefsPerBlock, maxRefsPerBlock;
    _mesh->blockStatistics(avgRefsPerBlock, minRefsPerBlock, maxRefsPerBlock);
    log->info("Created grid to accelerate which-cell operations:");
    log->info("  Number of cells                  : " + std::to_string(Ncells));
    log->info("  Number of blocks                 : " + std::to_string(nblocks*nblocks*nblocks) +
              " (" + std::to_string(nblocks) + " in each dimension)");
    log->info("  Average number of cells per block: " + StringUtils::toString(avgRefsPerBlock,'f',1));
    log->info("  Minimum number of cells per block: " + std::to_string(minRefsPerBlock));
    log->info("  Maximum number of cells per block: " + std::to_string(maxRefsPerBlock));

    // Log statistics on the search trees
    double avgRefsPerTree;
    int nTrees, minRefsPerTree, maxRefsPerTree;
    _mesh->treeStatistics(nTrees, avgRefsPerTree, minRefsPerTree, maxRefsPerTree);
    log->info("Created search trees to accelerate which-cell operations:");
    log->info("  Number of trees                  : " + std::to_string(nTrees) +
              " (" + StringUtils::toString(100.*nTrees/(nblocks*nblocks*nblocks),'f',1) + "% of blocks)");
    log->info("  Average number of cells per tree : " + StringUtils::toString(avgRefsPerTree,'f',1));
    log->info("  Minimum number of cells per tree : " + std::to_string(minRefsPerTree));
    log->info("  Maximum number of cells per tree : " + std::to_string(maxRefsPerTree));
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::performWriteGrid() const
{
    // create the plot files
    DustGridPlotFile plotxy(this, "ds_gridxy");
    DustGridPlotFile plotxz(this, "ds_gridxz");
    DustGridPlotFile plotyz(this, "ds_gridyz");
    DustGridPlotFile plotxyz(this, "ds_gridxyz");

    // load all particles in a Voro container
    int Ncells = _mesh->numCells();
    int nb = max(3, min(1000, static_cast<int>(pow(Ncells/5.,1./3.)) ));
    voro::container con(xmin(), xmax(), ymin(), ymax(), zmin(), zmax(), nb, nb, nb, false,false,false, 8);
    for (int m=0; m<Ncells; m++)
    {
        Vec r = _mesh->particlePosition(m);
        con.put(m, r.x(),r.y(),r.z());
    }

    // loop over all Voro cells
    voro::c_loop_all loop(con);
    if (loop.start()) do
    {
        // Compute the cell
        voro::voronoicell fullcell;
        con.compute_cell(fullcell, loop);

        // Get the edges of the cell
        double x,y,z;
        loop.pos(x,y,z);
        vector<double> coords;
        fullcell.vertices(x,y,z, coords);
        vector<int> indices;
        fullcell.face_vertices(indices);

        // Write the edges of the cell to the plot files
        Box bounds = _mesh->extent(loop.pid());
        if (bounds.zmin()<=0 && bounds.zmax()>=0) plotxy.writePolyhedron(coords, indices);
        if (bounds.ymin()<=0 && bounds.ymax()>=0) plotxz.writePolyhedron(coords, indices);
        if (bounds.xmin()<=0 && bounds.xmax()>=0) plotyz.writePolyhedron(coords, indices);
        if (loop.pid() <= 1000) plotxyz.writePolyhedron(coords, indices);
    }
    while (loop.inc());
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustGrid::volume(int m) const
{
    return _mesh->volume(m);
}

//////////////////////////////////////////////////////////////////////

int VoronoiDustGrid::numCells() const
{
    return _mesh->numCells();
}

//////////////////////////////////////////////////////////////////////

int VoronoiDustGrid::whichCell(Position bfr) const
{
    return _mesh->cellIndex(bfr);
}

//////////////////////////////////////////////////////////////////////

Position VoronoiDustGrid::centralPositionInCell(int m) const
{
    return _mesh->centralPosition(m);
}

//////////////////////////////////////////////////////////////////////

Position VoronoiDustGrid::randomPositionInCell(int m) const
{
    return _mesh->randomPosition(_random, m);
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::path(DustGridPath* path) const
{
    _mesh->path(path);
}

//////////////////////////////////////////////////////////////////////
