/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustGrid.hpp"
#include "DustGridPlotFile.hpp"
#include "FatalError.hpp"
#include "MonteCarloSimulation.hpp"

//////////////////////////////////////////////////////////////////////

void DustGrid::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // verify that the grid has at least the dimension of the simulation's star and dust geometry
    int dimGeometry = find<MonteCarloSimulation>()->dimension();
    int dimGrid = dimension();
    if (dimGeometry > dimGrid)
        throw FATALERROR("The grid dimension " + std::to_string(dimGrid) +
                         " is lower than the geometry dimension " + std::to_string(dimGeometry));
}

//////////////////////////////////////////////////////////////////////

void DustGrid::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    if (_writeGrid) performWriteGrid();
}

//////////////////////////////////////////////////////////////////////

void DustGrid::performWriteGrid() const
{
    int dimension = find<MonteCarloSimulation>()->dimension();

    // For the xy plane (always)
    {
        DustGridPlotFile outfile(this, "ds_gridxy");
        write_xy(&outfile);
    }

    // For the xz plane (only if dimension is at least 2)
    if (dimension >= 2)
    {
        DustGridPlotFile outfile(this, "ds_gridxz");
        write_xz(&outfile);
    }

    // For the yz plane (only if dimension is 3)
    if (dimension == 3)
    {
        DustGridPlotFile outfile(this, "ds_gridyz");
        write_yz(&outfile);
    }

    // Full 3D coordinates (only if dimension is 3)
    if (dimension == 3)
    {
        DustGridPlotFile outfile(this, "ds_gridxyz");
        write_xyz(&outfile);
    }
}

//////////////////////////////////////////////////////////////////////

double DustGrid::weight(int m) const
{
    return (m==-1) ? 0.0 : 1.0;
}

//////////////////////////////////////////////////////////////////////

void DustGrid::write_xy(DustGridPlotFile* /*outfile*/) const
{
}

//////////////////////////////////////////////////////////////////////

void DustGrid::write_xz(DustGridPlotFile* /*outfile*/) const
{
}

//////////////////////////////////////////////////////////////////////

void DustGrid::write_yz(DustGridPlotFile* /*outfile*/) const
{
}

//////////////////////////////////////////////////////////////////////

void DustGrid::write_xyz(DustGridPlotFile* /*outfile*/) const
{
}

//////////////////////////////////////////////////////////////////////
