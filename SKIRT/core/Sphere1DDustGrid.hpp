/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHERE1DDUSTGRID_HPP
#define SPHERE1DDUSTGRID_HPP

#include "SphereDustGrid.hpp"
#include "Array.hpp"
#include "Mesh.hpp"
class Random;

////////////////////////////////////////////////////////////////////

/** The Sphere1DDustGrid class is a subclass of the SphereDustGrid class, and represents
    one-dimensional, spherically symmetric dust grids. Each cell in such a grid is a spherical
    shell. Internally, a spherical dust grid is specified through a set of \f$N_r+1\f$ radial grid
    points \f$r_i\f$ (with \f$i=0,\ldots,N_r\f$). */
class Sphere1DDustGrid : public SphereDustGrid
{
    ITEM_CONCRETE(Sphere1DDustGrid, SphereDustGrid, "a spherically symmetric dust grid")
        ATTRIBUTE_ALLOWED_IF(Sphere1DDustGrid,
                             "(!AxGeometry)&(!GenGeometry)&(!NetzerAccretionDiskGeometry)&(!SolarPatchGeometry)"
                             "&(!StellarSurfaceGeometry)&(!OffsetGeometryDecorator)&(!AdaptiveMeshStellarComp)"
                             "&(!SPHStellarComp)&(!VoronoiStellarComp)&(!AdaptiveMeshDustDistribution)"
                             "&(!SPHDustDistribution)&(!VoronoiDustDistribution)")

    PROPERTY_ITEM(meshRadial, Mesh, "the bin distribution in the radial direction")
        ATTRIBUTE_DEFAULT_VALUE(meshRadial, "LinMesh")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function sets up a number of data members that depend on the Mesh object configured
        for this grid. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is 1 for this class. */
    int dimension() const override;

    /** This function returns the number of cells in the dust grid. */
    int numCells() const override;

    /** This function returns the volume of the dust cell with cell number \f$m\f$. For a spherical
        dust grid, cell number \f$m\f$ corresponds to the radial bin with lower border index
        \f$i=m\f$, and the volume is easily calculated as \f[V = \frac{4\pi}{3}\,
        (r_{i+1}^3-r_i^3),\f] with \f$r_i\f$ and \f$r_{i+1}\f$ the inner and outer radius of the
        shell respectively. */
    double volume(int m) const override;

    /** This function returns the number of the dust cell that contains the position
        \f${\bf{r}}\f$. It just determines the radial bin index and returns that number. */
    int whichCell(Position bfr) const override;

    /** This function returns the central location from the dust cell with cell number \f$m\f$. For
        a spherical dust grid, cell number \f$m\f$ corresponds to the radial bin with lower border
        index \f$i=m\f$, and the central radius is determined using \f[ r = \frac{r_i +
        r_{i+1}}{2}. \f] The returned position is arbitrarily located on the x-axis. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the dust cell with cell number \f$m\f$. For a
        spherical dust grid, cell number \f$m\f$ corresponds to the radial bin with lower border
        index \f$i=m\f$, and a random radius is determined using \f[ r = r_i +
        {\cal{X}}\,(r_{i+1}-r_i) \f] with \f${\cal{X}}\f$ a random deviate. This random radius is
        combined with a random position on the unit sphere to generate a random position from the
        cell. */
    Position randomPositionInCell(int m) const override;

    /** This function calculates a path through the grid. The DustGridPath object passed as an
        argument specifies the starting position \f${\bf{r}}\f$ and the direction \f${\bf{k}}\f$
        for the path. The data on the calculated path are added back into the same object. */
    void path(DustGridPath* path) const override;

private:
    /** This function writes the intersection of the dust grid with the xy plane to the specified
        DustGridPlotFile object. */
    void write_xy(DustGridPlotFile* outfile) const override;

    //======================== Data Members ========================

private:
    Random* _random{nullptr};
    int _Nr{0};
    Array _rv;
};

////////////////////////////////////////////////////////////////////

#endif
