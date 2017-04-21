/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VORONOIDUSTGRID_HPP
#define VORONOIDUSTGRID_HPP

#include "BoxDustGrid.hpp"
#include "VoronoiMeshFile.hpp"
class Random;
class VoronoiMesh;

//////////////////////////////////////////////////////////////////////

/** VoronoiDustGrid is a concrete subclass of the BoxDustGrid class. It is used for representing a
    three-dimensional dust grid based on a Voronoi tesselation of the cuboid containing
    substantially all of the dust. See the VoronoiMesh class for more information on Voronoi
    tesselations. The class offers several options for determining the locations of the particles
    generating the Voronoi tesselation. A specified number of particles can be distributed randomly
    over the domain, either uniformly or with the same overall density distribution as the dust.
    Alternatively, the locations can be copied from the particles in an SPH dust distribution. This
    class uses the Voro++ code written by Chris H. Rycroft (LBL / UC Berkeley) to generate output
    files for plotting the Voronoi grid. */
class VoronoiDustGrid : public BoxDustGrid
{
    /** The enumeration type indicating the probability distribution used for generating the random
        particles. */
    ENUM_DEF(Distribution, Uniform, CentralPeak, DustDensity, DustTesselation, SPHParticles, File)
    ENUM_VAL(Distribution, Uniform, "uniform")
    ENUM_VAL(Distribution, CentralPeak, "with a steep central peak")
    ENUM_VAL(Distribution, DustDensity, "same as dust density")
    ENUM_VAL(Distribution, DustTesselation, "copy Voronoi tesselation in dust distribution")
    ENUM_VAL(Distribution, SPHParticles, "exact locations of SPH particles in dust distribution")
    ENUM_VAL(Distribution, File, "particle locations in a specified data file")
    ENUM_END()

    ITEM_CONCRETE(VoronoiDustGrid, BoxDustGrid, "a Voronoi dust grid")

    PROPERTY_INT(numParticles, "the number of random particles (or cells in the grid)")
        ATTRIBUTE_MIN_VALUE(numParticles, "5")
        ATTRIBUTE_DEFAULT_VALUE(numParticles, "500")

    PROPERTY_ENUM(distribution, Distribution, "the probablity distribution for the particles")
        ATTRIBUTE_DEFAULT_VALUE(distribution, "DustDensity")
        ATTRIBUTE_TRUE_IF(distribution, "File")

    PROPERTY_ITEM(voronoiMeshFile, VoronoiMeshFile, "the Voronoi mesh data file")
        ATTRIBUTE_RELEVANT_IF(voronoiMeshFile, "distribution")
        ATTRIBUTE_DEFAULT_VALUE(voronoiMeshFile, "VoronoiMeshAsciiFile")
        ATTRIBUTE_OPTIONAL(voronoiMeshFile)

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor cleans up all memory allocated in this class. */
    ~VoronoiDustGrid();

protected:
    /** This function verifies that the attributes have been appropriately set, selects the
        requested particles for generating the Voronoi tesselation, and finally constructs the
        Voronoi tesselation through an instance of the VoronoiMesh class. If requested, it also
        outputs files that can be used for plotting the grid. To generate the random particles
        according to the simulation's dust density distribution, the setup function partitions the
        domain using a three-dimensional cuboidal cell structure, called the foam. The distribution
        of the foam cells is determined automatically from the dust density distribution. The foam
        allows to efficiently generate random points drawn from this probability distribution. Once
        the particles have been generated, the foam is discarded. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function outputs the plot files; it is overriden here with a special version because
        the regular mechanism does not apply. The function reconstructs the Voronoi tesselation in
        order to produce the coordinates of the Voronoi cell vertices. */
    void performWriteGrid() const override;

    /** This function returns the volume of the dust cell with cell number \f$m\f$. */
    double volume(int m) const override;

    /** This function returns the number of cells in the dust grid. */
    int numCells() const override;

    /** This function returns the number of the dust cell that contains the position
        \f${\bf{r}}\f$. See the VoronoiMesh class for more information. */
    int whichCell(Position bfr) const override;

    /** This function returns the central location of the dust cell with cell number \f$m\f$. In
        this class the function returns the centroid of the Voronoi cell. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the dust cell with cell number \f$m\f$. */
    Position randomPositionInCell(int m) const override;

    /** This function calculates a path through the grid. The DustGridPath object passed as an
        argument specifies the starting position \f${\bf{r}}\f$ and the direction \f${\bf{k}}\f$
        for the path. The data on the calculated path are added back into the same object. See the
        VoronoiMesh class for more information. */
    void path(DustGridPath* path) const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Random* _random{nullptr};
    VoronoiMesh* _mesh{nullptr};
    bool _meshOwned{true};
};

//////////////////////////////////////////////////////////////////////

#endif
