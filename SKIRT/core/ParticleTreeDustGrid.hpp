/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARTICLETREEDUSTGRID_HPP
#define PARTICLETREEDUSTGRID_HPP

#include "BoxDustGrid.hpp"
#include "DustGridDensityInterface.hpp"
class DustMassInBoxInterface;
class Random;
class TreeNode;

//////////////////////////////////////////////////////////////////////

/** ParticleTreeDustGrid is a concrete subclass of the BoxDustGrid class, representing
    three-dimensional dust grids with cuboidal cells organized in a tree. This class is similar to
    the TreeDustGrid class, except that the structure of the tree is determined by a set of
    particle locations retrieved from the dust distribution through the DustParticleInterface
    interface. In a first step the tree is subdivided in such a way that each leaf cell contains at
    most one particle. Subsequently each leaf node can be further subdivided a fixed number of
    times. Depending on the type of TreeNode, the tree can become an octtree (8 children per node)
    or a kd-tree (2 children per node). Other node types could be implemented, as long as they are
    cuboids lined up with the axes. The current implementation always uses regular (geocentric)
    subdivision and top-down search for locating cells. */
class ParticleTreeDustGrid : public BoxDustGrid, public DustGridDensityInterface
{
    /** The enumeration type indicating the type of tree to be constructed: an octtree (8 children
        per node) or a kd-tree (2 children per node). */
    ENUM_DEF(TreeType, OctTree, BinTree)
    ENUM_VAL(TreeType, OctTree, "an octtree (8 children per node)")
    ENUM_VAL(TreeType, BinTree, "a kd-tree (2 children per node)")
    ENUM_END()

    ITEM_CONCRETE(ParticleTreeDustGrid, BoxDustGrid, "a tree dust grid derived from a set of particles")

    PROPERTY_ENUM(treeType, TreeType, "the type of tree")
        ATTRIBUTE_DEFAULT_VALUE(treeType, "OctTree")

    PROPERTY_INT(numExtraLevels, "the number of additional subdivision levels")
        ATTRIBUTE_MIN_VALUE(numExtraLevels, "0")
        ATTRIBUTE_MAX_VALUE(numExtraLevels, "30")
        ATTRIBUTE_DEFAULT_VALUE(numExtraLevels, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor deletes all nodes from the tree vector created during setup. */
    ~ParticleTreeDustGrid();

protected:
    /** This function constructs the tree. The particle locations are retrieved from the dust
        distribution through the DustParticleInterface interface, and the tree nodes are subdivided
        (using regular subdivision) until each leaf cell contains at most one particle. If
        requested, each leaf node is further subdivided by a fixed number of levels. When this task
        is accomplished, the function creates a vector that contains the node IDs of all leaves.
        This is the actual dust cell vector (only the leaf nodes are the actual dust cells). The
        function also creates a vector with the cell numbers of all the nodes, i.e. the rank
        \f$m\f$ of the node in the ID vector if the node is a leaf, and the number -1 if the node
        is not a leaf (and hence not a dust cell). Finally, the function logs some details on the
        number of nodes and the number of cells, and if writeFlag() returns true, it writes the
        distribution of the grid cells to a file. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the volume of the dust cell with cell number \f$m\f$. For a tree dust
        grid, it determines the node ID corresponding to the cell number \f$m\f$, and then simply
        calculates the volume of the corresponding cuboidal cell using \f$V = \Delta x\, \Delta y\,
        \Delta z\f$. */
    double volume(int m) const override;

    /** This function returns the number of cells in the dust grid. */
    int numCells() const override;

    /** This function returns the number of the dust cell that contains the position
        \f${\bf{r}}\f$. For a tree dust grid, the search algorithm starts at the root node and
        selects the child node that contains the position. This procedure is repeated until the
        node is childless, i.e. until it is a leaf node that corresponds to an actual dust cell. */
    int whichCell(Position bfr) const override;

    /** This function returns the central location of the dust cell with cell number \f$m\f$. For a
        tree dust grid, it determines the node ID corresponding to the cell number \f$m\f$, and
        then calculates the central position in that cell through \f[ \begin{split} x &=
        x_{\text{min}} + \frac12\, \Delta x \\ y &= y_{\text{min}} + \frac12\, \Delta y \\ z &=
        z_{\text{min}} + \frac12\, \Delta z \end{split} \f] */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the dust cell with cell number \f$m\f$. For a
        tree dust grid, it determines the node ID corresponding to the cell number \f$m\f$, and
        then calculates a random position in that cell through \f[ \begin{split} x &=
        x_{\text{min}} + {\cal{X}}_1\, \Delta x \\ y &= y_{\text{min}} + {\cal{X}}_2\, \Delta y \\
        z &= z_{\text{min}} + {\cal{X}}_3\, \Delta z \end{split} \f] with \f${\cal{X}}_1\f$,
        \f${\cal{X}}_2\f$ and \f${\cal{X}}_3\f$ three uniform deviates. */
    Position randomPositionInCell(int m) const override;

    /** This function calculates a path through the grid. The DustGridPath object passed as an
        argument specifies the starting position \f${\bf{r}}\f$ and the direction \f${\bf{k}}\f$
        for the path. The data on the calculated path are added back into the same object. For a
        tree dust grid, the function uses a rather straighforward algorithm. It determines the dust
        cell that contains the starting position, and calculates the first wall of the cell that
        will be crossed. The pathlength \f$\Delta s\f$ is determined and the current position is
        moved to a new position along this path, a tiny fraction further than \f$\Delta s\f$, \f[
        \begin{split} x_{\text{new}} &= x_{\text{current}} + (\Delta s + \epsilon)\,k_x \\
        y_{\text{new}} &= y_{\text{current}} + (\Delta s + \epsilon)\,k_y \\ z_{\text{new}} &=
        z_{\text{current}} + (\Delta s + \epsilon)\,k_z \end{split} \f] where \f[ \epsilon =
        10^{-12} \sqrt{x_{\text{max}}^2 + y_{\text{max}}^2 + z_{\text{max}}^2} \f] By adding this
        small extra bit, we ensure that the new position is now within the next cell, and we can
        repeat this exercise. This loop is terminated when the next position is outside the dust
        grid. */
    void path(DustGridPath* path) const override;

    /** This function is used by the interface() template function in the SimulationItem class. It
        returns a list of simulation items that should be considered in the search for an item that
        implements the requested interface. The implementation in this class returns the default
        list (i.e. the receiving object) except in the following case. If the requested interface
        is DustGridDensityInterface (which is implemented by this class) and the dust distribution
        for this simulation does \em not offer the DustMassInBoxInterface interface, the returned
        list is empty (because it is then impossible to provide DustGridDensityInterface). */
    vector<SimulationItem*> interfaceCandidates(const std::type_info& interfaceTypeInfo) override;

    /** This function implements the DustGridDensityInterface interface. It returns the density for
        the dust component \em h in the dust grid cell with index \em m. The implementation relies
        on the DustMassInBoxInterface interface in the dust distribution for this simulation. */
    double density(int h, int m) const override;

protected:
    /** This function writes the intersection of the dust grid with the xy plane to the specified
        DustGridPlotFile object. */
    void write_xy(DustGridPlotFile* outfile) const override;

    /** This function writes the intersection of the dust grid with the xz plane to the specified
        DustGridPlotFile object. */
    void write_xz(DustGridPlotFile* outfile) const override;

    /** This function writes the intersection of the dust grid with the yz plane to the specified
        DustGridPlotFile object. */
    void write_yz(DustGridPlotFile* outfile) const override;

    /** This function writes 3D information for the cells up to a certain level in the dust grid
        structure to the specified DustGridPlotFile object. The output is restricted to limit the
        number of cells being written (to keep the line density in the output plot within reason).
        */
    void write_xyz(DustGridPlotFile* outfile) const override;

private:
    /** This function returns a pointer to the root node of the tree. */
    TreeNode* root() const;

    /** This function returns a pointer to the node corresponding to cell number \f$m\f$. It just
        reads the node ID number of the \f$m\f$'th leave cell and returns the corresponding pointer
        of the tree vector. */
    TreeNode* getNode(int m) const;

    /** This function returns the cell number \f$m\f$ of a node in the tree. It just reads the ID
        of the node and determines the corresponding cell number from the internal cell number
        vector. */
    int cellNumber(const TreeNode* node) const;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Random* _random{nullptr};
    DustMassInBoxInterface* _dmib{nullptr};
    double _eps{0.};
    std::vector<TreeNode*> _tree;
    std::vector<int> _cellnumberv;
    std::vector<int> _idv;
    int _highestWriteLevel{0};
};

//////////////////////////////////////////////////////////////////////

#endif
