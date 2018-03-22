/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILETREEDUSTGRID_HPP
#define FILETREEDUSTGRID_HPP

#include "BoxDustGrid.hpp"
#include "DustGridDensityInterface.hpp"
class DustMassInBoxInterface;
class TreeNode;
class Random;

//////////////////////////////////////////////////////////////////////

/** FileTreeDustGrid represents hierarchical tree dust grids for which the structure is imported
    from a file previously written by an instance of the TreeDustGrid class. Depending on the
    incoming data, the tree can be an octtree (8 children per node) or a kd-tree (2 children per
    node). See TreeDustGrid for more information.

    Note: the current implementation involves a substantial amount of code duplication from
          TreeDustGrid because of complications with the inheritance structure.
*/
class FileTreeDustGrid : public DustGrid, public Box, public DustGridDensityInterface
{
    /** The enumeration type indicating the search method to be used for finding the subsequent
        node while traversing the tree grid. The TopDown method (the default) always starts at the
        root node and recursively finds the child node containing the new position. The Neighbor
        method constructs a neighbor list for each node (at each of the six walls) during setup,
        and then uses this list to locate the neighboring node containing the new position. The
        Bookkeeping method relies on the order in which the occtree nodes are created and stored to
        derive the appropriate neighbor solely through the respective node indices. */
    ENUM_DEF(SearchMethod, TopDown, Neighbor, Bookkeeping)
    ENUM_VAL(SearchMethod, TopDown, "top-down (start at root and recursively find appropriate child node)")
    ENUM_VAL(SearchMethod, Neighbor, "neighbor (construct and use neighbor list for each node wall) ")
    ENUM_VAL(SearchMethod, Bookkeeping, "bookkeeping (derive appropriate neighbor through node indices)")
    ENUM_END()

    ITEM_CONCRETE(FileTreeDustGrid, DustGrid, "an octtree or k-d tree dust grid imported from file")

    PROPERTY_STRING(filename, "the name of the input file representing the tree grid")

    PROPERTY_ENUM(searchMethod, SearchMethod, "the search method used for traversing the tree grid")
        ATTRIBUTE_DEFAULT_VALUE(searchMethod, "Neighbor")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor deletes all nodes from the tree vector created during setup. */
    ~FileTreeDustGrid();

protected:
    /** This function . */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is 3 for all subclasses of this class. */
    int dimension() const override;

    /** This function returns the bounding box that encloses the dust grid. */
    Box boundingBox() const override;

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
        grid. To determine the cell numbers in this algorithm, the function uses the method
        configured with setSearchMethod(). */
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
    vector<TreeNode*> _tree;
    vector<int> _cellnumberv;
    vector<int> _idv;
    int _Nnodes{0};
    int _highestWriteLevel{0};
    double _eps{0.};
    Random* _random{nullptr};
    DustMassInBoxInterface* _dmib{nullptr};
};

//////////////////////////////////////////////////////////////////////

#endif
