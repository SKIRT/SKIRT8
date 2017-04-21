/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TREEDUSTGRID_HPP
#define TREEDUSTGRID_HPP

#include "BoxDustGrid.hpp"
#include "DustGridDensityInterface.hpp"
class DustDistribution;
class DustMassInBoxInterface;
class TreeNode;
class Parallel;
class ProcessAssigner;
class Random;

//////////////////////////////////////////////////////////////////////

/** TreeDustGrid is an abstract subclass of the BoxDustGrid class, and represents three-dimensional
    dust grids with cuboidal cells organized in a tree. The tree's root node encloses the complete
    spatial domain, and nodes on subsequent levels recursively divide space into ever finer nodes.
    The depth of the tree can vary from place to place. The leaf cells (those that are not further
    subdivided) are the actual dust cells. The type of TreeNode used by the TreeDustGrid
    is decided in each subclass through a factory method. Depending on the type of TreeNode, the
    tree can become an octtree (8 children per node) or a kd-tree (2 children per node). Other node
    types could be implemented, as long as they are cuboids lined up with the axes. */
class TreeDustGrid : public BoxDustGrid, public DustGridDensityInterface
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

    ITEM_ABSTRACT(TreeDustGrid, BoxDustGrid, "a tree dust grid")

    PROPERTY_INT(minLevel, "the minimum level of grid refinement (typically 2 to 3)")
        ATTRIBUTE_MIN_VALUE(minLevel, "0")
        ATTRIBUTE_MAX_VALUE(minLevel, "50")
        ATTRIBUTE_DEFAULT_VALUE(minLevel, "2")

    PROPERTY_INT(maxLevel, "the maximum level of grid refinement (typically 6 to 10)")
        ATTRIBUTE_MIN_VALUE(maxLevel, "2")
        ATTRIBUTE_MAX_VALUE(maxLevel, "50")
        ATTRIBUTE_DEFAULT_VALUE(maxLevel, "6")

    PROPERTY_ENUM(searchMethod, SearchMethod, "the search method used for traversing the tree grid")
        ATTRIBUTE_DEFAULT_VALUE(searchMethod, "Neighbor")

    PROPERTY_INT(numSamples, "the number of random density samples for determining cell subdivision")
        ATTRIBUTE_MIN_VALUE(numSamples, "10")
        ATTRIBUTE_MAX_VALUE(numSamples, "1000")
        ATTRIBUTE_DEFAULT_VALUE(numSamples, "100")

    PROPERTY_DOUBLE(maxOpticalDepth, "the maximum mean optical depth for each dust cell")
        ATTRIBUTE_MIN_VALUE(maxOpticalDepth, "[0")
        ATTRIBUTE_MAX_VALUE(maxOpticalDepth, "100]")
        ATTRIBUTE_DEFAULT_VALUE(maxOpticalDepth, "0")

    PROPERTY_DOUBLE(maxMassFraction, "the maximum fraction of dust mass contained in each dust cell")
        ATTRIBUTE_MIN_VALUE(maxMassFraction, "[0")
        ATTRIBUTE_MAX_VALUE(maxMassFraction, "1e-2]")
        ATTRIBUTE_DEFAULT_VALUE(maxMassFraction, "1e-6")

    PROPERTY_DOUBLE(maxDensityDispersion,
                    "the maximum density dispersion in each dust cell, as fraction of the reference density")
        ATTRIBUTE_MIN_VALUE(maxDensityDispersion, "[0")
        ATTRIBUTE_MAX_VALUE(maxDensityDispersion, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxDensityDispersion, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor deletes all nodes from the tree vector created during setup. */
    ~TreeDustGrid();

protected:
    /** This function verifies that all attribute values have been appropriately set and actually
        constructs the tree. The first step is to create the root node (through the factory method
        createRoot() to be implemented in each subclass), and store it in the tree vector, which is
        just a list of pointers to nodes). The second phase is to recursively subdivide the root
        node and add the children at the end of the tree vector, until all nodes satisfy the
        criteria for no further subdivision. When this task is accomplished, the function creates a
        vector that contains the node IDs of all leaves. This is the actual dust cell vector (only
        the leaf nodes are the actual dust cells). The function also creates a vector with the cell
        numbers of all the nodes, i.e. the rank \f$m\f$ of the node in the ID vector if the node is
        a leaf, and the number -1 if the node is not a leaf (and hence not a dust cell). Finally,
        the function logs some details on the number of nodes and the number of cells, and if
        writeFlag() returns true, it writes the distribution of the grid cells to a file. */
    void setupSelfBefore() override;

private:
    /** This function, only to be called during the construction phase, investigates whether a node
        should be further subdivided and also takes care of the actual subdivision. There are
        several criteria for subdivision. The simplest criterion is the level of subdivision of the
        node: if it is less then a minimum level, the node is always subdivided, if it higher then
        a maximum level, there is no subdivision (these levels are input parameters). In the
        general case, whether or not we subdivide depends on the following criterion: if the ratio
        of the dust mass in the cell and the total dust mass, corresponding to the dust density
        distribution \f$\rho({\bf{r}})\f$, is larger than a preset threshold (an input parameter as
        well), the node is subdivided. The dust mass in the cell is calculated by generating
        \f$N_{\text{random}}\f$ random positions \f${\bf{r}}_n\f$ in the node, so that the total
        mass using the density in these points is given by \f[ M \approx \frac{\Delta x\, \Delta
        y\, \Delta z}{N_{\text{random}}} \sum_{n=0}^{N_{\text{random}}} \rho({\bf{r}}_n), \f] If
        the mass criterion is not satisfied, the node is not subdivided. If the conditions for
        subdivision are met, the first task is to calculate the division point. There are two
        options: geocentric division and barycentric division (an input flag again). In the former
        case, the division point of the cell is just the geometric centre, i.e. \f[ {\bf{r}}_c =
        (x_c,y_c,z_c) = \left( \frac{x_{\text{min}}+x_{\text{max}}}{2},
        \frac{y_{\text{min}}+y_{\text{max}}}{2}, \frac{z_{\text{min}}+z_{\text{max}}}{2} \right).
        \f] In the latter case the division point is the centre of mass, which we estimate using
        the \f$N_{\text{random}}\f$ points generated before, \f[ {\bf{r}}_c = \frac{ \sum_n
        \rho({\bf{r}}_n)\, {\bf{r}}_n}{ \sum_n \rho({\bf{r}}_n) }. \f] The last task is to actually
        create the eight child nodes of the node and add them to the tree. */
    void subdivide(TreeNode* node);

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

protected:
    /** This pure virtual function, to be implemented in each subclass, creates a root node of the
        appropriate type, using a node identifier of zero and the specified spatial extent, and
        returns a pointer to it. The caller takes ownership of the newly created object. */
    virtual TreeNode* createRoot(const Box& extent) = 0;

    /** This pure virtual function, to be implemented in each subclass, returns true if this type
        of TreeDustGrid instance allows the use of the DustMassInBoxInterface for determining
        subdivisions, and false if not. */
    virtual bool canUseDmibForSubdivide() = 0;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Random* _random{nullptr};
    Parallel* _parallel{nullptr};
    DustDistribution* _dd{nullptr};
    DustMassInBoxInterface* _dmib{nullptr};
    double _totalmass{0.};
    double _eps{0.};
    int _Nnodes{0};
    vector<TreeNode*> _tree;
    vector<int> _cellnumberv;
    vector<int> _idv;
    int _highestWriteLevel{0};
    bool _useDmibForSubdivide{false};
};

//////////////////////////////////////////////////////////////////////

#endif
