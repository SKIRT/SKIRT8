/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHAMRVACFILE_HPP
#define ADAPTIVEMESHAMRVACFILE_HPP

#include "AdaptiveMeshFile.hpp"
#include <fstream>

////////////////////////////////////////////////////////////////////

/** The AdaptiveMeshAmrvacFile class can read the relevant information on a cartesian
    three-dimensional (3D) Adaptive Mesh Refinement (AMR) grid from a file in the format produced
    by the MPI-AMRVAC code developed in Leuven; see for example Keppens et al 2012. The web pages
    http://homes.esat.kuleuven.be/~keppens/amrstructure.html and
    http://homes.esat.kuleuven.be/~keppens/fileformat.html provide an overview of the data format.
    Specific information was obtained from a python script kindly written for this purpose by Tom
    Hendrix.

    This class only supports cartesian grids. It converts 1D and 2D grids to 3D grids by assuming
    a thickness of 1 cell in the missing directions.

    For some reason the MPI-AMRVAC data file does not contain the size of the mesh at the coarsest
    level. This information must be provided separately as part of the ski file.
*/
class AdaptiveMeshAmrvacFile : public AdaptiveMeshFile
{
    ITEM_CONCRETE(AdaptiveMeshAmrvacFile, AdaptiveMeshFile, "an adaptive mesh data file in MPI-AMRVAC format")

    PROPERTY_INT(levelOneX, "the number of mesh cells at the coarsest level, in the X direction")
        ATTRIBUTE_MIN_VALUE(levelOneX, "1")
        ATTRIBUTE_MAX_VALUE(levelOneX, "10000")

    PROPERTY_INT(levelOneY, "the number of mesh cells at the coarsest level, in the Y direction")
        ATTRIBUTE_MIN_VALUE(levelOneY, "1")
        ATTRIBUTE_MAX_VALUE(levelOneY, "10000")

    PROPERTY_INT(levelOneZ, "the number of mesh cells at the coarsest level, in the Z direction")
        ATTRIBUTE_MIN_VALUE(levelOneZ, "1")
        ATTRIBUTE_MAX_VALUE(levelOneZ, "10000")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function opens the adaptive mesh data file, or throws a fatal error if the file can't
        be opened. It does not yet read any records. */
    void open() override;

    /** This function closes the adaptive mesh data file. */
    void close() override;

    /** This function reads the next record from the file, and holds its information ready for
        inspection through the other functions of this class. The function returns true if a record
        was successfully read, or false if the end of the file was reached or another error
        occurred. */
    bool read() override;

    /** This function returns true if the current record represents a nonleaf node, or false if
        the current record represents a leaf node. If there is no current record, the result is
        undefined. */
    bool isNonLeaf() const override;

    /** If the current record represents a nonleaf node, this function returns \f$N_x,N_y,N_z\f$,
        i.e. the number of child nodes carried by the node in each spatial direction. If the
        current record represents a leaf node or if there is no current record, the result is
        undefined. */
    void numChildNodes(int& nx, int& ny, int& nz) const override;

    /** If the current record represents a leaf node, this function returns the value \f$F_g\f$ of
        the field with given zero-based index \f$0\le g \le N_{fields}-1\f$. If the index is out of
        range, a fatal error is thrown. If the current record represents a nonleaf node or if there
        is no current record, the result is undefined. */
    double value(int g) const override;

    //========================= Data members =======================

private:
    // information about the mesh provided by the user
    int _nxlone[3];    // number of mesh cells at the coarsest level, in each direction

    // information about the mesh read from the input file
    int _nblocks;      // total number of blocks in the mesh
    int _ndims;        // dimension of the mesh (1D, 2D or 3D)
    int _nvars;        // number of variables in each cell
    int _nx[3];        // number of mesh cells in each block, in each direction
    int _ng[3];        // number of blocks at the coarsest level, in each direction
    int _nr[3];        // refinement factor for nested level, in each direction
    int _ncells;       // number of cells in a block
    int _blocksize;    // the size of a block in bytes
    vector<bool> _forest; // the forest representing the grid structure
    int _forestsize;   // the size of the forest vector

    // the input file and the current record
    std::ifstream _infile; // the input file
    int _state;            // indication of the current record (maintained by the "state machine" in read())
    vector<double> _block; // the values for all cells in the current block, if there is one
    int _forestindex;      // the index of the next forest item to be read
};

////////////////////////////////////////////////////////////////////

#endif
