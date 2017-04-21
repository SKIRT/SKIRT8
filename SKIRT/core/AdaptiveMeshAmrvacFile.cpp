/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshAmrvacFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    const int INT_SIZE = sizeof(int32_t);
    const int DOUBLE_SIZE = sizeof(double);

    // read a 32-bit integer from the data stream
    // !!! this function assumes that the incoming data has native byte order !!!
    int32_t readInt(std::istream& in)
    {
        union { int32_t number; char bytes[INT_SIZE]; } result;
        in.read(result.bytes, INT_SIZE);
        if (!in) throw FATALERROR("File error while reading integer value");
        return result.number;
    }

    // state constants (zero or positive values indicate cell index in current block)
    const int STATE_BEFORE_DATA = -1;
    const int STATE_TOPLEVEL_NONLEAF = -2;
    const int STATE_REFINEMENT_NONLEAF = -3;
    const int STATE_BLOCK_NONLEAF = -4;
    const int STATE_AFTER_DATA = -5;
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshAmrvacFile::open()
{
    // copy the number of mesh cells at the coarsest level from the user configuration
    _nxlone[0] = _levelOneX;
    _nxlone[1] = _levelOneY;
    _nxlone[2] = _levelOneZ;

    // open the data file
    string filepath = find<FilePaths>()->input(filename());
    _infile = System::ifstream(filepath);
    if (!_infile.is_open()) throw FATALERROR("Could not open the adaptive mesh data file " + filepath);
    find<Log>()->info("Reading adaptive mesh data from MPI-AMRVAC file " + filepath + "...");

    // read the parameters at the end of the file (starting at EOF minus 7 integers and 1 double)
    _infile.seekg(- 7*INT_SIZE - 1*DOUBLE_SIZE, std::ios_base::end);
    _nblocks = readInt(_infile);    // number of active tree leafs 'nleafs' (= #blocks)
               readInt(_infile);    // maximal refinement level present 'levmax'
    _ndims   = readInt(_infile);    // dimensionality 'NDIM'
               readInt(_infile);    // number of vector components 'NDIR'
    _nvars   = readInt(_infile);    // number of variables 'nw'
    int pars = readInt(_infile);    // number of equation-specific variables 'neqpar+nspecialpar'

    // read the block size in each dimension (before equation-specific variables)
    _infile.seekg(- 7*INT_SIZE - 1*DOUBLE_SIZE - _ndims*INT_SIZE - pars*DOUBLE_SIZE, std::ios_base::end);
    _nx[0] = _nx[1] = _nx[2] = 1;   // provide default of one for missing dimensions
    for (int i=0; i<_ndims; i++) _nx[i] = readInt(_infile);

    // calculate some handy grid characteristics:
    //    number of blocks at the coarsest level
    for (int i=0; i<3; i++)
    {
        if (_nxlone[i]%_nx[i])
            throw FATALERROR("Number of cells at the coarsest level is not a multiple of block size");
        _ng[i] = _nxlone[i]/_nx[i];
    }
    //    refinement factor: always 2 except for missing dimensions
    _nr[0] = _nr[1] = _nr[2] = 1;
    for (int i=0; i<_ndims; i++) _nr[i] = 2;
    //    number of cells in a block
    _ncells = _nx[0]*_nx[1]*_nx[2];
    //    block size in bytes
    _blocksize = _ncells*_nvars*DOUBLE_SIZE;

    // read the forest representing the grid structure (just after the data blocks)
    // there are are exactly _nblocks "true" values in the forest, but there can be many additional "false" values
    _infile.seekg(static_cast<int64_t>(_nblocks)*static_cast<int64_t>(_blocksize));
    _forest.clear();
    _forest.reserve(_nblocks);
    for (int i=0; i<_nblocks; i++)
    {
        while (true)
        {
            bool leaf = readInt(_infile) ? true : false;
            _forest.push_back(leaf);
            if (leaf) break;
        }
    }
    _forestsize = _forest.size();

    // set the context at the beginning of the file with no current record
    _infile.seekg(0);
    _state = STATE_BEFORE_DATA;
    _block.resize(_blocksize);
    _forestindex = 0;
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshAmrvacFile::close()
{
    _infile.close();
    _forest.clear();
    _block.clear();
    _state = STATE_AFTER_DATA;
}

//////////////////////////////////////////////////////////////////////

bool AdaptiveMeshAmrvacFile::read()
{
    if (_state == STATE_BEFORE_DATA)
    {
        _state = STATE_TOPLEVEL_NONLEAF;
        return true;
    }
    if (_state == STATE_TOPLEVEL_NONLEAF || _state == STATE_REFINEMENT_NONLEAF)
    {
        if (_forestindex < _forestsize)
        {
            _state = _forest[_forestindex++] ? STATE_BLOCK_NONLEAF : STATE_REFINEMENT_NONLEAF;
            return true;
        }
        else
        {
            _state = STATE_AFTER_DATA;
            return false;
        }
    }
    if (_state == STATE_BLOCK_NONLEAF)
    {
        // !!! this statement assumes that the incoming data has native byte order !!!
        _infile.read(reinterpret_cast<char*>(&_block[0]), _blocksize);
        if (!_infile) throw FATALERROR("File error while reading cell data");
        _state = 0;
        return true;
    }
    if (_state >= 0)
    {
        _state++;
        if (_state >= _ncells)
        {
            if (_forestindex < _forestsize)
            {
                _state = _forest[_forestindex++] ? STATE_BLOCK_NONLEAF : STATE_REFINEMENT_NONLEAF;
                return true;
            }
            else
            {
                _state = STATE_AFTER_DATA;
                return false;
            }
        }
        return true;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////

bool AdaptiveMeshAmrvacFile::isNonLeaf() const
{
    return _state < 0;
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshAmrvacFile::numChildNodes(int &nx, int &ny, int &nz) const
{
    if (_state == STATE_TOPLEVEL_NONLEAF)
    {
        nx = _ng[0]; ny = _ng[1]; nz = _ng[2];
    }
    else if (_state == STATE_REFINEMENT_NONLEAF)
    {
        nx = _nr[0]; ny = _nr[1]; nz = _nr[2];
    }
    else if (_state == STATE_BLOCK_NONLEAF)
    {
        nx = _nx[0]; ny = _nx[1]; nz = _nx[2];
    }

    // we expect three positive integers
    if (nx<1 || ny<1 || nz<1) throw FATALERROR("Invalid nonleaf information in mesh data");
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshAmrvacFile::value(int g) const
{
    // verify index range
    if (g < 0) throw FATALERROR("Field index out of range");
    if (g >= _nvars) throw FATALERROR("Insufficient number of field values in mesh data");
    if (_state < 0) throw FATALERROR("Invocation of value function for nonleaf node");

    return _block[g*_ncells+_state];
}

//////////////////////////////////////////////////////////////////////
