/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshAsciiFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "System.hpp"
#include "StringUtils.hpp"

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshAsciiFile::open()
{
    // open the data file
    string filepath = find<FilePaths>()->input(filename());
    _infile = System::ifstream(filepath);
    if (!_infile.is_open()) throw FATALERROR("Could not open the adaptive mesh data file " + filepath);
    find<Log>()->info("Reading adaptive mesh data from ASCII file " + filepath + "...");
    _columns.clear();
    _isNonLeaf = false;
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshAsciiFile::close()
{
    _infile.close();
    _columns.clear();
    _isNonLeaf = false;
}

//////////////////////////////////////////////////////////////////////

bool AdaptiveMeshAsciiFile::read()
{
    // read the next line, splitting it in columns, and skip empty and comment lines
    while (_infile)
    {
        string line;
        std::getline(_infile, line);
        _columns = StringUtils::split(StringUtils::squeeze(line), " ");
        if (!_columns[0].empty() && !StringUtils::startsWith(_columns[0], "#"))
        {
            // remember the node type
            _isNonLeaf = StringUtils::startsWith(_columns[0], "!");
            if (_isNonLeaf)
            {
                // remove the exclamation mark
                if (_columns[0].size() > 1) _columns[0].erase(0,1);
                else _columns.erase(_columns.begin());
            }
            return true;
        }
    }
    _columns.clear();
    _isNonLeaf = false;
    return false;
}

//////////////////////////////////////////////////////////////////////

bool AdaptiveMeshAsciiFile::isNonLeaf() const
{
    return _isNonLeaf;
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshAsciiFile::numChildNodes(int &nx, int &ny, int &nz) const
{
    if (_columns.size() < 3) throw FATALERROR("Invalid nonleaf line in mesh data");

    // get the column values; illegal values default to zero
    nx = StringUtils::toInt(_columns[0]);
    ny = StringUtils::toInt(_columns[1]);
    nz = StringUtils::toInt(_columns[2]);

    // we expect three positive integers
    if (nx<1 || ny<1 || nz<1) throw FATALERROR("Invalid nonleaf line in mesh data");
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshAsciiFile::value(int g) const
{
    // verify index range
    if (g < 0) throw FATALERROR("Field index out of range");
    if (static_cast<size_t>(g) >= _columns.size()) throw FATALERROR("Insufficient number of field values in mesh data");

    // get the appropriate column value
    bool ok;
    double value = StringUtils::toDouble(_columns[g], &ok);
    if (!ok) throw FATALERROR("Invalid leaf line in mesh data");
    return value;
}

//////////////////////////////////////////////////////////////////////
