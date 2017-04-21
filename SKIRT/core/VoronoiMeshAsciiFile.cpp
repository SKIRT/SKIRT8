/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoronoiMeshAsciiFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

void VoronoiMeshAsciiFile::open()
{
    // open the data file
    string filepath = find<FilePaths>()->input(filename());
    _infile = System::ifstream(filepath);
    if (!_infile.is_open()) throw FATALERROR("Could not open the Voronoi mesh data file " + filepath);
    find<Log>()->info("Reading Voronoi mesh data from ASCII file " + filepath + "...");
    _columns.clear();
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshAsciiFile::close()
{
    _infile.close();
    _columns.clear();
}

////////////////////////////////////////////////////////////////////

bool VoronoiMeshAsciiFile::read()
{
    // read the next line, splitting it in columns, and skip empty and comment lines
    while (_infile)
    {
        string line;
        std::getline(_infile, line);
        _columns = StringUtils::split(StringUtils::squeeze(line), " ");
        if (!_columns[0].empty() && !StringUtils::startsWith(_columns[0], "#"))
        {
            return true;
        }
    }
    _columns.clear();
    return false;
}

////////////////////////////////////////////////////////////////////

Vec VoronoiMeshAsciiFile::particle() const
{
    // verify index range
    if (_columns.size() < 3) throw FATALERROR("Insufficient number of particle coordinates in Voronoi mesh data");

    // get the coordinate values
    bool okx, oky, okz;
    double x = StringUtils::toDouble(_columns[0], &okx);
    double y = StringUtils::toDouble(_columns[1], &oky);
    double z = StringUtils::toDouble(_columns[2], &okz);
    if (!okx || !oky || !okz) throw FATALERROR("Invalid particle coordinate(s) in Voronoi mesh data");

    // convert to SI units
    return Vec(x*_coordinateUnits, y*_coordinateUnits, z*_coordinateUnits);
}

////////////////////////////////////////////////////////////////////

double VoronoiMeshAsciiFile::value(int g) const
{
    // verify index range
    if (g < 0) throw FATALERROR("Field index out of range");
    if (static_cast<size_t>(g+3) >= _columns.size())
        throw FATALERROR("Insufficient number of field values in Voronoi mesh data");

    // get the appropriate column value
    bool ok;
    double value = StringUtils::toDouble(_columns[g+3], &ok);
    if (!ok) throw FATALERROR("Invalid field value in Voronoi mesh data");
    return value;
}

////////////////////////////////////////////////////////////////////
