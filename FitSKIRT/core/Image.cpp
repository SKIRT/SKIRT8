/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Image.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "FITSInOut.hpp"
#include "Log.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "SimulationItem.hpp"
#include "StringUtils.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

Image::Image()
{
}

////////////////////////////////////////////////////////////////////

Image::Image(const SimulationItem* item, string filename, string directory, bool log)
{
    load(item, filename, directory, log);
}

////////////////////////////////////////////////////////////////////

void Image::load(const SimulationItem* item, string filename, string directory, bool log)
{
    // Determine the path of the input FITS file
    filename = StringUtils::addExtension(filename, "fits");
    string filepath;
    if (directory.empty())
    {
        filepath = item->find<FilePaths>()->input(filename);
    }
    else
    {
        if (!System::isDir(directory)) throw FATALERROR("Import directory does not exist: " + directory);
        filepath = StringUtils::joinPaths(directory, filename);
    }

    // Read the input data
    if (log) item->find<Log>()->info("Reading FITS file " + filepath);
    int nframes;
    FITSInOut::read(filepath, _data, _nx, _ny, nframes);

    // Verify that the FITS file contains only one frame
    if (nframes != 1) throw FATALERROR("FITS image contains multiple frames");

    // Log the dimensions of the image
    if (log) item->find<Log>()->info("Frame dimensions: " + std::to_string(_nx) + " x " + std::to_string(_ny));
}

////////////////////////////////////////////////////////////////////

void Image::resize(int nx, int ny)
{
    // Set the size in the x and y direction
    _nx = nx;
    _ny = ny;

    // Resize the pixel data
    _data.resize(nx*ny);
}

////////////////////////////////////////////////////////////////////

void Image::setData(const Array& data)
{
    _data = data;
}

////////////////////////////////////////////////////////////////////

void Image::moveData(Array&& data)
{
    _data = std::move(data);
}

////////////////////////////////////////////////////////////////////

void Image::save(const SimulationItem* item, string description, string filename,
                 double incx, double incy, string dataUnits, string xyUnits)
{
    // Try to find a PeerToPeerCommunicator object and ensure it is setup
    PeerToPeerCommunicator* comm = item->find<PeerToPeerCommunicator>(false);
    if (comm) comm->setup();

    // Only write the FITS file if this process is the root or no PeerToPeerCommunicator was found
    if (!comm || comm->isRoot())
    {
        // Determine the path of the output FITS file
        string filepath = item->find<FilePaths>()->output(StringUtils::addExtension(filename, "fits"));

        // Write the FITS file
        item->find<Log>()->info("Writing " + description + " to " + filepath + "...");
        FITSInOut::write(filepath, _data, _nx, _ny, 1, incx, incy, 0., 0., dataUnits, xyUnits);
    }
}

////////////////////////////////////////////////////////////////////
