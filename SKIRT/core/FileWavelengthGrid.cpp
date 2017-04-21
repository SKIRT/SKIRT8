/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

void FileWavelengthGrid::setupSelfBefore()
{
    PanWavelengthGrid::setupSelfBefore();

    string filename = find<FilePaths>()->input(_filename);
    std::ifstream file = System::ifstream(filename);
    if (! file.is_open()) throw FATALERROR("Could not open the data file " + filename);
    find<Log>()->info("Reading wavelength grid data from file " + filename + "...");

    int Nlambda;
    file >> Nlambda;
    Array lambdav(Nlambda);
    for (int k=0; k<Nlambda; k++)
    {
        file >> lambdav[k];
        lambdav[k] /= 1e6;   // conversion from micron to m
    }
    file.close();
    find<Log>()->info("File " + filename + " closed.");

    NR::sort(lambdav);
    setWavelengths(lambdav);
}

//////////////////////////////////////////////////////////////////////
