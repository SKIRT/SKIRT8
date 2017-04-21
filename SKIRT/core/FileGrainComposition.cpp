/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

void FileGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    setBulkDensity(_bulkMassDensity);
    loadOpticalGrid(false, _opticalFilename, false, false, false, false);
    loadEnthalpyGrid(false, _calorimetricFilename);
}

//////////////////////////////////////////////////////////////////////

string FileGrainComposition::name() const
{
    return "File_" + _opticalFilename + "_" + _calorimetricFilename;
}

//////////////////////////////////////////////////////////////////////
