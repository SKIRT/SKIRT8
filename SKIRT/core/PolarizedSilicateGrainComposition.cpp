/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PolarizedSilicateGrainComposition.hpp"
#include "DraineSilicateGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

void PolarizedSilicateGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    setBulkDensity(3.0e3);
    calculateEnthalpyGrid(DraineSilicateGrainComposition::enthalpyFunction);
    loadPolarizedOpticalGrid(true, "Silicate_STOKES_Sxx.DAT");
}

//////////////////////////////////////////////////////////////////////

string PolarizedSilicateGrainComposition::name() const
{
    return "Polarized_Draine_Silicate";
}

//////////////////////////////////////////////////////////////////////
