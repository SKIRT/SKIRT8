/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MinSilicateGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

void MinSilicateGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    setBulkDensity(3.09e3);
    loadLogHeatCapacityGrid("GrainComposition/DustEM/hcap/C_aSil.DAT");
    loadOpticalGrid(true, "GrainComposition/Min/aSil_Min2007.dat", false, false, false, false);
}

//////////////////////////////////////////////////////////////////////

string MinSilicateGrainComposition::name() const
{
    return "Min_Silicate";
}

//////////////////////////////////////////////////////////////////////
