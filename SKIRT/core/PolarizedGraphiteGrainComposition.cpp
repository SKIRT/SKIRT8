/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PolarizedGraphiteGrainComposition.hpp"
#include "DraineGraphiteGrainComposition.hpp"


//////////////////////////////////////////////////////////////////////

void PolarizedGraphiteGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    setBulkDensity(2.24e3);
    calculateEnthalpyGrid(DraineGraphiteGrainComposition::enthalpyFunction);
    loadPolarizedOpticalGrid(true, "Graphite_STOKES_Sxx_001.DAT");
}

//////////////////////////////////////////////////////////////////////

string PolarizedGraphiteGrainComposition::name() const
{
    return "Polarized_Draine_Graphite";
}

//////////////////////////////////////////////////////////////////////
