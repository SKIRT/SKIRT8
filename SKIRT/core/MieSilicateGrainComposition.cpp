/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MieSilicateGrainComposition.hpp"
#include "DraineSilicateGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

MieSilicateGrainComposition::MieSilicateGrainComposition(SimulationItem *parent)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

void MieSilicateGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    loadOpticalGrid(true, "GrainComposition/Other/MieAmorphousSilicate.dat", false, false, false, false);
    calculateEnthalpyGrid(DraineSilicateGrainComposition::enthalpyFunction);
    setBulkDensity(3.0e3);
}

//////////////////////////////////////////////////////////////////////

string MieSilicateGrainComposition::name() const
{
    return "Mie_Silicate";
}

//////////////////////////////////////////////////////////////////////
