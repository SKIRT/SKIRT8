/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ForsteriteGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

ForsteriteGrainComposition::ForsteriteGrainComposition(SimulationItem *parent, GrainType type)
{
    parent->addChild(this);
    _grainType = type;
    setup();
}

//////////////////////////////////////////////////////////////////////

void ForsteriteGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    // determine the bulk density and resource filenames based on the grain type
    double density = 0.;
    string heatfile;
    string opticalfile;
    switch (_grainType)
    {
    case GrainType::Crystalline:
        density = 3330.;
        heatfile = "GrainComposition/Min/C_aSil.DAT";
        opticalfile = "GrainComposition/Min/Forsterite_Suto2006.dat";
        break;
    case GrainType::Amorphous:
        density = 2190.;
        heatfile = "GrainComposition/Themis/C_CM_amFo10Fe30FeS.DAT";
        opticalfile = "GrainComposition/Themis/CM_amFo10Fe30FeS_Jones2013_SKIRT.dat";
        break;
    }
    setBulkDensity(density);
    loadLogHeatCapacityGrid(heatfile);
    loadOpticalGrid(true, opticalfile, false, false, false, false);
}

//////////////////////////////////////////////////////////////////////

string ForsteriteGrainComposition::name() const
{
    switch (_grainType)
    {
    case GrainType::Crystalline:
        return "Crystalline_Forsterite";
    case GrainType::Amorphous:
        return "Amorphous_Forsterite";
    }
    return string();
}

//////////////////////////////////////////////////////////////////////
