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

    // set the bulk density and load the appropriate resources based on the grain type
    switch (_grainType)
    {
    case GrainType::Crystalline:
        setBulkDensity(3330.);
        loadLogHeatCapacityGrid("GrainComposition/Min/C_aSil.DAT");
        loadOpticalGrid(true, "GrainComposition/Min/Forsterite_Suto2006.dat", false, false, false, false);
        break;

    case GrainType::Amorphous:
        setBulkDensity(2190.);
        loadLogHeatCapacityGrid("GrainComposition/DustEM/hcap/C_aOLM5.DAT");
        loadOpticalGrid("GrainComposition/DustEM/oprop/LAMBDA.DAT",
                        "GrainComposition/DustEM/oprop/Q_aOLM5.DAT",
                        "GrainComposition/DustEM/oprop/G_aOLM5.DAT");
        break;
    }
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
