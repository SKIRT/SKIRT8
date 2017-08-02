/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EnstatiteGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

EnstatiteGrainComposition::EnstatiteGrainComposition(SimulationItem *parent, GrainType type)
{
    parent->addChild(this);
    _grainType = type;
    setup();
}

//////////////////////////////////////////////////////////////////////

void EnstatiteGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    // set the bulk density and load the appropriate resources based on the grain type
    switch (_grainType)
    {
    case GrainType::Crystalline:
        setBulkDensity(2800.);
        loadLogHeatCapacityGrid("GrainComposition/Min/C_aSil.DAT");
        loadOpticalGrid(true, "GrainComposition/Min/Enstatite_Jaeger1998.dat", false, false, false, false);
        break;

    case GrainType::Amorphous:
        setBulkDensity(2190.);
        loadLogHeatCapacityGrid("GrainComposition/DustEM/hcap/C_aPyM5.DAT");
        loadOpticalGrid("GrainComposition/DustEM/oprop/LAMBDA.DAT",
                        "GrainComposition/DustEM/oprop/Q_aPyM5.DAT",
                        "GrainComposition/DustEM/oprop/G_aPyM5.DAT");
        break;
    }
}

//////////////////////////////////////////////////////////////////////

string EnstatiteGrainComposition::name() const
{
    switch (_grainType)
    {
    case GrainType::Crystalline:
        return "Crystalline_Enstatite";
    case GrainType::Amorphous:
        return "Amorphous_Enstatite";
    }
    return string();
}

//////////////////////////////////////////////////////////////////////
