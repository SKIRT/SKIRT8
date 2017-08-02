/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AmHydrocarbonGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

AmHydrocarbonGrainComposition::AmHydrocarbonGrainComposition(SimulationItem *parent, double bulkDensity)
{
    parent->addChild(this);
    setup();
    setBulkDensity(bulkDensity);
}

//////////////////////////////////////////////////////////////////////

void AmHydrocarbonGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    setBulkDensity(1600.);
    loadLogHeatCapacityGrid("GrainComposition/DustEM/hcap/C_CM20.DAT");
    loadOpticalGrid("GrainComposition/DustEM/oprop/LAMBDA.DAT",
                    "GrainComposition/DustEM/oprop/Q_CM20.DAT",
                    "GrainComposition/DustEM/oprop/G_CM20.DAT");
}

//////////////////////////////////////////////////////////////////////

string AmHydrocarbonGrainComposition::name() const
{
    return "Amorphous_Hydrocarbon";
}

//////////////////////////////////////////////////////////////////////
