/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AmHydrocarbonGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

AmHydrocarbonGrainComposition::AmHydrocarbonGrainComposition(SimulationItem *parent, double bulkDensity)
{
    parent->addChild(this);
    setBulkDensity(bulkDensity);  // must be set before calling loadLogHeatCapacityGrid() !!
    setup();
}

//////////////////////////////////////////////////////////////////////

void AmHydrocarbonGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    if (!bulkDensity()) setBulkDensity(1600.);
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
