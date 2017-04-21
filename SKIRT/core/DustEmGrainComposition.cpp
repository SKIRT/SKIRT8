/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustEmGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

DustEmGrainComposition::DustEmGrainComposition(SimulationItem *parent, string graintype, double rhobulk)
    : _grainType(graintype), _bulkMassDensity(rhobulk)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

void DustEmGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    setBulkDensity(_bulkMassDensity);
    loadLogHeatCapacityGrid("GrainComposition/DustEM/hcap/C_" + _grainType + ".DAT");
    loadOpticalGrid("GrainComposition/DustEM/oprop/LAMBDA.DAT",
                    "GrainComposition/DustEM/oprop/Q_" + _grainType + ".DAT",
                    "GrainComposition/DustEM/oprop/G_" + _grainType + ".DAT");
}

//////////////////////////////////////////////////////////////////////

string DustEmGrainComposition::name() const
{
    return "DustEM_" + _grainType;
}

//////////////////////////////////////////////////////////////////////
