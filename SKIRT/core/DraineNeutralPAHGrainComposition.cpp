/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DraineNeutralPAHGrainComposition.hpp"
#include "DraineGraphiteGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

DraineNeutralPAHGrainComposition::DraineNeutralPAHGrainComposition(SimulationItem *parent)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

void DraineNeutralPAHGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    loadOpticalGrid(true, "GrainComposition/Draine/PAHneu_30.dat", true, false, true, false);
    calculateEnthalpyGrid(DraineGraphiteGrainComposition::enthalpyFunction);
    setBulkDensity(2.24e3);
}

//////////////////////////////////////////////////////////////////////

string DraineNeutralPAHGrainComposition::name() const
{
    return "Draine_Neutral_PAH";
}

//////////////////////////////////////////////////////////////////////
