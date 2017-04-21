/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DraineIonizedPAHGrainComposition.hpp"
#include "DraineGraphiteGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

DraineIonizedPAHGrainComposition::DraineIonizedPAHGrainComposition(SimulationItem *parent)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

void DraineIonizedPAHGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    loadOpticalGrid(true, "GrainComposition/Draine/PAHion_30.dat", true, false, true, false);
    calculateEnthalpyGrid(DraineGraphiteGrainComposition::enthalpyFunction);
    setBulkDensity(2.24e3);
}

//////////////////////////////////////////////////////////////////////

string DraineIonizedPAHGrainComposition::name() const
{
    return "Draine_Ionized_PAH";
}

//////////////////////////////////////////////////////////////////////
