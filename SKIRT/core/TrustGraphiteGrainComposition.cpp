/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TrustGraphiteGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

TrustGraphiteGrainComposition::TrustGraphiteGrainComposition(SimulationItem *parent)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

void TrustGraphiteGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    loadOpticalGrid(true, "GrainComposition/Trust/Gra_121_1201.dat", false, true, false, true);
    loadEnthalpyGrid(true, "GrainComposition/Trust/Graphitic_Calorimetry_1000.dat");
    setBulkDensity(2.24e3);
}

//////////////////////////////////////////////////////////////////////

string TrustGraphiteGrainComposition::name() const
{
    return "Trust_Graphite";
}

//////////////////////////////////////////////////////////////////////
