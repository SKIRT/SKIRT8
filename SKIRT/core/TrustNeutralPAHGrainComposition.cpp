/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TrustNeutralPAHGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

TrustNeutralPAHGrainComposition::TrustNeutralPAHGrainComposition(SimulationItem *parent)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

void TrustNeutralPAHGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    loadOpticalGrid(true, "GrainComposition/Trust/PAH_28_1201_neu.dat", false, true, false, true);
    loadEnthalpyGrid(true, "GrainComposition/Trust/Graphitic_Calorimetry_1000.dat");
    setBulkDensity(2.24e3);
}

//////////////////////////////////////////////////////////////////////

string TrustNeutralPAHGrainComposition::name() const
{
    return "Trust_Neutral_PAH";
}

//////////////////////////////////////////////////////////////////////
