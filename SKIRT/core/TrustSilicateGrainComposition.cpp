/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TrustSilicateGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

TrustSilicateGrainComposition::TrustSilicateGrainComposition(SimulationItem *parent)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

void TrustSilicateGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    loadOpticalGrid(true, "GrainComposition/Trust/suvSil_121_1201.dat", false, true, false, true);
    loadEnthalpyGrid(true, "GrainComposition/Trust/Silicate_Calorimetry_1000.dat");
    setBulkDensity(3.5e3);
}

//////////////////////////////////////////////////////////////////////

string TrustSilicateGrainComposition::name() const
{
    return "Trust_Silicate";
}

//////////////////////////////////////////////////////////////////////
