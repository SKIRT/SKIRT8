/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FitSchemeItemRegistry.hpp"
#include "ItemRegistry.hpp"
#include "SkirtUnitDef.hpp"

// ---> add new items below in alphabetical order

#include "AdjustableSkirtSimulation.hpp"
#include "ExtragalacticUnits.hpp"
#include "FitsKernel.hpp"
#include "GaussianKernel.hpp"
#include "OligoFitScheme.hpp"
#include "Optimization.hpp"
#include "ParameterRange.hpp"
#include "ParameterRanges.hpp"
#include "ReferenceImage.hpp"
#include "ReferenceImages.hpp"
#include "SIUnits.hpp"
#include "StellarUnits.hpp"

////////////////////////////////////////////////////////////////////

FitSchemeItemRegistry::FitSchemeItemRegistry(string version, string format)
{
    // start a new schema
    ItemRegistry::beginSchema("FitSKIRT", "a FitSKIRT parameter file", version, "fski",
                              "skirt-simulation-hierarchy", "FitScheme", format,
                              "http://www.skirt.ugent.be/skirt");

    // add the FitSKIRT unit definitions
    ItemRegistry::addUnitDef<SkirtUnitDef>();

    // add the SKIRT units classes
    ItemRegistry::add<SimulationItem>();
    ItemRegistry::add<Units>();
    ItemRegistry::add<SIUnits>();
    ItemRegistry::add<StellarUnits>();
    ItemRegistry::add<ExtragalacticUnits>();

    // add the FitSKIRT fit scheme items
    // ---> add new items in the order you want them to appear in choice lists for the user

    ItemRegistry::add<FitScheme>();
    ItemRegistry::add<OligoFitScheme>();
    ItemRegistry::add<AdjustableSkirtSimulation>();

    ItemRegistry::add<Optimization>();
    ItemRegistry::add<ParameterRanges>();
    ItemRegistry::add<ParameterRange>();
    ItemRegistry::add<ReferenceImages>();
    ItemRegistry::add<ReferenceImage>();

    ItemRegistry::add<ConvolutionKernel>();
    ItemRegistry::add<GaussianKernel>();
    ItemRegistry::add<FitsKernel>();
}

////////////////////////////////////////////////////////////////////

const SchemaDef* FitSchemeItemRegistry::getSchemaDef()
{
    return ItemRegistry::getSchemaDef("FitSKIRT");
}

////////////////////////////////////////////////////////////////////

FitSchemeItemRegistry::~FitSchemeItemRegistry()
{
    ItemRegistry::finalize();
}

////////////////////////////////////////////////////////////////////

