/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BuildInfo.hpp"
#include "FitSchemeItemRegistry.hpp"
#include "FitSkirtCommandLineHandler.hpp"
#include "ProcessManager.hpp"
#include "SignalHandler.hpp"
#include "SimulationItemRegistry.hpp"
#include "System.hpp"

//////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // Initialize inter-process communication capability, if present
    ProcessManager pm(&argc, &argv);

    // Initialize the system and install signal handlers
    System system(argc, argv);
    SignalHandler::InstallSignalHandlers();

    // Add all SKIRT simulation items and FitScheme items to the item registry
    string version = BuildInfo::projectVersion();
    SimulationItemRegistry registry1(version, "6.1");
    FitSchemeItemRegistry registry2(version, "6.1");

    // handle the command line arguments
    FitSkirtCommandLineHandler handler;
    return handler.perform();
}

//////////////////////////////////////////////////////////////////////
