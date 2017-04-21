/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FitScheme.hpp"
#include "FatalError.hpp"
#include "TimeLogger.hpp"

////////////////////////////////////////////////////////////////////

void FitScheme::setup()
{
    if (setupStarted()) return;

    _log->setup();
    TimeLogger logger(_log, "setup");
    SimulationItem::setup();
}

////////////////////////////////////////////////////////////////////

void FitScheme::run()
{
    // verify setup
    if (!setupStarted()) throw FATALERROR("Fit scheme has not been setup before being run");

    _communicator->acquireSlaves();
    if(_communicator->isMaster())
    {
        TimeLogger logger(_log, "fitting");
        runSelf();
    }
    _communicator->releaseSlaves();
}

////////////////////////////////////////////////////////////////////

void FitScheme::setupAndRun()
{
    _communicator->setup();
    string processInfo = _communicator->isMultiProc() ?
                             " with " + std::to_string(_communicator->size()) + " processes" : "";

    _log->setup();
    TimeLogger logger(_log, "fit scheme " + _paths->outputPrefix() + processInfo);

    setup();
    run();
}

////////////////////////////////////////////////////////////////////

MasterSlaveCommunicator* FitScheme::communicator() const
{
    return _communicator;
}

////////////////////////////////////////////////////////////////////

Log* FitScheme::log() const
{
    return _log;
}

////////////////////////////////////////////////////////////////////

FilePaths* FitScheme::filePaths() const
{
    return _paths;
}

////////////////////////////////////////////////////////////////////

void FitScheme::setParallelSimulationCount(int value)
{
    _parallelSimulations = value;
}

////////////////////////////////////////////////////////////////////

int FitScheme::parallelSimulationCount() const
{
    return _parallelSimulations;
}

////////////////////////////////////////////////////////////////////

void FitScheme::setParallelThreadCount(int value)
{
    _parallelThreads = value;
}

////////////////////////////////////////////////////////////////////

int FitScheme::parallelThreadCount() const
{
    return _parallelThreads;
}

////////////////////////////////////////////////////////////////////
