/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FITSCHEME_HPP
#define FITSCHEME_HPP

#include "SimulationItem.hpp"
#include "ConsoleLog.hpp"
#include "FilePaths.hpp"
#include "MasterSlaveCommunicator.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

/** FitScheme is the abstract base class that represents a complete FitSKIRT fit scheme. The class
    inherits SimulationItem so that all functionality built for SKIRT simulation hierarchies
    (including the discovery mechanism) can also be used for fit schemes. A FitScheme instance sits
    at the top of a run-time fit scheme hierarchy (i.e. it has no parent). It holds a number of
    essential fit-scheme-wide property instances. Some of these (a system of units) are
    discoverable and hence fully user-configurable. The other properties (a file paths object, a
    logging mechanism, a master-slave communicator) are not discoverable. When a FitScheme instance
    is constructed, a default instance is created for each of these properties. A reference to
    these property instances can be retrieved through the corresponding getter, and in some cases,
    the property can be further configured under program control (e.g., to set the input and output
    file paths for the fit scheme). */
class FitScheme : public SimulationItem
{
    ITEM_ABSTRACT(FitScheme, SimulationItem, "the fit scheme")

    PROPERTY_ITEM(units, Units, "the units system")
        ATTRIBUTE_DEFAULT_VALUE(units, "ExtragalacticUnits")

    ITEM_END()

    //======== Construction - Setup - Run - Destruction  ===========

public:
    /** This function performs setup for the complete fit scheme hierarchy. It invokes the setup()
        function defined in the SimulationItem base class, surrounded by start/finish log messages.
        It is recommended to use the setupAndRun() function rather than setup() and run()
        separately. */
    void setup();

    /** This function performs the fit scheme by invoking the runSelf() function to be defined in a
        subclass, surrounded by start/finish log messages. The setup() function must have been
        called before invoking run(). It is recommended to use the setupAndRun() function rather
        than setup() and run() separately. */
    void run();

    /** This function performs setup and executes the fit scheme by invoking setup() and run() in
        succession. */
    void setupAndRun();

protected:
    /** This function actually runs the fit scheme, assuming that setup() has been already
        performed. Its implementation must be provided by a subclass. */
    virtual void runSelf() = 0;

    //======== Getters and Setters for Non-Discoverable Attributes =======

public:
    /** Returns the MasterSlaveCommunicator of this fit scheme. */
    MasterSlaveCommunicator* communicator() const;

    /** Returns the logging mechanism for this fit scheme. */
    Log* log() const;

    /** Returns the input/output file paths object for this fit scheme. */
    FilePaths* filePaths() const;

    /** Sets the number of SKIRT simulations performed in parallel by this fit scheme. By default,
        the number of parallel simulations is equal to one (i.e. simulations are serialized). */
    void setParallelSimulationCount(int value);

    /** Returns the number of SKIRT simulations performed in parallel by this fit scheme. */
    int parallelSimulationCount() const;

    /** Sets the number of parallel threads for each SKIRT simulation performed by this fit scheme.
        By default, the number of parallel threads per simulation is equal to one. */
    void setParallelThreadCount(int value);

    /** Returns the number of parallel threads for each SKIRT simulation performed by this fit
        scheme. */
    int parallelThreadCount() const;

    //======================== Data Members ========================

private:
    MasterSlaveCommunicator* _communicator{ new MasterSlaveCommunicator(this) };
    Log* _log{ new ConsoleLog(this) };
    FilePaths* _paths{ new FilePaths(this) };
    int _parallelSimulations{1};           // the number of parallel simulations
    int _parallelThreads{1};               // the number of parallel theads per simulation
};

////////////////////////////////////////////////////////////////////

#endif
