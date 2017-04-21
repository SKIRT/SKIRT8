/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FITSCHEMEITEMREGISTRY_HPP
#define FITSCHEMEITEMREGISTRY_HPP

#include "Basics.hpp"
class SchemaDef;

////////////////////////////////////////////////////////////////////

/** The FitSchemeItemRegistry class manages the registration of all discoverable SimulationItem
    subclasses that can reside in a fit scheme used by the \c FitSKIRT program and defined in the
    \c FitSKIRT/core project subdirectory). A single instance of the FitSchemeItemRegistry class
    must be constructed just after program startup, and certainly before any parallel threads are
    started. A good place is early in the main() function. The program should not use the exit() or
    abort() functions, but rather let the main() function run to normal completion and return an
    exit code. */
class FitSchemeItemRegistry final
{
public:
    /** The constructor registers all discoverable SimulationItem subclasses that can reside in a
        fit scheme with the item registry, including them in a single SMILE schema called
        'FitSKIRT'. The first argument \em version specifies the version of this schema definition.
        The second argument \em format specifies the version of the described data format (listed
        on the root element). The constructor is \em not thread-safe and may be called only during
        program startup from a single thread. */
    FitSchemeItemRegistry(string version, string format);

    /** This static function returns a pointer to the 'FitSKIRT' schema definition. Ownership
        remains with the registry. This function is thread-safe and may called at any time after
        construction and before destruction of the SimulationItemRegistry instance. */
    static const SchemaDef* getSchemaDef();

    /** The destructor releases the global memory managed by the item registry. */
    ~FitSchemeItemRegistry();
};

////////////////////////////////////////////////////////////////////

#endif
