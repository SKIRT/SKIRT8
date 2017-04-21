/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INSTRUMENTSYSTEM_HPP
#define INSTRUMENTSYSTEM_HPP

#include "Instrument.hpp"

//////////////////////////////////////////////////////////////////////

/** An InstrumentSystem instance keeps a list of zero or more instruments. The instruments can be
    of various nature (e.g. photometric, spectroscopic,...) and do not need to be located at the
    same observing position. */
class InstrumentSystem : public SimulationItem
{
    ITEM_CONCRETE(InstrumentSystem, SimulationItem, "an instrument system")

    PROPERTY_ITEM_LIST(instruments, Instrument, "the instruments")
        ATTRIBUTE_DEFAULT_VALUE(instruments, "SimpleInstrument")
        ATTRIBUTE_OPTIONAL(instruments)

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function writes down the results of the instrument system. It calls the write()
        function for each of the instruments. */
    void write();
};

////////////////////////////////////////////////////////////////////

#endif
