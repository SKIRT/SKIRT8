/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Instrument.hpp"
#include "DustSystem.hpp"
#include "Log.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "PhotonPackage.hpp"
#include "TimeLogger.hpp"

////////////////////////////////////////////////////////////////////

void Instrument::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // get a pointer to the dust system, if present, without performing setup
    _ds = find<DustSystem>(false);
}

////////////////////////////////////////////////////////////////////

void Instrument::sumResults(const vector<Array*>& arrays)
{
    PeerToPeerCommunicator* comm = find<PeerToPeerCommunicator>();

    Log* log = find<Log>();
    TimeLogger logger(log->verbose() && comm->isMultiProc() ? log : 0, "communication of the observed fluxes");

    for (Array* arr : arrays) if (arr->size()) comm->sum(*arr);
}

////////////////////////////////////////////////////////////////////

double Instrument::opticalDepth(PhotonPackage* pp, double distance) const
{
    return _ds ? _ds->opticalDepth(pp,distance) : 0;
}

////////////////////////////////////////////////////////////////////
