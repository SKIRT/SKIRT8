/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiFrameInstrument.hpp"
#include "FatalError.hpp"
#include "InstrumentFrame.hpp"
#include "PhotonPackage.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void MultiFrameInstrument::setupSelfBefore()
{
    DistantInstrument::setupSelfBefore();

    // verify attribute values
    if (_frames.size() != static_cast<size_t>(find<WavelengthGrid>()->numWavelengths()))
        throw FATALERROR("Number of instrument frames must equal number of wavelengths");
}

////////////////////////////////////////////////////////////////////

void MultiFrameInstrument::detect(PhotonPackage* pp)
{
    _frames[pp->ell()]->detect(pp);
}

////////////////////////////////////////////////////////////////////

void MultiFrameInstrument::write()
{
    int Nlambda = _frames.size();
    for (int ell=0; ell<Nlambda; ell++)
    {
        _frames[ell]->calibrateAndWriteData(ell);
    }
}

////////////////////////////////////////////////////////////////////
