/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FrameInstrument.hpp"
#include "LockFree.hpp"
#include "PhotonPackage.hpp"

////////////////////////////////////////////////////////////////////

void FrameInstrument::setupSelfBefore()
{
    SingleFrameInstrument::setupSelfBefore();

    _distftotv.initialize("Instrument " + instrumentName() + " total flux", _Nframep, this);
}

////////////////////////////////////////////////////////////////////

void FrameInstrument::detect(PhotonPackage* pp)
{
    int l = pixelOnDetector(pp);
    if (l >= 0)
    {
        int ell = pp->ell();
        double L = pp->luminosity();
        double taupath = opticalDepth(pp);
        double extf = exp(-taupath);
        double Lextf = L*extf;

        LockFree::add(_distftotv(ell,l), Lextf);
    }
}

////////////////////////////////////////////////////////////////////

void FrameInstrument::write()
{
    // sum the flux data cubes element-wise across the different processes
    std::shared_ptr<Array> completeCube = _distftotv.constructCompleteCube();

    // construct a list of data cube pointers and the corresponding file names
    vector<Array*> farrays({ completeCube.get() });
    vector<string> fnames({ "total" });

    // calibrate and output the result
    calibrateAndWriteDataCubes(farrays, fnames);
}

////////////////////////////////////////////////////////////////////
