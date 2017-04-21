/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SimpleInstrument.hpp"
#include "LockFree.hpp"
#include "PhotonPackage.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void SimpleInstrument::setupSelfBefore()
{
    SingleFrameInstrument::setupSelfBefore();

    WavelengthGrid* wavelengthGrid = find<WavelengthGrid>();
    int Nlambda = wavelengthGrid->numWavelengths();

    _Ftotv.resize(Nlambda);
    _ftotv.initialize("Instrument " + instrumentName() + " total flux", _Nframep, this);
}

////////////////////////////////////////////////////////////////////

void SimpleInstrument::detect(PhotonPackage* pp)
{
    int l = pixelOnDetector(pp);
    int ell = pp->ell();
    double L = pp->luminosity();
    double taupath = opticalDepth(pp);
    double extf = exp(-taupath);
    double Lextf = L*extf;

    LockFree::add(_Ftotv[ell], Lextf);
    if (l>=0)
    {
        LockFree::add(_ftotv(ell,l), Lextf);
    }
}

////////////////////////////////////////////////////////////////////

void SimpleInstrument::write()
{
    // construct a list of SED array pointers and the corresponding column names
    vector<Array*> Farrays({ &_Ftotv });
    vector<string> Fnames({ "total flux" });

    // sum the SED arrays element-wise across the different processes, and calibrate and output the result
    sumResults(Farrays);
    calibrateAndWriteSEDs(Farrays, Fnames);

    // sum the flux data cubes element-wise across the different processes
    std::shared_ptr<Array> completeCube = _ftotv.constructCompleteCube();

    // construct a list of data cube pointers and the corresponding file names
    vector<Array*> farrays({ completeCube.get() });
    vector<string> fnames({ "total" });

    // calibrate and output the result
    calibrateAndWriteDataCubes(farrays, fnames);
}

////////////////////////////////////////////////////////////////////
