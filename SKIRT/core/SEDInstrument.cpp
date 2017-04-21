/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SEDInstrument.hpp"
#include "LockFree.hpp"
#include "PhotonPackage.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void SEDInstrument::setupSelfBefore()
{
    DistantInstrument::setupSelfBefore();

    int Nlambda = find<WavelengthGrid>()->numWavelengths();
    _Ftotv.resize(Nlambda);
}

////////////////////////////////////////////////////////////////////

void SEDInstrument::detect(PhotonPackage* pp)
{
    int ell = pp->ell();
    double L = pp->luminosity();
    double taupath = opticalDepth(pp);
    double extf = exp(-taupath);
    double Lextf = L*extf;

    LockFree::add(_Ftotv[ell], Lextf);
}

////////////////////////////////////////////////////////////////////

void SEDInstrument::write()
{
    // construct a list of SED array pointers and the corresponding column names
    vector<Array*> Farrays({ &_Ftotv });
    vector<string> Fnames({ "total flux" });

    // sum the flux arrays element-wise across the different processes
    sumResults(Farrays);

    // calibrate and output the arrays
    calibrateAndWriteSEDs(Farrays, Fnames);
}

////////////////////////////////////////////////////////////////////
