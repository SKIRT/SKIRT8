/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FullInstrument.hpp"
#include "DustEmissivity.hpp"
#include "LockFree.hpp"
#include "PanDustSystem.hpp"
#include "PhotonPackage.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void FullInstrument::setupSelfBefore()
{
    SingleFrameInstrument::setupSelfBefore();

    // determine whether the simulation contains a dust system
    _dustsystem = find<DustSystem>(false) != nullptr;

    if (_dustsystem)
    {
        // determine whether the simulation includes polarization;
        // a dust mix knows whether it supports polarization only after it has been setup,
        // so we need to fully setup the dust system before querying it
        _polarization = find<DustSystem>()->polarization();

        // determine whether the simulation includes dust emission;
        // which is supported (and known) only by a panchromatic dust system
        auto pds = find<PanDustSystem>(false);
        _dustemission = pds && pds->hasDustEmission();
    }

    // resize the detector arrays only when meaningful
    int Nlambda = find<WavelengthGrid>()->numWavelengths();
    string name = "Instrument " + instrumentName() + " ";

    _ftrav.initialize(name + "transparent flux", _Nframep, this);
    _Ftrav.resize(Nlambda);
    if (_dustsystem)
    {
        _fstrdirv.initialize(name + "direct stellar flux", _Nframep, this);
        _Fstrdirv.resize(Nlambda);
        _fstrscav.initialize(name + "scattered stellar flux", _Nframep, this);
        _Fstrscav.resize(Nlambda);
        if (_dustemission)
        {
            _fdusdirv.initialize(name + "direct dust flux", _Nframep, this);
            _Fdusdirv.resize(Nlambda);
            _fdusscav.initialize(name + "scattered dust flux", _Nframep, this);
            _Fdusscav.resize(Nlambda);
        }
        if (_numScatteringLevels > 0)
        {
            _fstrscavv.resize(_numScatteringLevels);
            for (auto& cube : _fstrscavv) cube.initialize(name + "scattered level flux", _Nframep, this);
            _Fstrscavv.resize(_numScatteringLevels, Nlambda);
        }
        if (_polarization)
        {
            _ftotQv.initialize(name + "Stokes Q", _Nframep, this);
            _FtotQv.resize(Nlambda);
            _ftotUv.initialize(name + "Stokes U", _Nframep, this);
            _FtotUv.resize(Nlambda);
            _ftotVv.initialize(name + "Stokes V", _Nframep, this);
            _FtotVv.resize(Nlambda);
        }
    }
}

////////////////////////////////////////////////////////////////////

void FullInstrument::detect(PhotonPackage* pp)
{
    int nscatt = pp->numScatt();
    int l = pixelOnDetector(pp);
    int ell = pp->ell();
    double L = pp->luminosity();
    double taupath = opticalDepth(pp);
    double extf = exp(-taupath);
    double Lextf = L*extf;

    // SEDs
    if (pp->isStellar())
    {
        if (nscatt==0)
        {
            LockFree::add(_Ftrav[ell], L);
            if (_dustsystem) LockFree::add(_Fstrdirv[ell], Lextf);
        }
        else
        {
            LockFree::add(_Fstrscav[ell], Lextf);
            if (nscatt<=_numScatteringLevels) LockFree::add(_Fstrscavv[nscatt-1][ell], Lextf);
        }
    }
    else
    {
        if (nscatt==0) LockFree::add(_Fdusdirv[ell], Lextf);
        else LockFree::add(_Fdusscav[ell], Lextf);
    }
    if (_polarization)
    {
        LockFree::add(_FtotQv[ell], Lextf*pp->stokesQ());
        LockFree::add(_FtotUv[ell], Lextf*pp->stokesU());
        LockFree::add(_FtotVv[ell], Lextf*pp->stokesV());
    }

    // frames
    if (l>=0)
    {
        if (pp->isStellar())
        {
            if (nscatt==0)
            {
                LockFree::add(_ftrav(ell,l), L);
                if (_dustsystem)
                    LockFree::add(_fstrdirv(ell,l), Lextf);
            }
            else
            {
                LockFree::add(_fstrscav(ell,l), Lextf);
                if (nscatt<=_numScatteringLevels)
                    LockFree::add(_fstrscavv[nscatt-1](ell,l), Lextf);
            }
        }
        else
        {
            if (nscatt==0)
                LockFree::add(_fdusdirv(ell,l), Lextf);
            else
                LockFree::add(_fdusscav(ell,l), Lextf);
        }
        if (_polarization)
        {
            LockFree::add(_ftotQv(ell,l), Lextf*pp->stokesQ());
            LockFree::add(_ftotUv(ell,l), Lextf*pp->stokesU());
            LockFree::add(_ftotVv(ell,l), Lextf*pp->stokesV());
        }
    }
}

////////////////////////////////////////////////////////////////////

void FullInstrument::write()
{
    // collect all the (or only the necessary) cubes
    std::shared_ptr<Array> ftravComp = _ftrav.constructCompleteCube();
    std::shared_ptr<Array> fstrdirvComp;
    std::shared_ptr<Array> fstrscavComp;
    std::shared_ptr<Array> fdusdirvComp;
    std::shared_ptr<Array> fdusscavComp;
    std::vector<std::shared_ptr<Array>> fstrscavvComp;
    std::shared_ptr<Array> ftotQvComp;
    std::shared_ptr<Array> ftotUvComp;
    std::shared_ptr<Array> ftotVvComp;

    if (_dustsystem)
    {
        fstrdirvComp = _fstrdirv.constructCompleteCube();
        fstrscavComp = _fstrscav.constructCompleteCube();
        if (_dustemission)
        {
            fdusdirvComp = _fdusdirv.constructCompleteCube();
            fdusscavComp = _fdusscav.constructCompleteCube();
        }
        if (_numScatteringLevels > 0)
        {
            fstrscavvComp.reserve(_numScatteringLevels);
            for (auto& cube : _fstrscavv) fstrscavvComp.push_back(cube.constructCompleteCube());
        }
        if (_polarization)
        {
            ftotQvComp = _ftotQv.constructCompleteCube();
            ftotUvComp = _ftotUv.constructCompleteCube();
            ftotVvComp = _ftotVv.constructCompleteCube();
        }
    }

    // compute the total flux and the total dust flux in temporary arrays
    Array ftotv;
    Array Ftotv;
    Array ftotdusv;
    Array Ftotdusv;
    if (_dustemission)
    {
        ftotv = *fstrdirvComp + *fstrscavComp + *fdusdirvComp + *fdusscavComp;
        Ftotv = _Fstrdirv + _Fstrscav + _Fdusdirv + _Fdusscav;
        ftotdusv = *fdusdirvComp + *fdusscavComp;
        Ftotdusv = _Fdusdirv + _Fdusscav;
    }
    else if (_dustsystem)
    {
        ftotv = *fstrdirvComp + *fstrscavComp;
        Ftotv = _Fstrdirv + _Fstrscav;
    }
    else
    {
        // don't output transparent frame separately because it is identical to the total frame
        ftotv = *ftravComp;
        ftravComp->resize(0);
        // do output integrated fluxes to avoid confusing zeros
        Ftotv = _Ftrav;
        _Fstrdirv = _Ftrav;
    }

    // construct list of SED array pointers and the corresponding column names
    vector<Array*> Farrays({ &Ftotv, &_Fstrdirv, &_Fstrscav, &Ftotdusv, &_Fdusscav, &_Ftrav });
    vector<string> Fnames({"total flux", "direct stellar flux", "scattered stellar flux",
                            "total dust emission flux", "dust emission scattered flux", "transparent flux" });
    if (_polarization)
    {
        Farrays.push_back(&_FtotQv);  Fnames.push_back("total Stokes Q");
        Farrays.push_back(&_FtotUv);  Fnames.push_back("total Stokes U");
        Farrays.push_back(&_FtotVv);  Fnames.push_back("total Stokes V");
    }
    if (_dustsystem)
    {
        for (int nscatt=0; nscatt<_numScatteringLevels; nscatt++)
        {
            Farrays.push_back( &(_Fstrscavv[nscatt]) );
            Fnames.push_back(std::to_string(nscatt+1) + "-times scattered flux");
        }
    }

    // sum the SED arrays element-wise across the different processes, and calibrate and output the result
    sumResults(Farrays);
    calibrateAndWriteSEDs(Farrays, Fnames);

    // construct list of data cube pointers and the corresponding file names
    vector<Array*> farrays({ &ftotv, ftravComp.get() });
    vector<string> fnames({ "total", "transparent" });
    if (_dustsystem)
    {
        farrays.push_back(fstrdirvComp.get());  fnames.push_back("direct");
        farrays.push_back(fstrscavComp.get());  fnames.push_back("scattered");
        if (_dustemission)
        {
            farrays.push_back(&ftotdusv);           fnames.push_back("dust");
            farrays.push_back(fdusscavComp.get());  fnames.push_back("dustscattered");
        }
    }
    if (_polarization)
    {
        farrays.push_back(ftotQvComp.get());  fnames.push_back("stokesQ");
        farrays.push_back(ftotUvComp.get());  fnames.push_back("stokesU");
        farrays.push_back(ftotVvComp.get());  fnames.push_back("stokesV");
    }
    if (_dustsystem)
    {
        for (int nscatt=0; nscatt<_numScatteringLevels; nscatt++)
        {
            farrays.push_back(fstrscavvComp[nscatt].get());
            fnames.push_back("scatteringlevel" + std::to_string(nscatt+1));
        }
    }

    // calibrate and output the data cubes
    calibrateAndWriteDataCubes(farrays, fnames);
}

////////////////////////////////////////////////////////////////////
