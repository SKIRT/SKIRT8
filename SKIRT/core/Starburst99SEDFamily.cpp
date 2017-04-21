/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Starburst99SEDFamily.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "FITSInOut.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // number of items in the library read by the constructor
    const int Nlambda = 1221;
    const int NZ = 25;
    const int Nt = 308;
}

///////////////////////////////////////////////////////////////////

Starburst99SEDFamily::Starburst99SEDFamily(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

///////////////////////////////////////////////////////////////////

void Starburst99SEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    // get the resource file path
    string filepath = FilePaths::resource("SED/Starburst99/Patrik-imfKroupa-Zmulti-ml.fits.gz");
    Log* log = find<Log>();
    log->info("Reading Starburst99 data from file " + filepath + "...");

    // read the wavelength, metallicity, and age vectors from the library
    FITSInOut::readColumn(filepath+"[AXES][col lambda]", _lambdav, Nlambda);
    FITSInOut::readColumn(filepath+"[AXES][col metallicity]", _Zv, NZ);
    FITSInOut::readColumn(filepath+"[AXES][col time]", _tv, Nt);

    // read the emissivity data cube from the library
    Array data;
    int nx, ny, nz;
    FITSInOut::read(filepath+"[SED]", data, nx, ny, nz);
    if (nx!=Nt || ny!=NZ || nz!=Nlambda) throw FATALERROR("Starburst99 library data size does not match expectations");
    log->info("File " + filepath + " closed.");

    // copy the emissivity data into the table
    _jvv.resize(Nt, NZ, Nlambda);
    size_t i = 0;
    for (size_t k = 0; k!=Nlambda; ++k)
    {
        for (size_t m = 0; m!=NZ; ++m)
        {
            for (size_t p = 0; p!=Nt; ++p)
            {
                _jvv(p, m, k) = pow(10., data[i++]);  // stored in file as log10
            }
        }
    }

    // cache the simulation's wavelength grid
    _lambdagrid = find<WavelengthGrid>();
}

//////////////////////////////////////////////////////////////////////

Array Starburst99SEDFamily::luminosities(double M, double Z, double t, double z) const
{
    // find the appropriate SED from interpolating in the library
    int mL, mR;
    double hZ = 0.0;
    if (Z<=_Zv[0])
        mL = mR = 0;
    else if (Z>=_Zv[NZ-1])
        mL = mR = NZ-1;
    else
    {
        mL = NR::locateClip(_Zv,Z);
        mR = mL+1;
        double ZL = _Zv[mL];
        double ZR = _Zv[mR];
        hZ = (Z-ZL)/(ZR-ZL);
    }
    int pL, pR;
    double ht = 0.0;
    if (t<=_tv[0])
        pL = pR = 0;
    else if (t>=_tv[Nt-1])
        pL = pR = Nt-1;
    else
    {
        pL = NR::locateClip(_tv,t);
        pR = pL+1;
        double tL = _tv[pL];
        double tR = _tv[pR];
        ht = (t-tL)/(tR-tL);
    }
    const Array& jLLv = _jvv(pL,mL);
    const Array& jLRv = _jvv(pL,mR);
    const Array& jRLv = _jvv(pR,mL);
    const Array& jRRv = _jvv(pR,mR);
    Array jv(Nlambda);
    for (int k=0; k<Nlambda; k++)
        jv[k] = (1.0-ht)*(1.0-hZ)*jLLv[k]
                + (1.0-ht)*hZ*jLRv[k]
                + ht*(1.0-hZ)*jRLv[k]
                + ht*hZ*jRRv[k];

    // resample to the possibly redshifted simulation wavelength grid,
    // convert emissivities to luminosities (i.e. multiply by the wavelength bins),
    // multiply by the mass of the population (in solar masses),
    // and return the result
    return NR::resample<NR::interpolateLogLog>(_lambdagrid->lambdav()*(1-z), _lambdav, jv)
                                    * _lambdagrid->dlambdav() * M;
}

//////////////////////////////////////////////////////////////////////

int Starburst99SEDFamily::numParams() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

Array Starburst99SEDFamily::luminositiesGeneric(const Array& params, int skipvals, double z) const
{
    return luminosities(params[skipvals], params[skipvals+1], params[skipvals+2], z);
}

//////////////////////////////////////////////////////////////////////

double Starburst99SEDFamily::massGeneric(const Array& params, int skipvals) const
{
    return params[skipvals];
}

//////////////////////////////////////////////////////////////////////

string Starburst99SEDFamily::sourceName() const
{
    return "star";
}

//////////////////////////////////////////////////////////////////////

string Starburst99SEDFamily::sourceDescription() const
{
    return "star";
}

//////////////////////////////////////////////////////////////////////
