/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FSPSVariableIMFSEDFamily.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // number of items in the library along each dimension
    const size_t Nalpha = 100;      // index s
    const size_t NZ = 22;           // index m
    const size_t Nt = 94;           // index p
    const size_t Nlambda = 1963;    // index k

    // key combining indices s and m
    size_t key(size_t s, size_t m) { return (s << 32) | m; }
}

///////////////////////////////////////////////////////////////////

FSPSVariableIMFSEDFamily::FSPSVariableIMFSEDFamily(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

///////////////////////////////////////////////////////////////////

namespace
{
    // returns "LoM" or "HiM" depending on the specified IMF type (as an enumeration value)
    string imfString(FSPSVariableIMFSEDFamily::InitialMassFunction type)
    {
        switch (type)
        {
        case FSPSVariableIMFSEDFamily::InitialMassFunction::VaryLowMass: return "LoM";
        case FSPSVariableIMFSEDFamily::InitialMassFunction::VaryHighMass: return "HiM";
        }
        return string();
    }
}

///////////////////////////////////////////////////////////////////

void FSPSVariableIMFSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    // prepare the vectors for the coordinate grid points along each dimension
    _alphav.resize(Nalpha);     // index s
    _Zv.resize(NZ);             // index m
    _tv.resize(Nt);             // index p
    _lambdav.resize(Nlambda);   // index k

    // load the coordinate grid points along each dimension
    find<Log>()->info("Reading grid points for variable IMF SED family...");
    string imf = imfString(initialMassFunction());

    // --> IMF slope (alpha)
    {
        string filename = FilePaths::externalResource("imf_slopes_" + imf + ".dat");
        std::ifstream file = System::ifstream(filename);
        if (!file.is_open()) throw FATALERROR("Could not open the external resource file " + filename);
        for (size_t s=0; s<Nalpha; ++s) file >> _alphav[s];
        // if needed, reverse the array so that the values are in increasing order
        _alphaReversed = (_alphav[0] > _alphav[1]);
        if (_alphaReversed) std::reverse(begin(_alphav), end(_alphav));
    }
    // --> metallicity (Z)
    {
        string filename = FilePaths::externalResource("zlegend.dat");
        std::ifstream file = System::ifstream(filename);
        if (!file.is_open()) throw FATALERROR("Could not open the external resource file " + filename);
        for (size_t m=0; m<NZ; ++m) file >> _Zv[m];
    }
    // --> wavelength (lambda) and age (t)
    {
        string filename = FilePaths::externalResource("SSP_" + imf + "_imf0000_z00.out.spec");
        std::ifstream file = System::ifstream(filename);
        if (!file.is_open()) throw FATALERROR("Could not open the external resource file " + filename);

        // skip header comments and one intro line
        string line;
        while (file.peek() == '#') getline(file,line);
        getline(file,line);

        // read the wavelengths (all on a single line) including the line ending
        for (size_t k=0; k<Nlambda; ++k)
        {
            double lambda;
            file >> lambda;
            _lambdav[k] = lambda * 1e-10;   // wavelength in file in A, we want in m
        }
        getline(file,line);

        // read the ages (first number on the line with a line in between)
        for (size_t p=0; p<Nt; ++p)
        {
            double logt;
            file >> logt;
            _tv[p] = pow(10,logt); // age in file in log10(yr), we want in yr
            getline(file,line);
            getline(file,line);
        }
    }
    find<Log>()->info("Done reading grid points for variable IMF SED family");

    // cache the simulation's wavelength grid
    _lambdagrid = find<WavelengthGrid>();
}

//////////////////////////////////////////////////////////////////////

namespace
{
    // returns properly clipped left and right grid indices and fractional distance for a given coordinate value;
    // if reversed is true; the indices are considered to be in reverse order
    void locateInGrid(const Array& grid, size_t N, bool reversed, double value, size_t& iL, size_t& iR, double& h)
    {
        if (value <= grid[0])
        {
            iL = iR = 0;
            h = 0;
        }
        else if (value >= grid[N-1])
        {
            iL = iR = N-1;
            h = 0;
        }
        else
        {
            iL = NR::locate(grid, value);
            iR = iL+1;
            h = (value-grid[iL])/(grid[iR]-grid[iL]);
        }

        if (reversed)
        {
            iL = N - 1 - iL;
            iR = N - 1 - iR;
        }
    }
}

//////////////////////////////////////////////////////////////////////

Array FSPSVariableIMFSEDFamily::luminosities(double M, double Z, double t, double alpha, double z) const
{
    // find the appropriate indices in the respective grids for each parameter dimension
    size_t sL, sR, mL, mR, pL, pR;
    double ha, hZ, ht;
    locateInGrid(_alphav, Nalpha, _alphaReversed, alpha, sL, sR, ha);
    locateInGrid(_Zv,     NZ,     false,          Z,     mL, mR, hZ);
    locateInGrid(_tv,     Nt,     false,          t,     pL, pR, ht);

    // get the emissivity tables bracketing the requested slope and metallicity
    const ArrayTable<2>& jLLvv = getEmissivityTable(sL, mL);
    const ArrayTable<2>& jLRvv = getEmissivityTable(sL, mR);
    const ArrayTable<2>& jRLvv = getEmissivityTable(sR, mL);
    const ArrayTable<2>& jRRvv = getEmissivityTable(sR, mR);

    // get the emissivity tables bracketing the requested slope, metallicity, and age
    const Array& jLLLv = jLLvv[pL];
    const Array& jLLRv = jLLvv[pR];
    const Array& jLRLv = jLRvv[pL];
    const Array& jLRRv = jLRvv[pR];
    const Array& jRLLv = jRLvv[pL];
    const Array& jRLRv = jRLvv[pR];
    const Array& jRRLv = jRRvv[pL];
    const Array& jRRRv = jRRvv[pR];

    // interpolate between the tables (in 3D)
    Array jv =    (1-ha)*(1-hZ)*(1-ht) * jLLLv
                + (1-ha)*(1-hZ)*(  ht) * jLLRv
                + (1-ha)*(  hZ)*(1-ht) * jLRLv
                + (1-ha)*(  hZ)*(  ht) * jLRRv
                + (  ha)*(1-hZ)*(1-ht) * jRLLv
                + (  ha)*(1-hZ)*(  ht) * jRLRv
                + (  ha)*(  hZ)*(1-ht) * jRRLv
                + (  ha)*(  hZ)*(  ht) * jRRRv;

    // resample to the possibly redshifted simulation wavelength grid,
    // convert emissivities to luminosities (i.e. multiply by the wavelength bins),
    // multiply by the mass of the population (in solar masses),
    // and return the result
    return NR::resample<NR::interpolateLogLog>(_lambdagrid->lambdav()*(1-z), _lambdav, jv)
                                    * _lambdagrid->dlambdav() * M;
}

//////////////////////////////////////////////////////////////////////

int FSPSVariableIMFSEDFamily::numParams() const
{
    return 4;
}

//////////////////////////////////////////////////////////////////////

Array FSPSVariableIMFSEDFamily::luminositiesGeneric(const Array& params, int skipvals, double z) const
{
    return luminosities(params[skipvals], params[skipvals+1], params[skipvals+2], params[skipvals+3], z);
}

//////////////////////////////////////////////////////////////////////

double FSPSVariableIMFSEDFamily::massGeneric(const Array& params, int skipvals) const
{
    return params[skipvals];
}

//////////////////////////////////////////////////////////////////////

string FSPSVariableIMFSEDFamily::sourceName() const
{
    return "star";
}

//////////////////////////////////////////////////////////////////////

string FSPSVariableIMFSEDFamily::sourceDescription() const
{
    return "star";
}

//////////////////////////////////////////////////////////////////////

const ArrayTable<2>& FSPSVariableIMFSEDFamily::getEmissivityTable(size_t s, size_t m) const
{
    // if the requested table is in the cache, simply return it
    size_t sm = key(s,m);
    if (_jmap.count(sm)) return _jmap.at(sm);

    // otherwise, create a new table and insert it with the appropriate key
    ArrayTable<2>& jvv = _jmap[sm];

    // allocate appropriate table dimensions (indices p,k)
    jvv.resize(Nt, Nlambda);

    // local constant for units
    const double Lsun = Constants::Lsun();
    const double c = Constants::c();

    // open the appropriate resource file
    string imf = imfString(initialMassFunction());
    string slope = StringUtils::toString(s, 'd', 0, 2, '0');
    string metal = StringUtils::toString(m, 'd', 0, 2, '0');
    string filename = FilePaths::externalResource("SSP_" + imf + "_imf00" + slope + "_z" + metal + ".out.spec");
    std::ifstream file = System::ifstream(filename);
    if (!file.is_open()) throw FATALERROR("Could not open the external resource file " + filename);
    find<Log>()->info("Reading variable IMF SED template from " + filename);

    // skip header comments and two intro lines
    string line;
    while (file.peek() == '#') getline(file,line);
    getline(file,line);
    getline(file,line);

    // for each age, skip the intro line and read the line with emissivities for each wavelength
    for (size_t p=0; p<Nt; ++p)
    {
        Array& jv = jvv[p];

        getline(file,line); // skip intro line
        for (size_t k=0; k<Nlambda; ++k)
        {
            double jnu;
            file >> jnu;
            jv[k] = jnu * Lsun * c/(_lambdav[k]*_lambdav[k]);   // emissivity in file in Lsun/Hz, we want in W/m
        }
        getline(file,line); // skip end of wavelength line
    }

    return jvv;
}

//////////////////////////////////////////////////////////////////////
