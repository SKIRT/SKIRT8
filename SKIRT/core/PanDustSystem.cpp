/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PanDustSystem.hpp"
#include "ArrayTable.hpp"
#include "DustGrid.hpp"
#include "DustMix.hpp"
#include "FatalError.hpp"
#include "FITSInOut.hpp"
#include "ISRF.hpp"
#include "LockFree.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "StaggeredAssigner.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void PanDustSystem::setupSelfBefore()
{
    DustSystem::setupSelfBefore();

    // if there is no dust emission, turn off the irrelevant flags
    if (!_dustEmissivity)
    {
        _dustLib = nullptr;
        _includeSelfAbsorption = false;
        _writeEmissivity = false;
        _writeTemperature = false;
        _writeISRF = false;
    }
    // if there is dust emission, make sure that there is a dust library as well
    else
    {
        if (!_dustLib) throw FATALERROR("There should be a dust library when dust emission is turned on");
    }

    // verify that the wavelength range includes the V-band center 0.55 micron (needed for normalization of dust)
    WavelengthGrid* lambdagrid = find<WavelengthGrid>();
    if (lambdagrid->nearest(0.55e-6) < 0)
        throw FATALERROR("Wavelength range should include 0.55 micron for a panchromatic simulation with dust");

    // cache size of wavelength grid
    _Nlambda = lambdagrid->numWavelengths();
}

////////////////////////////////////////////////////////////////////

namespace
{
    void writeEmissivitiesForField(PanDustSystem* ds, const Array& Jv, string filebody, string title)
    {
        WavelengthGrid* lambdagrid = ds->find<WavelengthGrid>();
        Units* units = ds->find<Units>();

        // Create a text file
        TextOutFile file(ds, filebody, "emissivities for " + title);

        // Get the emissivity for each dust mix
        int Ncomp = ds->numComponents();
        ArrayTable<2> evv(Ncomp,0);
        for (int h=0; h<Ncomp; h++) evv(h) = ds->dustEmissivity()->emissivity(ds->mix(h), Jv);

        // Write the header
        file.writeLine("# Dust emissivities for " + title);
        file.addColumn("lambda (" + units->uwavelength() + ")");
        file.addColumn("embedding field mean intensity -- J_lambda (W/m3/sr)");
        for (int h=0; h<Ncomp; h++)
            file.addColumn("dust mix " + std::to_string(h) + " -- lambda*j_lambda (W/sr/H)");

        // Write the input field and the emissivity for each dust mix to file
        int Nlambda = lambdagrid->numWavelengths();
        for (int ell=0; ell<Nlambda; ell++)
        {
            double lambda = lambdagrid->lambda(ell);
            vector<double> values;
            values.push_back(units->owavelength(lambda));
            values.push_back(Jv[ell]);
            for (int h=0; h<Ncomp; h++)
                values.push_back(ds->mix(h)->mu()*lambda*evv(h,ell));
            file.writeRow(values);
        }
    }
}

////////////////////////////////////////////////////////////////////

void PanDustSystem::setupSelfAfter()
{
    DustSystem::setupSelfAfter();
    int Ncells = dustGrid()->numCells();

    PeerToPeerCommunicator* comm = find<PeerToPeerCommunicator>();
    bool dataParallel = comm->dataParallel();

    if (dataParallel) // assign this process to work with a subset of dust cells
        _assigner = new StaggeredAssigner(Ncells, this);

    WavelengthGrid* wg = find<WavelengthGrid>();

    // resize the tables that hold the absorbed energies for each dust cell and wavelength
    // - absorbed stellar emission is relevant for calculating dust emission
    // - absorbed dust emission is relevant for calculating dust self-absorption
    _haveLabsStel = false;
    _haveLabsDust = false;
    if (hasDustEmission())
    {
        string tableName = "Absorbed stellar luminosity";
        if (dataParallel)
            _LabsStelvv.initialize(tableName, ParallelTable::WriteState::COLUMN, wg->assigner(), _assigner, comm);
        else
            _LabsStelvv.initialize(tableName, ParallelTable::WriteState::COLUMN, _Nlambda, Ncells, comm);
        _haveLabsStel = true;

        if (includeSelfAbsorption())
        {
            string tableName = "Absorbed dust luminosity";
            if (dataParallel)
                _LabsDustvv.initialize(tableName, ParallelTable::WriteState::COLUMN, wg->assigner(), _assigner, comm);
            else
                _LabsDustvv.initialize(tableName, ParallelTable::WriteState::COLUMN, _Nlambda, Ncells, comm);
            _haveLabsDust = true;
        }
    }

    // write emissivities if so requested
    if (writeEmissivity())
    {
        // write emissivities for a range of scaled Mathis ISRF input fields
        Array Jv = ISRF::mathis(this);
        for (int i=-4; i<7; i++)
        {
            double U = pow(10.,i);
            writeEmissivitiesForField(this, U*Jv, "Mathis_U_" + StringUtils::toString(U,'e',0),
                                      StringUtils::toString(U, 'g') + " * Mathis ISRF");
        }

        // write emissivities for a range of diluted Black Body input fields
        const int Tv[] = { 3000, 6000, 9000, 12000, 15000, 18000 };
        const double Dv[] = { 8.28e-12, 2.23e-13, 2.99e-14, 7.23e-15, 2.36e-15, 9.42e-16 };
        for (int i=0; i<6; i++)
        {
            Jv = Dv[i] * ISRF::blackbody(this, Tv[i]);
            writeEmissivitiesForField(this, Jv,
                "BlackBody_T_" + StringUtils::toString(Tv[i],'d',0,5,'0'),
                StringUtils::toString(Dv[i],'e',2) + " * B(" + StringUtils::toString(Tv[i]) + "K)" );
        }

        find<Log>()->info("Done writing emissivities.");
    }
}

////////////////////////////////////////////////////////////////////

void PanDustSystem::setEmulationMode()
{
    _emulationMode = true;
    // In emulation mode, always perform exactly one self-absorption iteration
    if (_dustEmissivity && _includeSelfAbsorption) _minIterations = _maxIterations = 1;
}

////////////////////////////////////////////////////////////////////

bool PanDustSystem::emulationMode()
{
    return _emulationMode;
}

////////////////////////////////////////////////////////////////////

const ProcessAssigner* PanDustSystem::assigner() const
{
    return _assigner;
}

//////////////////////////////////////////////////////////////////////

bool PanDustSystem::hasDustEmission() const
{
    return _dustEmissivity!=0;
}

////////////////////////////////////////////////////////////////////

bool PanDustSystem::hasDustAbsorption() const
{
    return hasDustEmission();
}

//////////////////////////////////////////////////////////////////////

void PanDustSystem::absorb(int m, int ell, double DeltaL, bool ynstellar)
{
    if (ynstellar)
    {
        if (!_haveLabsStel) throw FATALERROR("This dust system does not support absorption of stellar emission");
        LockFree::add(_LabsStelvv(m,ell), DeltaL);
    }
    else
    {
        if (!_haveLabsDust) throw FATALERROR("This dust system does not support absorption of dust emission");
        LockFree::add(_LabsDustvv(m,ell), DeltaL);
    }
}

//////////////////////////////////////////////////////////////////////

void PanDustSystem::resetDustAbsorption()
{
    _LabsDustvv.reset();
}

//////////////////////////////////////////////////////////////////////

double PanDustSystem::absorbedLuminosity(int m, int ell) const
{
    // Only callable on cells assigned to this process, and after sumResults
    double sum = 0;
    if (_haveLabsStel) sum += _LabsStelvv(m,ell);
    if (_haveLabsDust) sum += _LabsDustvv(m,ell);
    return sum;
}

//////////////////////////////////////////////////////////////////////

double PanDustSystem::absorbedLuminosity(int m) const
{
    // Only callable on cells assigned to this process, and after sumResults
    double sum = 0;

    if (_haveLabsStel)
        sum += _LabsStelvv.sumRow(m);
    if (_haveLabsDust)
        sum += _LabsDustvv.sumRow(m);

    return sum;
}

//////////////////////////////////////////////////////////////////////

Array PanDustSystem::absorbedLuminosity() const
{
    // Only callable on cells assigned to this process, and after sumResults
    Array sum(dustGrid()->numCells());

    if (_haveLabsStel)
        sum += _LabsStelvv.stackColumns();
    if (_haveLabsDust)
        sum += _LabsDustvv.stackColumns();

    return sum;
}

//////////////////////////////////////////////////////////////////////

double PanDustSystem::absorbedStellarLuminosity() const
{
    return _LabsStelvv.sumEverything();
}

//////////////////////////////////////////////////////////////////////

double PanDustSystem::absorbedDustLuminosity() const
{
    return _LabsDustvv.sumEverything();
}

////////////////////////////////////////////////////////////////////

void PanDustSystem::calculateDustEmission()
{
    if (_dustEmissivity)
    {
        sumResults();
        _dustLib->calculate();
    }
}

////////////////////////////////////////////////////////////////////

void PanDustSystem::sumResults()
{
    if (_haveLabsStel) _LabsStelvv.switchScheme();
    if (_haveLabsDust) _LabsDustvv.switchScheme();
}

////////////////////////////////////////////////////////////////////

double PanDustSystem::emittedDustLuminosity(int m, int ell) const
{
    // Only callable on wavelengths assigned to this process.
    return _dustEmissivity ? _dustLib->luminosity(m,ell) : 0.;
}

////////////////////////////////////////////////////////////////////

// Private class to output a FITS file with the mean dust temperatures
// in one of the coordinate planes (xy, xz, or yz).
namespace
{
    // The image size in each direction, in pixels
    const int Np = 1024;

    class WriteTempCut : public ParallelTarget
    {
    private:
        // cached values initialized in constructor
        const PanDustSystem* _ds;
        bool _dataParallel;
        const ProcessAssigner* _cellAssigner;
        DustGrid* _grid;
        Units* _units;
        Log* _log;
        PeerToPeerCommunicator* _comm;
        double xbase, ybase, zbase, xpsize, ypsize, zpsize, xcenter, ycenter, zcenter;
        int Nmaps;

        // data members initialized in setup()
        bool xd, yd, zd;  // direction of coordinate plane (110, 101, 011)
        string plane;    // name of the coordinate plane (xy, xz, yz)

        // results vector, properly sized in constructor and initialized to zero in setup()
        Array tempv;

    public:
        // constructor
        WriteTempCut(const PanDustSystem* ds)
        {
            _ds = ds;
            _cellAssigner = ds->assigner();
            _grid = ds->dustGrid();
            _units = ds->find<Units>();
            _log = ds->find<Log>();
            _comm = ds->find<PeerToPeerCommunicator>();
            _dataParallel = _comm->dataParallel();

            double xmin, ymin, zmin, xmax, ymax, zmax;
            _grid->boundingBox().extent(xmin,ymin,zmin,xmax,ymax,zmax);
            xpsize = (xmax-xmin)/Np;
            ypsize = (ymax-ymin)/Np;
            zpsize = (zmax-zmin)/Np;
            xbase = xmin + 0.5*xpsize;
            ybase = ymin + 0.5*ypsize;
            zbase = zmin + 0.5*zpsize;
            xcenter = (xmin+xmax)/2.0;
            ycenter = (ymin+ymax)/2.0;
            zcenter = (zmin+zmax)/2.0;

            Nmaps = 0;
            for (int h=0; h<_ds->numComponents(); h++) Nmaps += _ds->mix(h)->numPopulations();

            tempv.resize(Np*Np*Nmaps);
        }

        // setup for calculating a specific coordinate plane
        void setup(bool xdir, bool ydir, bool zdir)
        {
            xd = xdir;
            yd = ydir;
            zd = zdir;
            plane = "";
            if (xd) plane += "x";
            if (yd) plane += "y";
            if (zd) plane += "z";
            _log->info("Calculating dust temperatures in the " + plane + " plane...");

            tempv = 0.0;  // initialize all values to zero to facilitate the code in body()
        }

        // the parallized loop body; calculates the results for a single line in the images
        void body(size_t j)
        {
            double z = zd ? (zbase + j*zpsize) : 0.;
            for (int i=0; i<Np; i++)
            {
                double x = xd ? (xbase + i*xpsize) : 0.;
                double y = yd ? (ybase + (zd ? i : j)*ypsize) : 0.;
                Position bfr(x,y,z);
                int m = _grid->whichCell(bfr);

                bool m_is_available = !_dataParallel || _cellAssigner->validIndex(m);
                if (m_is_available && m!=-1 && _ds->absorbedLuminosity(m)>0.0)
                {
                    const Array& Jv = _ds->meanIntensity(m);
                    int p = 0;
                    for (int h=0; h<_ds->numComponents(); h++)
                    {
                        double rho = _ds->density(m,h);
                        int Npop = _ds->mix(h)->numPopulations();
                        for (int c=0; c<Npop; c++)
                        {
                            if (rho>0.0)
                            {
                                double T = _ds->mix(h)->equilibrium(Jv,c);
                                int l = i + Np*j + Np*Np*p;
                                tempv[l] = _units->otemperature(T);
                            }
                            p++;
                        }
                    }
                }
            }
        }

        // Write the results to a FITS file with an appropriate name
        void write()
        {
            // If we didn't have all the cells, sum the results first
            if (_dataParallel) _comm->sum(tempv);

            string filename = "ds_temp" + plane;
            string description = "dust temperatures";
            FITSInOut::write(_ds, description, filename, tempv, Np, Np, Nmaps,
                             _units->olength(xd?xpsize:ypsize), _units->olength(zd?zpsize:ypsize),
                             _units->olength(xd?xcenter:ycenter), _units->olength(zd?zcenter:ycenter),
                             _units->utemperature(), _units->ulength());
        }
    };
}

////////////////////////////////////////////////////////////////////

// Private class to output a text file with an indicative temperature for each dust cell
namespace
{
    class WriteTempData : public ParallelTarget
    {
    private:
        // cached values initialized in constructor
        const PanDustSystem* _ds;
        DustGrid* _grid;
        Units* _units;
        int _Ncells;

        // results vectors, properly sized in constructor
        Array _Mv, _Tv;

    public:
        // constructor
        WriteTempData(const PanDustSystem* ds)
        {
            _ds = ds;
            _grid = ds->dustGrid();
            _units = ds->find<Units>();
            _Ncells = ds->numCells();
            _Mv.resize(_Ncells);
            _Tv.resize(_Ncells);
        }

        // the parallized loop body; calculates the results for a single dust cell
        void body(size_t m)
        {
            // dust mass in cell
            _Mv[m] = _ds->density(m) * _ds->volume(m);

            // indicative temperature = average population equilibrium temperature weighed by population mass fraction
            if (_ds->absorbedLuminosity(m)>0.0)
            {
                const Array& Jv = _ds->meanIntensity(m);

                // average over dust components
                double sumRho_h = 0;
                double sumRhoT_h = 0;
                for (int h=0; h<_ds->numComponents(); h++)
                {
                    double rho_h = _ds->density(m,h);
                    if (rho_h>0.0)
                    {
                        // average over dust populations within component
                        double sumMu_c = 0;
                        double sumMuT_c = 0;
                        for (int c=0; c<_ds->mix(h)->numPopulations(); c++)
                        {
                            double mu_c = _ds->mix(h)->mu(c);
                            double T_c = _ds->mix(h)->equilibrium(Jv,c);
                            sumMu_c += mu_c;
                            sumMuT_c += mu_c * T_c;
                        }
                        double T_h = sumMuT_c / sumMu_c;

                        sumRho_h += rho_h;
                        sumRhoT_h += rho_h * T_h;
                    }
                }
                _Tv[m] = sumRhoT_h / sumRho_h;
            }
        }

        // Write the results to a text file with an appropriate name
        void write()
        {
            PeerToPeerCommunicator* comm = _ds->find<PeerToPeerCommunicator>();
            // Sum the calculated results if necessary
            if (comm->dataParallel())
            {
                comm->sum(_Tv);
                comm->sum(_Mv);
            }

            // Create a text file
            TextOutFile file(_ds, "ds_celltemps", "dust cell temperatures");

            // Write the header
            file.addColumn("dust cell index", 'd');
            file.addColumn("indicative temperature in cell (" + _units->utemperature() + ")", 'g');

            // Write a line for each cell
            for (int m=0; m<_Ncells; m++)
            {
                file.writeRow(vector<double>({ static_cast<double>(m), _units->otemperature(_Tv[m]) }));
            }
        }
    };
}

////////////////////////////////////////////////////////////////////

void PanDustSystem::write() const
{
    DustSystem::write();

    PeerToPeerCommunicator* comm = find<PeerToPeerCommunicator>();
    bool dataParallel = comm->dataParallel();

    // If requested, output temperature map(s) along coordinate axes and temperature data for each dust cell
    if (_writeTemperature)
    {
        // Parallelize the calculation over the threads
        Parallel* parallel = find<ParallelFactory>()->parallel();
        // If the necessary data is distributed over the processes, do the calculation on all processes.
        // Else, let the root do everything.
        bool isRoot = comm->isRoot();

        // Output temperature map(s) along coordinate axes
        {
            // Construct a private class instance to do the work (parallelized)
            WriteTempCut wt(this);

            // Get the dimension of the dust grid
            int dimDust = dustGrid()->dimension();

            // For the xy plane (always)
            {
                wt.setup(1,1,0);
                if (dataParallel) parallel->call(&wt, Np);
                else if (isRoot) parallel->call(&wt, Np);
                wt.write();
            }

            // For the xz plane (only if dimension is at least 2)
            if (dimDust >= 2)
            {
                wt.setup(1,0,1);
                if (dataParallel) parallel->call(&wt, Np);
                else if (isRoot) parallel->call(&wt, Np);
                wt.write();
            }

            // For the yz plane (only if dimension is 3)
            if (dimDust == 3)
            {
                wt.setup(0,1,1);
                if (dataParallel) parallel->call(&wt, Np);
                else if (isRoot) parallel->call(&wt, Np);
                wt.write();
            }
        }

        // Output a text file with temperature data for each dust cell
        {
            find<Log>()->info("Calculating indicative dust temperatures for each cell...");

            // Construct a private class instance to do the work (parallelized)
            WriteTempData wt(this);

            // Call the body on the right cells. If everything is available, no unnecessary communication will be done.
            if (dataParallel)
            {
                // Calculate the temperature for the cells owned by this process
                parallel->call(&wt, _assigner);
            }
            else if (isRoot)
            {
                // Let root calculate it for everything
                parallel->call(&wt, dustGrid()->numCells());
            }
            wt.write();
        }
    }

    // If requested, output the interstellar radiation field in every dust cell to a data file
    if (_writeISRF)
    {
        WavelengthGrid* lambdagrid = find<WavelengthGrid>();
        Units* units = find<Units>();

        // Create a text file
        TextOutFile file(this, "ds_isrf", "ISRF");

        // Write the header
        file.writeLine("# Mean field intensities for all dust cells");
        file.addColumn("dust cell index", 'd');
        file.addColumn("bolometric luminosity absorbed in cell (" + units->ubolluminosity() + ")", 'g');
        for (int ell=0; ell<_Nlambda; ell++)
            file.addColumn("J_lambda (W/m3/sr) for lambda = "
                           + StringUtils::toString(units->owavelength(lambdagrid->lambda(ell)), 'g')
                           + " " + units->uwavelength(), 'g');

        // Write one line for each dust cell
        int Ncells = dustGrid()->numCells();
        for (int m=0; m<Ncells; m++)
        {
            if (!dataParallel)
            {
                vector<double> values({ static_cast<double>(m), units->obolluminosity(absorbedLuminosity(m)) });
                for (auto J : meanIntensity(m)) values.push_back(J);
                file.writeRow(values);
            }
            // for distributed mode
            else
            {
                // the correct process gets Labs and Jv and broadcasts them
                Array Labs(1);
                Array Jv(_Nlambda);
                if (_assigner->validIndex(m))
                {
                    Labs[0] = absorbedLuminosity(m);
                    Jv = meanIntensity(m);
                }
                int sender = _assigner->rankForIndex(m);
                comm->broadcast(Labs, sender);
                comm->broadcast(Jv, sender);

                // write the line
                vector<double> values({ static_cast<double>(m), units->obolluminosity(Labs[0]) });
                for (auto J : Jv) values.push_back(J);
                file.writeRow(values);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
