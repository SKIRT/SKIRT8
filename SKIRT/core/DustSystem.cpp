/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustSystem.hpp"
#include "Constants.hpp"
#include "DustGridDensityInterface.hpp"
#include "DustGridPath.hpp"
#include "DustMix.hpp"
#include "DustSystemDensityCalculator.hpp"
#include "DustSystemDepthCalculator.hpp"
#include "FatalError.hpp"
#include "FITSInOut.hpp"
#include "Geometry.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "PhotonPackage.hpp"
#include "Random.hpp"
#include "ShortArray.hpp"
#include "StaggeredAssigner.hpp"
#include "StellarSystem.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "TimeLogger.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

void DustSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // Copy some basic properties
    _Ncomp = _dd->numComponents();
    _Ncells = _grid->numCells();

    // Make sure that all dust mixes support polarization, or none of them do
    for (int h=1; h<_Ncomp; h++)
    {
        if (_dd->mix(h)->polarization() != _dd->mix(0)->polarization())
            throw FATALERROR("All dust mixes must consistenly support polarization, or not support polarization");
    }

    // Resize the tables that hold essential dust cell properties
    _volumev.resize(_Ncells);
    _rhovv.resize(_Ncells,_Ncomp);

    // Set the volume of the cells (parallelized over different threads, except when multiprocessing is enabled)
    find<Log>()->info("Calculating the volume of the cells...");
    PeerToPeerCommunicator* comm = find<PeerToPeerCommunicator>();
    ParallelFactory* pfactory = find<ParallelFactory>();
    size_t nthreads = comm->isMultiProc() ? 1 : pfactory->maxThreadCount();
    pfactory->parallel(nthreads)->call(this, &DustSystem::setVolumeBody, _Ncells);

    // use a StaggeredAssigner to calculate the densities
    _setupAssigner = new StaggeredAssigner(_Ncells, this);

    // Calculate and set the density of the cells that are assigned to this process
    _gdi = _grid->interface<DustGridDensityInterface>();
    if (_gdi)
    {
        // if the dust grid offers a special interface, use it
        find<Log>()->info("Setting the value of the density in the cells using grid interface...");
        find<ParallelFactory>()->parallel()->call(this, &DustSystem::setGridDensityBody, _setupAssigner);
    }
    else
    {
        // otherwise take an average of the density in 100 random positions in the cell (parallelized)
        find<Log>()->info("Setting the value of the density in the cells...");
        find<ParallelFactory>()->parallel()->call(this, &DustSystem::setSampleDensityBody, _setupAssigner);
    }

    // Wait for the other processes to reach this point
    comm->wait("the calculation of the dust cell densities");

    // Obtain the densities in all dust cells, if the calculation has been performed by parallel processes
    if (comm->isMultiProc()) assemble();

    // Perform a convergence check on the grid.
    if (_writeConvergence) doWriteConvergence();

    // Write the density in the xy plane, xz plane and yz plane to a file.
    if (_writeDensity) doWriteDensity();

    // Output optical depth map as seen from the center
    if (_writeDepthMap) doWriteDepthMap();

    // Calculate and output some quality metrics for the dust grid
    if (_writeQuality) doWriteQuality();

    // Output properties for all cells in the dust grid
    if (_writeCellProperties) doWriteCellProperties();

    // Output stellar distribution on the dust grid
    if (_writeStellarDensity) doWriteStellarDensity();
}

////////////////////////////////////////////////////////////////////

// parallelized body used above
void DustSystem::setVolumeBody(size_t m)
{
    _volumev[m] = (_grid->weight(m) > 0) ? _grid->volume(m) : 0;
}

////////////////////////////////////////////////////////////////////

// parallelized body used above
void DustSystem::setGridDensityBody(size_t m)
{
    for (int h=0; h<_Ncomp; h++)
        _rhovv(m,h) = _gdi->density(h,m);
}

////////////////////////////////////////////////////////////////////

// parallelized body used above
void DustSystem::setSampleDensityBody(size_t m)
{
    if (m%100000==0)
    {
        find<Log>()->info("  Computing density for cell " + std::to_string(m)
                          + " ("
                          + std::to_string(100*_setupAssigner->relativeIndex(m)/_setupAssigner->assigned())
                          + "%)");
    }
    double weight = _grid->weight(m);
    if (weight > 0)
    {
        Array sumv(_Ncomp);
        for (int n=0; n<_numSamples; n++)
        {
            Position bfr = _grid->randomPositionInCell(m);
            for (int h=0; h<_Ncomp; h++) sumv[h] += _dd->density(h,bfr);
        }
        for (int h=0; h<_Ncomp; h++)
        {
            _rhovv(m,h) = weight*sumv[h]/_numSamples;
        }
    }
    else
    {
        for (int h=0; h<_Ncomp; h++) _rhovv(m,h) = 0;
    }
}

////////////////////////////////////////////////////////////////////

void DustSystem::assemble()
{
    // Get a pointer to the PeerToPeerCommunicator of this simulation
    PeerToPeerCommunicator* comm = find<PeerToPeerCommunicator>();

    Log* log = find<Log>();
    TimeLogger logger(log->verbose() && comm->isMultiProc() ? log : 0, "communication of the dust densities");

    // Sum the densities array across all processes
    comm->sumAll(_rhovv.data());
}

////////////////////////////////////////////////////////////////////

void DustSystem::doWriteConvergence() const
{
    // Perform a convergence check on the grid. First calculate the total
    // mass and the principle axes surface densities by integrating
    // over the grid.

    find<Log>()->info("Performing a convergence check on the grid...");

    // calculation of the mass

    double M = 0.0;
    for (int m=0; m<_Ncells; m++)
    {
        M += density(m)*volume(m);
    }

    // calculation of the X-axis surface density

    double SigmaX = 0;
    DustGridPath dgp(Position(0.,0.,0.), Direction(1.,0.,0.));
    _grid->path(&dgp);
    SigmaX += dgp.opticalDepth([this](int m){ return density(m); });
    dgp.setDirection(Direction(-1.,0.,0.));
    _grid->path(&dgp);
    SigmaX += dgp.opticalDepth([this](int m){ return density(m); });

    // calculation of the Y-axis surface density

    double SigmaY = 0.0;
    dgp.setDirection(Direction(0.,1.,0.));
    _grid->path(&dgp);
    SigmaY += dgp.opticalDepth([this](int m){ return density(m); });
    dgp.setDirection(Direction(0.,-1.,0.));
    _grid->path(&dgp);
    SigmaY += dgp.opticalDepth([this](int m){ return density(m); });

    // calculation of the Z-axis surface density

    double SigmaZ = 0.0;
    dgp.setDirection(Direction(0.,0.,1.));
    _grid->path(&dgp);
    SigmaZ += dgp.opticalDepth([this](int m){ return density(m); });
    dgp.setDirection(Direction(0.,0.,-1.));
    _grid->path(&dgp);
    SigmaZ += dgp.opticalDepth([this](int m){ return density(m); });

    // Compare these values to the expected values and write the result to file

    Units* units = find<Units>();

    // Create a text file
    TextOutFile file(this, "ds_convergence", "convergence check on the dust system");
    file.writeLine("Convergence check on the grid: ");

    switch (_dd->dimension())
    {
        case 1:
        {
            double Sigmar = 0.5*SigmaX;
            double Sigmarref = 0.5*_dd->SigmaX();

            // Write the surface densities to file
            file.writeLine("   - radial (r-axis) surface density");
            file.writeLine("         expected value = " + StringUtils::toString(units->omasssurfacedensity(Sigmarref), 'g')
                                                        + " " + units->umasssurfacedensity());
            file.writeLine("         actual value =   " + StringUtils::toString(units->omasssurfacedensity(Sigmar), 'g')
                                                        + " " + units->umasssurfacedensity());
        }
        break;
        case 2:
        {
            double SigmaR = 0.5*SigmaX;
            double SigmaRref = 0.5*_dd->SigmaX();
            double SigmaZref = _dd->SigmaZ();

            // Write the surface densities to file
            file.writeLine("   - edge-on (R-axis) surface density");
            file.writeLine("         expected value = " + StringUtils::toString(units->omasssurfacedensity(SigmaRref), 'g')
                                                        + " " + units->umasssurfacedensity());
            file.writeLine("         actual value =   " + StringUtils::toString(units->omasssurfacedensity(SigmaR), 'g')
                                                        + " " + units->umasssurfacedensity());
            file.writeLine("   - face-on (Z-axis) surface density");
            file.writeLine("         expected value = " + StringUtils::toString(units->omasssurfacedensity(SigmaZref), 'g')
                                                        + " " + units->umasssurfacedensity());
            file.writeLine("         actual value =   " + StringUtils::toString(units->omasssurfacedensity(SigmaZ), 'g')
                                                        + " " + units->umasssurfacedensity());
        }
        break;
        case 3:
        {
            double SigmaXref = _dd->SigmaX();
            double SigmaYref = _dd->SigmaY();
            double SigmaZref = _dd->SigmaZ();

            // Write the surface densities to file
            file.writeLine("   - X-axis surface density");
            file.writeLine("         expected value = " + StringUtils::toString(units->omasssurfacedensity(SigmaXref), 'g')
                                                        + " " + units->umasssurfacedensity());
            file.writeLine("         actual value =   " + StringUtils::toString(units->omasssurfacedensity(SigmaX), 'g')
                                                        + " " + units->umasssurfacedensity());
            file.writeLine("   - Y-axis surface density");
            file.writeLine("         expected value = " + StringUtils::toString(units->omasssurfacedensity(SigmaYref), 'g')
                                                        + " " + units->umasssurfacedensity());
            file.writeLine("         actual value =   " + StringUtils::toString(units->omasssurfacedensity(SigmaY), 'g')
                                                        + " " + units->umasssurfacedensity());
            file.writeLine("   - Z-axis surface density");
            file.writeLine("         expected value = " + StringUtils::toString(units->omasssurfacedensity(SigmaZref), 'g')
                                                        + " " + units->umasssurfacedensity());
            file.writeLine("         actual value =   " + StringUtils::toString(units->omasssurfacedensity(SigmaZ), 'g')
                                                        + " " + units->umasssurfacedensity());
        }
        break;
        default:
            throw FATALERROR("Wrong dimension in dust distribution");
    }
    double Mref = _dd->mass();

    // Write the (expected and actual) total dust mass
    file.writeLine("   - total dust mass");
    file.writeLine("         expected value = " + StringUtils::toString(units->omass(Mref), 'g') + " " + units->umass());
    file.writeLine("         actual value =   " + StringUtils::toString(units->omass(M), 'g') + " " + units->umass());
}

////////////////////////////////////////////////////////////////////

// Private class to output FITS files with the theoretical dust density and the constructed grid density
// in one of the coordinate planes (xy, xz, or yz).
namespace
{
    // The image size in each direction, in pixels
    const int Np = 1024;

    class WriteDensity : public ParallelTarget
    {
    private:
        // results -- sized to fit in constructor
        Array trhov;
        Array grhov;

        // data members initialized in constructor
        const DustSystem* _ds;
        DustDistribution* _dd;
        DustGrid* _grid;
        Units* _units;
        Log* _log;
        double xbase, ybase, zbase, xpsize, ypsize, zpsize, xcenter, ycenter, zcenter;

        // data members initialized in setup()
        bool xd, yd, zd; // direction of coordinate plane (110, 101, 011)
        string plane;    // name of the coordinate plane (xy, xz, yz)

    public:
        // constructor
        WriteDensity(const DustSystem* ds)
            : trhov(Np*Np), grhov(Np*Np)
        {
            _ds = ds;
            _dd = ds->dustDistribution();
            _grid = ds->dustGrid();
            _units = ds->find<Units>();
            _log = ds->find<Log>();

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
            _log->info("Calculating density in the " + plane + " plane...");
        }

        // the parallized loop body; calculates the results for a single line in the images
        void body(size_t j)
        {
            double z = zd ? (zbase + j*zpsize) : 0.;
            for (int i=0; i<Np; i++)
            {
                int l = i + Np*j;
                double x = xd ? (xbase + i*xpsize) : 0.;
                double y = yd ? (ybase + (zd ? i : j)*ypsize) : 0.;
                Position bfr(x,y,z);
                trhov[l] = _units->omassvolumedensity(_dd->density(bfr));
                int m = _grid->whichCell(bfr);
                if (m==-1)
                    grhov[l] = 0.0;
                else
                    grhov[l] = _units->omassvolumedensity(_ds->density(m));
            }
        }

        // write the results to two FITS files with appropriate names
        void write()
        {
            write(trhov, "theoretical", "ds_trho");
            write(grhov, "grid", "ds_grho");
        }

    private:
        void write(const Array& rhov, string label, string prefix)
        {
            string filename = prefix + plane;
            FITSInOut::write(_ds, label, filename, rhov, Np, Np, 1,
                             _units->olength(xd?xpsize:ypsize), _units->olength(zd?zpsize:ypsize),
                             _units->olength(xd?xcenter:ycenter), _units->olength(zd?zcenter:ycenter),
                             _units->umassvolumedensity(), _units->ulength());
        }
    };
}

////////////////////////////////////////////////////////////////////

void DustSystem::doWriteDensity() const
{
    // construct a private class instance to do the work (parallelized)
    WriteDensity wd(this);
    Parallel* parallel = find<ParallelFactory>()->parallel();

    // get the dimension of the dust system
    int dimDust = dimension();

    // Only perform the calculation at the root processs
    bool isRoot = find<PeerToPeerCommunicator>()->isRoot();
    // For the xy plane (always)
    {
        wd.setup(1,1,0);
        if (isRoot) parallel->call(&wd, Np);
        wd.write();
    }

    // For the xz plane (only if dimension is at least 2)
    if (dimDust >= 2)
    {
        wd.setup(1,0,1);
        if (isRoot) parallel->call(&wd, Np);
        wd.write();
    }

    // For the yz plane (only if dimension is 3)
    if (dimDust == 3)
    {
        wd.setup(0,1,1);
        if (isRoot) parallel->call(&wd, Np);
        wd.write();
    }
}

////////////////////////////////////////////////////////////////////

// Private class to encapsulate the call-back function for calculating optical depths
namespace
{
    class KappaRho
    {
    private:
        // data members initialized in constructor
        const DustSystem* _ds;
        int _Ncomp;
        ShortArray<8> _kappaextv;

    public:
        // constructor
        // stores the extinction coefficients at the specified wavelength for all dust mixes
        KappaRho(const DustSystem* ds, int ell) : _ds(ds), _Ncomp(ds->numComponents()), _kappaextv(_Ncomp)
        {
            for (int h=0; h<_Ncomp; h++)
                _kappaextv[h] = _ds->mix(h)->kappaext(ell);
        }

        // call-back function
        // returns kappa*rho for the specified cell number (and for the wavelength-index bound in the constructor)
        double operator() (int m) const
        {
            double result = 0;
            for (int h=0; h<_Ncomp; h++)
                result += _kappaextv[h] * _ds->density(m,h);
            return result;
        }
    };
}

////////////////////////////////////////////////////////////////////

// Private class to output a FITS file with an optical depth map viewed from the center using Mollweide projection
namespace
{
    // The image size in each direction, in pixels
    const int Npx = 1600;
    const int Npy = 800;

    class WriteDepthMap : public ParallelTarget
    {
    private:
        // results -- sized to fit in constructor
        Array tauv;

        // data members initialized in constructor
        const DustSystem* _ds;
        DustGrid* _grid;
        Log* _log;
        int _ell;

    public:
        // constructor
        WriteDepthMap(const DustSystem* ds)
            : tauv(Npx*Npy)
        {
            _ds = ds;
            _grid = ds->dustGrid();
            _log = ds->find<Log>();
            _log->info("Calculating optical depth map viewed from the center...");
            _ell = max(0, ds->find<WavelengthGrid>()->nearest(Constants::lambdaV()));
        }

        // the parallized loop body; calculates the results for a single line in the image
        void body(size_t j)
        {
            double y = (j+0.5) / Npy;
            for (int i=0; i<Npx; i++)
            {
                double x = (i+0.5) / Npx;

                // perform the inverse Mollweide projection
                double alpha = asin(2*y-1);
                double theta = acos((2*alpha+sin(2*alpha))/M_PI);
                double phi = M_PI*(2*x-1)/cos(alpha);

                // if the deprojected direction is within range, compute the optical depth
                if (phi > -M_PI && phi < M_PI)
                    tauv[i+Npx*j] = opticaldepth(_ell, Position(), Direction(theta, phi));
            }
        }

        // write the results to a FITS file with an appropriate name
        void write()
        {
            Units* units = _ds->find<Units>();
            WavelengthGrid* lambdagrid = _ds->find<WavelengthGrid>();

            string filename = "ds_tau";
            string description = "optical depth map at λ = "
                                  + StringUtils::toString(units->owavelength(lambdagrid->lambda(_ell)))
                                  + " " + units->uwavelength();

            FITSInOut::write(_ds, description, filename, tauv, Npx, Npy, 1,
                             units->oposangle(2*M_PI/Npx), units->oposangle(M_PI/Npy), 0., 0.,
                             "", units->uposangle());
        }

    private:
        double opticaldepth(int ell, Position bfr, Direction bfk)
        {
            DustGridPath dgp(bfr, bfk);
            _grid->path(&dgp);
            return dgp.opticalDepth(KappaRho(_ds, ell));
        }
    };
}

////////////////////////////////////////////////////////////////////

void DustSystem::doWriteDepthMap() const
{
    // construct a private class instance to do the work (parallelized)
    WriteDepthMap wdm(this);

    // perform the calculation at the root process
    Parallel* parallel = find<ParallelFactory>()->parallel();
    bool isRoot = find<PeerToPeerCommunicator>()->isRoot();
    if (isRoot) parallel->call(&wdm, Npy);
    wdm.write();
}

////////////////////////////////////////////////////////////////////

void DustSystem::doWriteQuality() const
{
    Log* log = find<Log>();
    Units* units = find<Units>();
    Parallel* parallel = find<ParallelFactory>()->parallel();
    bool isRoot = find<PeerToPeerCommunicator>()->isRoot();

    // Density metric

    log->info("Calculating quality metric for the grid density...");
    DustSystemDensityCalculator calc1(this, _numSamples, _Ncells/5);
    if (isRoot) parallel->call(&calc1, _numSamples);

    log->info("  Mean value of density delta: "
              + StringUtils::toString(units->omassvolumedensity(calc1.meanDelta()))
              + units->umassvolumedensity());
    log->info("  Standard deviation of density delta: "
              + StringUtils::toString(units->omassvolumedensity(calc1.stdDevDelta()))
              + units->umassvolumedensity());

    // Optical depth metric

    log->info("Calculating quality metric for the optical depth in the grid...");
    DustSystemDepthCalculator calc2(this, _numSamples, _Ncells/50, _numSamples*10);
    if (isRoot) parallel->call(&calc2, _numSamples);

    log->info("  Mean value of optical depth delta: " + StringUtils::toString(calc2.meanDelta()));
    log->info("  Standard deviation of optical depth delta: " + StringUtils::toString(calc2.stdDevDelta()));

    // Create a text file
    TextOutFile file(this, "ds_quality", "quality metrics for the grid");

    // Write quality metrics
    file.writeLine("Mean value of density delta: "
                    + StringUtils::toString(units->omassvolumedensity(calc1.meanDelta()), 'g') + ' '
                    + units->umassvolumedensity());
    file.writeLine("Standard deviation of density delta: "
                    + StringUtils::toString(units->omassvolumedensity(calc1.stdDevDelta()), 'g') + ' '
                    + units->umassvolumedensity());
    file.writeLine("Mean value of optical depth delta: " + StringUtils::toString(calc2.meanDelta(), 'g'));
    file.writeLine("Standard deviation of optical depth delta: " + StringUtils::toString(calc2.stdDevDelta(), 'g'));
}

////////////////////////////////////////////////////////////////////

void DustSystem::doWriteCellProperties() const
{
    Log* log = find<Log>();
    Units* units = find<Units>();

    // Create a text file
    TextOutFile file(this, "ds_cellprops", "dust cell properties");

    // Write the header
    file.addColumn("dust cell index", 'd');
    file.addColumn("x coordinate of cell center (" + units->ulength() + ")");
    file.addColumn("y coordinate of cell center (" + units->ulength() + ")");
    file.addColumn("z coordinate of cell center (" + units->ulength() + ")");
    file.addColumn("cell volume (" + units->uvolume() + ")");
    file.addColumn("average dust density in cell (" + units->umassvolumedensity() + ")");
    file.addColumn("dust mass in cell (" + units->umass() + ")");
    file.addColumn("mass fraction, i.e. dust mass in cell over total dust mass (dimensionless fraction)");
    file.addColumn("V-band optical depth of cell diagonal (dimensionless)");

    // Write a line for each cell; remember the tau values so we can compute some statistics
    Array tauV(_Ncells);
    double totalmass = _dd->mass();
    for (int m=0; m<_Ncells; m++)
    {
        Position bfr = dustGrid()->centralPositionInCell(m);
        double V = volume(m);
        double rho = density(m);
        double M = V*rho;
        double delta = M/totalmass;
        double tau = Constants::kappaV()*rho*pow(V,1./3.);
        file.writeRow(vector<double>({ static_cast<double>(m),
                                       units->olength(bfr.x()), units->olength(bfr.y()), units->olength(bfr.z()),
                                       units->ovolume(V), units->omassvolumedensity(rho), units->omass(M),
                                       delta, tau }));
        tauV[m] = tau;
    }

    // Calculate some statistics on optical depth
    double tauavg = tauV.sum()/_Ncells;
    double taumin = tauV.min();
    double taumax = tauV.max();
    const int Nbins = 500;
    vector<int> countV(Nbins+1);
    for (int m=0; m<_Ncells; m++)
    {
        int index = max(0,min(Nbins, static_cast<int>((tauV[m]-taumin)/(taumax-taumin)*Nbins)));
        countV[index]++;
    }
    int count = 0;
    int index = 0;
    for (; index<Nbins; index++)
    {
        count += countV[index];
        if (count > 0.9*_Ncells) break;
    }
    double tau90 = taumin + index*(taumax-taumin)/Nbins;

    // write the statistics on optical depth to the file
    file.writeLine("# smallest optical depth: " + StringUtils::toString(taumin, 'g'));
    file.writeLine("# largest optical depth:  " + StringUtils::toString(taumax, 'g'));
    file.writeLine("# average optical depth:  " + StringUtils::toString(tauavg, 'g'));
    file.writeLine("# 90 % of the cells have optical depth smaller than: " + StringUtils::toString(tau90, 'g'));

    // report the statistics on optical depth to the console
    log->info("  Smallest optical depth: " + StringUtils::toString(taumin));
    log->info("  Largest optical depth:  " + StringUtils::toString(taumax));
    log->info("  Average optical depth:  " + StringUtils::toString(tauavg));
    log->info("  90 % of the cells have optical depth smaller than: " + StringUtils::toString(tau90));
}

////////////////////////////////////////////////////////////////////

void DustSystem::doWriteStellarDensity() const
{
    // Get the geometry for the single stellar component of type GeometricStellarComp
    StellarSystem* stellarSystem = find<StellarSystem>();
    if (stellarSystem->numComponents() != 1) throw FATALERROR("There must be a single stellar component");
    Geometry* stellarGeometry = stellarSystem->components()[0]->find<Geometry>();

    // Create a text file
    TextOutFile file(this, "ds_stellar", "stellar density");

    // Write the header
    Units* units = find<Units>();
    file.addColumn("dust cell index", 'd');
    file.addColumn("stellar density (1/" + units->uvolume() + ")");

    // Write a line for each cell
    for (int m=0; m<_Ncells; m++)
    {
        double rho = stellarGeometry->density(_grid->centralPositionInCell(m));
        file.writeRow(vector<double>({ static_cast<double>(m), 1./units->ovolume(1./rho) }));
    }
}

////////////////////////////////////////////////////////////////////

int DustSystem::dimension() const
{
    return _dd->dimension();
}

//////////////////////////////////////////////////////////////////////

int DustSystem::numCells() const
{
    return _Ncells;
}

//////////////////////////////////////////////////////////////////////

int DustSystem::numComponents() const
{
    return _Ncomp;
}

////////////////////////////////////////////////////////////////////

bool DustSystem::polarization() const
{
    return _Ncomp>0 ? _dd->mix(0)->polarization() : false;
}

//////////////////////////////////////////////////////////////////////

DustMix* DustSystem::mix(int h) const
{
    return _dd->mix(h);
}

//////////////////////////////////////////////////////////////////////

DustMix* DustSystem::randomMixForPosition(Position bfr, int ell) const
{
    int hmix = 0;
    if (_Ncomp>1)
    {
        int m = whichCell(bfr);
        if (m>=0)
        {
            Array Xv;
            NR::cdf(Xv, _Ncomp, [this,ell,m](int h){return mix(h)->kappasca(ell)*density(m,h);} );
            hmix = NR::locateClip(Xv, find<Random>()->uniform());
        }
    }
    return mix(hmix);
}

//////////////////////////////////////////////////////////////////////

int DustSystem::whichCell(Position bfr) const
{
    return _grid->whichCell(bfr);
}

//////////////////////////////////////////////////////////////////////

Position DustSystem::randomPositionInCell(int m) const
{
    return _grid->randomPositionInCell(m);
}

//////////////////////////////////////////////////////////////////////

double DustSystem::volume(int m) const
{
    return _volumev[m];
}

//////////////////////////////////////////////////////////////////////

double DustSystem::density(int m, int h) const
{
    return m >= 0 ? _rhovv(m,h) : 0;
}

//////////////////////////////////////////////////////////////////////

double DustSystem::density(int m) const
{
    double rho = 0;
    if (m >= 0)
        for (int h=0; h<_Ncomp; h++) rho += _rhovv(m,h);
    return rho;
}

//////////////////////////////////////////////////////////////////////

Array DustSystem::meanIntensity(int m) const
{
    WavelengthGrid* lambdagrid = find<WavelengthGrid>();
    int Nlambda = lambdagrid->numWavelengths();
    Array Jv(Nlambda);
    double fac = 4.0*M_PI*volume(m);
    for (int ell=0; ell<Nlambda; ell++)
    {
        double kappaabsrho = 0.0;
        for (int h=0; h<_Ncomp; h++)
        {
            double kappaabs = mix(h)->kappaabs(ell);
            double rho = density(m,h);
            kappaabsrho += kappaabs*rho;
        }
        double J = absorbedLuminosity(m,ell) / (kappaabsrho*fac) / lambdagrid->dlambda(ell);
        // guard against (rare) situations where both Labs and kappa*fac are zero
        Jv[ell] = std::isfinite(J) ? J : 0.0;
    }
    return Jv;
}

//////////////////////////////////////////////////////////////////////

void DustSystem::fillOpticalDepth(PhotonPackage* pp)
{
    // determine the path and store the geometric details in the photon package
    _grid->path(pp);

    // if such statistics are requested, keep track of the number of cells crossed
    if (_writeCellsCrossed)
    {
        std::unique_lock<std::mutex> lock(_crossedMutex);
        unsigned int index = pp->size();
        if (index >= _crossed.size()) _crossed.resize(index+1);
        _crossed[index] += 1;
    }

    // calculate and store the optical depth details in the photon package
    pp->fillOpticalDepth(KappaRho(this, pp->ell()));

    // verify that the result makes sense
    double tau = pp->tau();
    if (tau<0.0 || std::isnan(tau) || std::isinf(tau))
        throw FATALERROR("The optical depth along the path is not a positive number: tau = "
                         + StringUtils::toString(tau));
}

//////////////////////////////////////////////////////////////////////

double DustSystem::opticalDepth(PhotonPackage* pp, double distance)
{
    // determine the path and store the geometric details in the photon package
    _grid->path(pp);

    // if such statistics are requested, keep track of the number of cells crossed
    if (_writeCellsCrossed)
    {
        std::unique_lock<std::mutex> lock(_crossedMutex);
        unsigned int index = pp->size();
        if (index >= _crossed.size()) _crossed.resize(index+1);
        _crossed[index] += 1;
    }

    // calculate and return the optical depth at the specified distance
    return pp->opticalDepth(KappaRho(this, pp->ell()), distance);
}

////////////////////////////////////////////////////////////////////

void DustSystem::write() const
{
    // If requested, output statistics on the number of cells crossed
    if (_writeCellsCrossed)
    {
        // Create a text file
        TextOutFile file(this, "ds_crossed", "number of cells crossed");

        // Write the header
        file.writeLine("# total number of cells in grid: " + std::to_string(_Ncells));
        file.addColumn("number of cells crossed", 'd');
        file.addColumn("number of paths that crossed this number of cells", 'd');

        // Write the body
        int Nlines = _crossed.size();
        for (int index=0; index<Nlines; index++)
        {
            file.writeRow(vector<double>({ static_cast<double>(index), static_cast<double>(_crossed[index]) }));
        }
    }
}

//////////////////////////////////////////////////////////////////////
