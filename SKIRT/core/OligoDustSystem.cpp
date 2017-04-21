/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OligoDustSystem.hpp"
#include "FatalError.hpp"
#include "FITSInOut.hpp"
#include "LockFree.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void OligoDustSystem::setupSelfAfter()
{
    DustSystem::setupSelfAfter();
    if (_writeMeanIntensity) _Labsvv.resize(dustGrid()->numCells(), find<WavelengthGrid>()->numWavelengths());
}

//////////////////////////////////////////////////////////////////////

bool OligoDustSystem::hasDustEmission() const
{
    return false;
}

//////////////////////////////////////////////////////////////////////

bool OligoDustSystem::hasDustAbsorption() const
{
    return _writeMeanIntensity;
}

//////////////////////////////////////////////////////////////////////

void OligoDustSystem::absorb(int m, int ell, double DeltaL, bool ynstellar)
{
    if (!ynstellar)
        throw FATALERROR("It is impossible to absorb non-stellar radiation in an oligochromatic simulation");
    LockFree::add(_Labsvv(m,ell), DeltaL);
}

//////////////////////////////////////////////////////////////////////

double OligoDustSystem::absorbedLuminosity(int m, int ell) const
{
    return _Labsvv(m,ell);
}

////////////////////////////////////////////////////////////////////

// Private class to output a FITS file with the mean intensity of the radiation field
// in each of the coordinate planes (xy, xz, or yz).
namespace
{
    // The image size in each direction, in pixels
    const int Np = 1024;

    class WriteMeanIntensityCut : public ParallelTarget
    {
    private:
        // cached values initialized in constructor
        const OligoDustSystem* _ds;
        DustGrid* _grid;
        Units* _units;
        Log* _log;
        double xbase, ybase, zbase, xpsize, ypsize, zpsize, xcenter, ycenter, zcenter;
        Array lambdav;
        int Nlambda;

        // data members initialized in setup()
        bool xd, yd, zd;  // direction of coordinate plane (110, 101, 011)
        string plane;    // name of the coordinate plane (xy, xz, yz)

        // results vector, properly sized in constructor and initialized to zero in setup()
        Array Jv;

    public:
        // constructor
        WriteMeanIntensityCut(const OligoDustSystem* ds)
        {
            _ds = ds;
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

            lambdav = _ds->find<WavelengthGrid>()->lambdav();
            Nlambda = lambdav.size();
            Jv.resize(Np*Np*Nlambda);
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
            _log->info("Calculating mean intensity in the " + plane + " plane...");

            Jv = 0.0;  // initialize all values to zero to facilitate the code in body()
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
                if (m!=-1)
                {
                    const Array& JJv = _ds->meanIntensity(m);
                    for (int ell=0; ell<Nlambda; ell++)
                    {
                        int l = i + Np*j + Np*Np*ell;
                        Jv[l] = _units->osurfacebrightness(lambdav[ell], JJv[ell]);
                    }
                }
            }
        }

        // Write the results to a FITS file with an appropriate name
        void write()
        {
            string filename = "ds_J" + plane;
            string description = "mean intensity";
            FITSInOut::write(_ds, description, filename, Jv, Np, Np, Nlambda,
                             _units->olength(xd?xpsize:ypsize), _units->olength(zd?zpsize:ypsize),
                             _units->olength(xd?xcenter:ycenter), _units->olength(zd?zcenter:ycenter),
                             _units->usurfacebrightness(), _units->ulength());
        }
    };
}

////////////////////////////////////////////////////////////////////

void OligoDustSystem::write() const
{
    DustSystem::write();

    // Perform the calculations only at the root
    if (_writeMeanIntensity && find<PeerToPeerCommunicator>()->isRoot())
    {
        // Get the parallel engine
        Parallel* parallel = find<ParallelFactory>()->parallel();

        // Output map(s) along coordinate axes
        {
            // Construct a private class instance to do the work (parallelized)
            WriteMeanIntensityCut wt(this);

            // Get the dimension of the dust grid
            int dimDust = dustGrid()->dimension();

            // For the xy plane (always)
            {
                wt.setup(1,1,0);
                parallel->call(&wt, Np);
                wt.write();
            }

            // For the xz plane (only if dimension is at least 2)
            if (dimDust >= 2)
            {
                wt.setup(1,0,1);
                parallel->call(&wt, Np);
                wt.write();
            }

            // For the yz plane (only if dimension is 3)
            if (dimDust == 3)
            {
                wt.setup(0,1,1);
                parallel->call(&wt, Np);
                wt.write();
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
