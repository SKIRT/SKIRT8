/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Dim1DustLib.hpp"
#include "ISRF.hpp"
#include "Log.hpp"
#include "PanDustSystem.hpp"
#include "StringUtils.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

int Dim1DustLib::numEntries() const
{
    return _numFieldStrengths;
}

////////////////////////////////////////////////////////////////////

vector<int> Dim1DustLib::mapping() const
{
    // get basic information about the wavelength grid and the dust system
    WavelengthGrid* lambdagrid = find<WavelengthGrid>();
    PanDustSystem* ds = find<PanDustSystem>();
    int Ncells = ds->numCells();

    // calculate the properties of the ISRF in all cells of the dust system;
    // remember the minimum and maximum values of the strength of the ISRF
    double JtotMW = ( ISRF::mathis(lambdagrid) * lambdagrid->dlambdav() ).sum();
    double Umin = DBL_MAX;
    double Umax = 0.0;
    vector<double> Ucellv(Ncells);
    for (int m=0; m<Ncells; m++)
    {
        double Jtot = ( ds->meanIntensity(m) * lambdagrid->dlambdav() ).sum();
        double U = Jtot/JtotMW;
        // ignore cells with extremely small radiation fields (compared to the average in the Milky Way)
        // to avoid wasting library grid points on fields that won't change simulation results anyway
        if (U > 1e-6)
        {
            Ucellv[m] = U;
            Umin = min(Umin,U);
            Umax = max(Umax,U);
        }
    }
    find<Log>()->info("ISRF strengths vary from U = " + StringUtils::toString(Umin)
                                         + " to U = " + StringUtils::toString(Umax) + ".");

    // determine for every dust cell m the corresponding library entry n
    vector<int> nv(Ncells);
    for (int m=0; m<Ncells; m++)
    {
        double U = Ucellv[m];
        if (U>0.0)
        {
            double logU = log10(U);
            double logUmin = log10(Umin);
            double logUmax = log10(Umax);
            double dlogU = (logUmax-logUmin)/_numFieldStrengths;
            nv[m] = max(0, min(_numFieldStrengths-1, static_cast<int>((logU-logUmin)/dlogU) ));
        }
        else
        {
            nv[m] = -1;
        }
    }

    return nv;
}

////////////////////////////////////////////////////////////////////
