/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Dim2DustLib.hpp"
#include "DustMix.hpp"
#include "Log.hpp"
#include "PanDustSystem.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

int Dim2DustLib::numEntries() const
{
    return _numTemperatures * _numWavelengths;
}

////////////////////////////////////////////////////////////////////

vector<int> Dim2DustLib::mapping() const
{
    // get basic information about the wavelength grid and the dust system
    WavelengthGrid* lambdagrid = find<WavelengthGrid>();
    PanDustSystem* ds = find<PanDustSystem>();
    int Nlambda = lambdagrid->numWavelengths();
    int Ncells = ds->numCells();
    int Ncomp = ds->numComponents();
    Log* log = find<Log>();
    Units* units = find<Units>();

    // calculate the properties of the ISRF in all cells of the dust system;
    // determine the minimum and maximum values of the mean temperature and mean wavelength
    double Tmin = DBL_MAX;
    double Tmax = 0.0;
    double lambdamin = DBL_MAX;
    double lambdamax = 0.0;
    Array Tmeanv(Ncells);
    Array lambdameanv(Ncells);
    for (int m=0; m<Ncells; m++)
    {
        if (ds->absorbedLuminosity(m) > 0.0)
        {
            const Array& Jv = ds->meanIntensity(m);
            double sumrho = 0.;
            for (int h=0; h<Ncomp; h++)
            {
                double sum0 = 0.0;
                double sum1 = 0.0;
                for (int ell=0; ell<Nlambda; ell++)
                {
                    double lambda = lambdagrid->lambda(ell);
                    double dlambda = lambdagrid->dlambda(ell);
                    double sigmaJ = ds->mix(h)->sigmaabs(ell) * Jv[ell];
                    sum0 += sigmaJ * dlambda;
                    sum1 += sigmaJ * lambda * dlambda;
                }
                double rho = ds->density(m,h);
                Tmeanv[m] += rho * ds->mix(h)->invplanckabs(sum0);
                lambdameanv[m] += rho * (sum1/sum0);
                sumrho += rho;
            }
            Tmeanv[m] /= sumrho;
            lambdameanv[m] /= sumrho;

            Tmin = min(Tmin,Tmeanv[m]);
            Tmax = max(Tmax,Tmeanv[m]);
            lambdamin = min(lambdamin,lambdameanv[m]);
            lambdamax = max(lambdamax,lambdameanv[m]);
        }
    }
    log->info("Temperatures vary"
              " from T = " + StringUtils::toString(units->otemperature(Tmin)) + " " + units->utemperature() +
              " to T = " + StringUtils::toString(units->otemperature(Tmax)) + " " + units->utemperature() + ".");
    log->info("Mean wavelengths vary"
              " from λ = " + StringUtils::toString(units->owavelength(lambdamin)) + " " + units->uwavelength() +
              " to λ = " + StringUtils::toString(units->owavelength(lambdamax)) + " " + units->uwavelength() + ".");

    // determine for every dust cell m the corresponding library entry n
    vector<int> nv(Ncells);
    for (int m=0; m<Ncells; m++)
    {
        if (Tmeanv[m] > 0.0 && lambdameanv[m] > 0.0)
        {
            double T = Tmeanv[m];
            double dT = (Tmax-Tmin)/_numTemperatures;
            int i = max(0, min(_numTemperatures-1, static_cast<int>((T-Tmin)/dT) ));

            double lambdamean = lambdameanv[m];
            double loglambdamean = log10(lambdamean);
            double loglambdamin = log10(lambdamin);
            double loglambdamax = log10(lambdamax);
            double dloglambdamean = (loglambdamax-loglambdamin)/_numWavelengths;
            int j = max(0, min(_numWavelengths-1, static_cast<int>((loglambdamean-loglambdamin)/dloglambdamean) ));

            nv[m] = i + _numTemperatures*j;
        }
        else
        {
            nv[m] = -1;
        }
    }

    return nv;
}

////////////////////////////////////////////////////////////////////
