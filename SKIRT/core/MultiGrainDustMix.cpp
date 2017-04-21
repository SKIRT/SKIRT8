/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiGrainDustMix.hpp"
#include "DustDistribution.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "GrainComposition.hpp"
#include "GrainSizeDistributionInterface.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

//////////////////////////////////////////////////////////////////////

void MultiGrainDustMix::addPopulations(const GrainComposition *gc, const GrainSizeDistributionInterface *gs, int Nbins)
{
    Log* log = find<Log>();
    Units* units = find<Units>();

    // determine the grain size bins, one per requested population
    if (Nbins<1) throw FATALERROR("There must be at least one dust population bin");
    double amin = gs->amin();
    double amax = gs->amax();
    Array aminv(Nbins);
    Array amaxv(Nbins);
    {
        aminv[0] = amin;
        amaxv[Nbins-1] = amax;
        double logamin = log10(amin);
        double logamax = log10(amax);
        double deltaloga = (logamax-logamin)/Nbins;
        for (int c=1; c<Nbins; c++)
        {
            double logaminc = logamin + c*deltaloga;
            amaxv[c-1] = aminv[c] = pow(10,logaminc);
        }
    }

    // for each dust population (i.e. for each grain size bin)
    string gcname = gc->name(); // name of the grain composition class
    for (int c=0; c<Nbins; c++)
    {
        int popIndex = numPopulations();  // "global" zero-based index of this population

        // log the grain size range for this population
        log->info("Adding dust population #" + std::to_string(popIndex) + " based on " + gcname);
        log->info("  Grain sizes range from amin = "
                  + StringUtils::toString(units->ograinsize(aminv[c])) + " "
                  + units->ugrainsize()
                  + " to amax = "
                  + StringUtils::toString(units->ograinsize(amaxv[c])) + " "
                  + units->ugrainsize());

        // add the grain size information to an output text file, if so requested
        if (_writeSize)
        {
            int h = find<DustDistribution>()->indexForMix(this);

            // Create a text file
            TextOutFile file(this, "ds_mix_" + std::to_string(h) + "_size", "grain size information", !popIndex);

            if (!popIndex)
            {
                file.addColumn("dust mix population index");
                file.addColumn("grain composition type");
                file.addColumn("minimum dust grain radius (" + units->ugrainsize() + ")");
                file.addColumn("average dust grain radius (" + units->ugrainsize() + ")");
                file.addColumn("maximum dust grain radius (" + units->ugrainsize() + ")");
            }
            double aminc = aminv[c];
            double amaxc = amaxv[c];
            double aavec = pow(10.0,(log10(aminc)+log10(amaxc))/2.0);
            file.writeLine(StringUtils::toString(popIndex) + ' '
                            + gcname + ' '
                            + StringUtils::toString(units->ograinsize(aminc), 'g') + ' '
                            + StringUtils::toString(units->ograinsize(aavec), 'g') + ' '
                            + StringUtils::toString(units->ograinsize(amaxc), 'g'));
        }

        // create an integration grid over grain size within this bin
        int Na = 201;      // # points in grid  (must be > 2)
        Array av(Na);      // "a" for each point
        Array dav(Na);     // "da" for each point
        Array dndav(Na);   // "dnda" for each point
        Array weightv(Na); // integration weight for each point (1/2 or 1)
        {
            double logamin = log10(aminv[c]);
            double logamax = log10(amaxv[c]);
            double dloga = (logamax-logamin)/(Na-1);
            for (int i=0; i<Na; i++)
            {
                av[i] = pow(10, logamin + i*dloga);
                dav[i] = av[i] * M_LN10 * dloga;
                dndav[i] = gs->dnda(av[i]);
                weightv[i] = 1.;
            }
            weightv[0] = weightv[Na-1] = 0.5;
        }

        // get the simulation's wavelength grid
        const Array& lambdav = simlambdav();
        int Nlambda = lambdav.size();

        // calculate the optical properties for each wavelength
        Array sigmaabsv(Nlambda);
        Array sigmascav(Nlambda);
        Array asymmparv(Nlambda);
        for (int ell=0; ell<Nlambda; ell++)
        {
            double lamdba = lambdav[ell];
            double sumsigmaabs = 0.0;
            double sumsigmasca = 0.0;
            double sumgsigmasca = 0.0;
            for (int i=0; i<Na; i++)
            {
                double area = M_PI * av[i] * av[i];
                double sigmaabs = area * gc->Qabs(lamdba ,av[i]);
                double sigmasca = area * gc->Qsca(lamdba, av[i]);
                double gsigmasca = sigmasca * gc->asymmpar(lamdba, av[i]);
                sumsigmaabs += weightv[i] * dndav[i] * sigmaabs * dav[i];
                sumsigmasca += weightv[i] * dndav[i] * sigmasca * dav[i];
                sumgsigmasca += weightv[i] * dndav[i] * gsigmasca * dav[i];
            }
            sigmaabsv[ell] = sumsigmaabs;
            sigmascav[ell] = sumsigmasca;
            asymmparv[ell] = sumsigmasca ? sumgsigmasca/sumsigmasca : 0.;
        }

        // calculate the total mass per hydrogen atom, and the norm of the integration
        // (to calculate the mean mass of a single grain)
        double mu = 0.;
        double norm = 0.;
        double bulkdensity = gc->bulkDensity();
        for (int i=0; i<Na; i++)
        {
            double volume = 4.0*M_PI/3.0 * av[i] * av[i] * av[i];
            mu += weightv[i] * dndav[i] * volume * bulkdensity * dav[i];
            norm += weightv[i] * dndav[i] * dav[i];
        }

        // add a dust population with these properties (without resampling)
        addPopulation(mu, sigmaabsv, sigmascav, asymmparv);

        // remember the additional multi-grain properties needed for enthalpy calculations
        _gcv.push_back(gc);
        _meanmassv.push_back(mu/norm);

        // if the grain composition supports polarization, then calculate and add the polarization properties
        if (gc->polarization())
        {
            int Ntheta = 181;
            Table<2> S11vv(Nlambda,Ntheta);
            Table<2> S12vv(Nlambda,Ntheta);
            Table<2> S33vv(Nlambda,Ntheta);
            Table<2> S34vv(Nlambda,Ntheta);
            for (int ell=0; ell<Nlambda; ell++)
            {
                double lambda = lambdav[ell];
                for (int t=0; t<Ntheta; t++)
                {
                    double theta = t * M_PI/(Ntheta-1);
                    for (int i=0; i<Na; i++)
                    {
                        double w = weightv[i] * dndav[i] * dav[i];
                        double S11, S12, S33, S34;
                        gc->Sxx(lambda, av[i], theta, S11, S12, S33, S34);
                        S11vv(ell,t) += w * S11;
                        S12vv(ell,t) += w * S12;
                        S33vv(ell,t) += w * S33;
                        S34vv(ell,t) += w * S34;
                    }
                }
            }
            addPolarization(S11vv, S12vv, S33vv, S34vv);
        }
    }
}

//////////////////////////////////////////////////////////////////////

namespace
{
    class GrainSizeDistributionAdapter : public GrainSizeDistributionInterface
    {
    private:
        double _amin, _amax;
        double (*_dnda)(double a);
    public:
        GrainSizeDistributionAdapter(double amin, double amax, double (*dnda)(double a))
            : _amin(amin), _amax(amax), _dnda(dnda)  { }
        double amin() const { return _amin; }
        double amax() const { return _amax; }
        double dnda(double a) const { return _dnda(a); }
    };
}

//////////////////////////////////////////////////////////////////////

void MultiGrainDustMix::addPopulations(const GrainComposition *gc,
                                       double amin, double amax, double (*dnda)(double a),
                                       int Nbins)
{
    GrainSizeDistributionAdapter gsa(amin, amax, dnda);
    addPopulations(gc, &gsa, Nbins);
}

//////////////////////////////////////////////////////////////////////

string MultiGrainDustMix::grainCompositionName(int c) const
{
    return _gcv[c]->name();
}

//////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::upperTemperature() const
{
    double upper = 0.;
    int Np = numPopulations();
    for (int c=0; c<Np; c++)
    {
        upper = max(upper, _gcv[c]->upperTemperature());
    }
    return upper;
}

//////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::enthalpy(double T, int c) const
{
    return _meanmassv[c] * _gcv[c]->specificEnthalpy(T);
}

//////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::meanMass(int c) const
{
    return _meanmassv[c];
}

//////////////////////////////////////////////////////////////////////
