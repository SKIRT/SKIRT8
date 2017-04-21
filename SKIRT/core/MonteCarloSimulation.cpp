/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MonteCarloSimulation.hpp"
#include "DustDistribution.hpp"
#include "DustMix.hpp"
#include "DustSystem.hpp"
#include "FatalError.hpp"
#include "Instrument.hpp"
#include "InstrumentSystem.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "PhotonPackage.hpp"
#include "ProcessAssigner.hpp"
#include "Random.hpp"
#include "ShortArray.hpp"
#include "StellarSystem.hpp"
#include "StringUtils.hpp"
#include "TimeLogger.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setupSelfBefore()
{
    Simulation::setupSelfBefore();

    // just cache the pointers without requiring setup
    _lambdagrid = find<WavelengthGrid>(false);
    _ss = find<StellarSystem>(false);
    _ds = find<DustSystem>(false);  // dust system is optional
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setChunkParams(double packages)
{
    // Cache the number of wavelengths
    _Nlambda = _lambdagrid->numWavelengths();

    // Determine the number of chunks and the corresponding chunk size
    if (packages <= 0)
    {
        _Nchunks = 0;
        _chunksize = 0;
        _Npp = 0;
    }
    else
    {
        // Get the number of processes and threads per process
        int Nprocs = communicator()->size();
        int Nthreads = parallelFactory()->maxThreadCount();

        // Step 1: consider threading and determine the total number of chunks
        uint64_t totalChunks = 0;
        if (Nthreads == 1) totalChunks = 1;
        else totalChunks = static_cast<uint64_t>(ceil( max(10.*Nthreads*Nprocs/_Nlambda, packages/1e7) ));

        // Step 2: consider the work division and determine the number of chunks per process (_Nchunks)
        if (communicator()->dataParallel())  // Do some wavelengths for all chunks
        {
            _chunksize = static_cast<uint64_t>(ceil(packages/totalChunks));
            _Nchunks = totalChunks;
            _myTotalNpp = _lambdagrid->assigner()->assigned() * _Nchunks * _chunksize;
        }
        else                        // Do all wavelengths for some chunks
        {
            if ((totalChunks % Nprocs)) totalChunks = totalChunks + Nprocs - (totalChunks % Nprocs);
            _chunksize = static_cast<uint64_t>(ceil(packages/totalChunks));
            _Nchunks = totalChunks/Nprocs;
            _myTotalNpp = _Nlambda * _Nchunks * _chunksize;
        }

        // Calculate the the definitive number of photon packages per wavelength
        _Npp = totalChunks * _chunksize;

        log()->info("Using " + std::to_string(totalChunks) + " chunks per wavelength");
    }

    // Determine the log frequency; continuous scattering is much slower!
    _logchunksize = _continuousScattering ? 5000 : 50000;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setEmulationMode()
{
    _emulationMode = true;
    _numPackages = 0;
}

////////////////////////////////////////////////////////////////////

bool MonteCarloSimulation::emulationMode()
{
    return _emulationMode;
}

////////////////////////////////////////////////////////////////////

int MonteCarloSimulation::dimension() const
{
    return max(_ss->dimension(), _ds ? _ds->dimension() : 1);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::initProgress(string phase)
{
    _phase = phase;
    _Ndone = 0;

    log()->info(std::to_string(_Npp) + " photon packages for "
               + (_Nlambda==1 ? "a single wavelength" : "each of " + std::to_string(_Nlambda) + " wavelengths"));

    log()->infoSetElapsed(3);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::logProgress(uint64_t extraDone)
{
    // accumulate the work already done
    _Ndone.fetch_add(extraDone);

    // log message if the minimum time has elapsed
    double completed = _Ndone * 100. / _myTotalNpp;
    log()->infoIfElapsed("Launched " + _phase + " photon packages: " + StringUtils::toString(completed,'f',1) + "%");
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runStellarEmission()
{
    TimeLogger logger(log(), "the stellar emission phase");
    setChunkParams(_numPackages);
    initProgress("stellar emission");
    Parallel* parallel = find<ParallelFactory>()->parallel();

    if (_lambdagrid->assigner())
        parallel->call(this, &MonteCarloSimulation::doStellarEmissionChunk, _lambdagrid->assigner(), _Nchunks);
    else
        parallel->call(this, &MonteCarloSimulation::doStellarEmissionChunk, _Nlambda, _Nchunks);

    // Wait for the other processes to reach this point
    communicator()->wait("the stellar emission phase");
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::doStellarEmissionChunk(size_t index)
{
    int ell = index % _Nlambda;
    double L = _ss->luminosity(ell)/_Npp;
    if (L > 0)
    {
        double Lthreshold = L / minWeightReduction();
        PhotonPackage pp,ppp;

        uint64_t remaining = _chunksize;
        while (remaining > 0)
        {
            uint64_t count = min(remaining, _logchunksize);
            for (uint64_t i=0; i<count; i++)
            {
                _ss->launch(&pp,ell,L);
                if (pp.luminosity()>0)
                {
                    peelOffEmission(&pp,&ppp);
                    if (_ds) while (true)
                    {
                        _ds->fillOpticalDepth(&pp);
                        if (_continuousScattering) continuousPeelOffScattering(&pp,&ppp);
                        simulateEscapeAndAbsorption(&pp,_ds->hasDustAbsorption());
                        if (pp.luminosity()<=0 || (pp.luminosity()<=Lthreshold && pp.numScatt()>=_minScattEvents)) break;
                        simulatePropagation(&pp);
                        if (!_continuousScattering) peelOffScattering(&pp,&ppp);
                        simulateScattering(&pp);
                    }
                }
            }
            logProgress(count);
            remaining -= count;
        }
    }
    else logProgress(_chunksize);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::peelOffEmission(const PhotonPackage* pp, PhotonPackage* ppp)
{
    Position bfr = pp->position();

    for (Instrument* instr : _instrumentSystem->instruments())
    {
        Direction bfknew = instr->bfkobs(bfr);
        ppp->launchEmissionPeelOff(pp, bfknew);
        instr->detect(ppp);
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::peelOffScattering(const PhotonPackage* pp, PhotonPackage* ppp)
{
    int Ncomp = _ds->numComponents();
    int ell = pp->ell();
    Position bfr = pp->position();

    // Determine the weighting factors of the phase functions corresponding to
    // the different dust components: each component h is weighted by kappasca(h)*rho(m,h)
    ShortArray<4> wv(Ncomp);
    if (Ncomp==1)
        wv[0] = 1.0;
    else
    {
        int m = _ds->whichCell(bfr);
        if (m==-1) return; // abort peel-off
        for (int h=0; h<Ncomp; h++) wv[h] = _ds->mix(h)->kappasca(ell) * _ds->density(m,h);
        double sum = 0;
        for (int h=0; h<Ncomp; h++) sum += wv[h];
        if (sum<=0) return; // abort peel-off
        for (int h=0; h<Ncomp; h++) wv[h] /= sum;
    }

    // Now do the actual peel-off
    for (Instrument* instr : _instrumentSystem->instruments())
    {
        Direction bfkobs = instr->bfkobs(bfr);
        Direction bfkx = instr->bfkx();
        Direction bfky = instr->bfky();
        double I = 0, Q = 0, U = 0, V = 0;
        for (int h=0; h<Ncomp; h++)
        {
            DustMix* mix = _ds->mix(h);
            double w = wv[h] * mix->phaseFunctionValue(pp, bfkobs);
            StokesVector sv;
            mix->scatteringPeelOffPolarization(&sv, pp, bfkobs, bfkx, bfky);
            I += w * sv.stokesI();
            Q += w * sv.stokesQ();
            U += w * sv.stokesU();
            V += w * sv.stokesV();
        }
        ppp->launchScatteringPeelOff(pp, bfkobs, I);
        ppp->setPolarized(I, Q, U, V, pp->normal());
        instr->detect(ppp);
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::continuousPeelOffScattering(const PhotonPackage *pp, PhotonPackage *ppp)
{
    int ell = pp->ell();
    Position bfr = pp->position();
    Direction bfk = pp->direction();

    int Ncomp = _ds->numComponents();
    ShortArray<4> kappascav(Ncomp);
    ShortArray<4> kappaextv(Ncomp);
    for (int h=0; h<Ncomp; h++)
    {
        DustMix* mix = _ds->mix(h);
        kappascav[h] = mix->kappasca(ell);
        kappaextv[h] = mix->kappaext(ell);
    }

    int Ncells = pp->size();
    for (int n=0; n<Ncells; n++)
    {
        int m = pp->m(n);
        if (m!=-1)
        {
            ShortArray<4> wv(Ncomp);
            double ksca = 0.0;
            double kext = 0.0;
            for (int h=0; h<Ncomp; h++)
            {
                double rho = _ds->density(m,h);
                wv[h] = rho*kappascav[h];
                ksca += rho*kappascav[h];
                kext += rho*kappaextv[h];
            }
            if (ksca>0.0)
            {
                for (int h=0; h<Ncomp; h++) wv[h] /= ksca;
                double albedo = ksca/kext;
                double tau0 = (n==0) ? 0.0 : pp->tau(n-1);
                double dtau = pp->dtau(n);
                double s0 = (n==0) ? 0.0 : pp->s(n-1);
                double ds = pp->ds(n);
                double factorm = albedo * exp(-tau0) * (-expm1(-dtau));
                double s = s0 + random()->uniform()*ds;
                Position bfrnew(bfr+s*bfk);
                for (Instrument* instr : _instrumentSystem->instruments())
                {
                    Direction bfkobs = instr->bfkobs(bfrnew);
                    Direction bfkx = instr->bfkx();
                    Direction bfky = instr->bfky();
                    double I = 0, Q = 0, U = 0, V = 0;
                    for (int h=0; h<Ncomp; h++)
                    {
                        DustMix* mix = _ds->mix(h);
                        double w = wv[h] * mix->phaseFunctionValue(pp, bfkobs);
                        StokesVector sv;
                        mix->scatteringPeelOffPolarization(&sv, pp, bfkobs, bfkx, bfky);
                        I += w * sv.stokesI();
                        Q += w * sv.stokesQ();
                        U += w * sv.stokesU();
                        V += w * sv.stokesV();
                    }
                    ppp->launchScatteringPeelOff(pp, bfrnew, bfkobs, factorm*I);
                    ppp->setPolarized(I, Q, U, V, pp->normal());
                    instr->detect(ppp);
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::simulateEscapeAndAbsorption(PhotonPackage* pp, bool storeabsorptionrates)
{
    double taupath = pp->tau();
    int ell = pp->ell();
    double L = pp->luminosity();
    bool ynstellar = pp->isStellar();
    int Ncomp = _ds->numComponents();

    // Easy case: there is only one dust component
    if (Ncomp==1)
    {
        double albedo = _ds->mix(0)->albedo(ell);
        double expfactor = -expm1(-taupath);
        if (storeabsorptionrates)
        {
            int Ncells = pp->size();
            for (int n=0; n<Ncells; n++)
            {
                int m = pp->m(n);
                if (m!=-1)
                {
                    double taustart = (n==0) ? 0.0 : pp->tau(n-1);
                    double dtau = pp->dtau(n);
                    double expfactorm = -expm1(-dtau);
                    double Lintm = L * exp(-taustart) * expfactorm;
                    double Labsm = (1.0-albedo) * Lintm;
                    _ds->absorb(m,ell,Labsm,ynstellar);
                }
            }
        }
        double Lsca = L * albedo * expfactor;
        pp->setLuminosity(Lsca);
    }

    // Difficult case: there are different dust components.
    // The absorption/scattering in each cell is weighted by the density contribution of the component.
    else
    {
        Array kappascav(Ncomp);
        Array kappaextv(Ncomp);
        for (int h=0; h<Ncomp; h++)
        {
            DustMix* mix = _ds->mix(h);
            kappascav[h] = mix->kappasca(ell);
            kappaextv[h] = mix->kappaext(ell);
        }
        int Ncells = pp->size();
        double Lsca = 0.0;
        for (int n=0; n<Ncells; n++)
        {
            int m = pp->m(n);
            if (m!=-1)
            {
                double ksca = 0.0;
                double kext = 0.0;
                for (int h=0; h<Ncomp; h++)
                {
                    double rho = _ds->density(m,h);
                    ksca += rho*kappascav[h];
                    kext += rho*kappaextv[h];
                }
                double albedo = (kext>0.0) ? ksca/kext : 0.0;
                double taustart = (n==0) ? 0.0 : pp->tau(n-1);
                double dtau = pp->dtau(n);
                double expfactorm = -expm1(-dtau);
                double Lintm = L * exp(-taustart) * expfactorm;
                double Lscam = albedo * Lintm;
                Lsca += Lscam;
                if (storeabsorptionrates)
                {
                    double Labsm = (1.0-albedo) * Lintm;
                    _ds->absorb(m,ell,Labsm,ynstellar);
                }
            }
        }
        pp->setLuminosity(Lsca);
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::simulatePropagation(PhotonPackage* pp)
{
    double taupath = pp->tau();
    if (taupath==0.0) return;
    double tau = 0.0;
    if (_scattBias==0.0)
        tau = random()->exponCutoff(taupath);
    else
    {
        double X = random()->uniform();
        tau = (X<_scattBias) ? random()->uniform()*taupath : random()->exponCutoff(taupath);
        double p = -exp(-tau)/expm1(-taupath);
        double q = (1.0-_scattBias)*p + _scattBias/taupath;
        double weight = p/q;
        pp->setLuminosity(pp->luminosity()*weight);
    }
    double s = pp->pathLength(tau);
    pp->propagate(s);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::simulateScattering(PhotonPackage* pp)
{
    // Randomly select a dust mix; the probability of each dust component h is weighted by kappasca(h)*rho(m,h)
    DustMix* mix = _ds->randomMixForPosition(pp->position(), pp->ell());

    // Now perform the scattering using this dust mix
    Direction bfknew = mix->scatteringDirectionAndPolarization(pp, pp);
    pp->scatter(bfknew);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::write()
{
    TimeLogger logger(log(), "writing results");
    _instrumentSystem->write();
    if (_ds) _ds->write();
}

////////////////////////////////////////////////////////////////////
