/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PanMonteCarloSimulation.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "PanDustSystem.hpp"
#include "PanWavelengthGrid.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "PhotonPackage.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "TimeLogger.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::setupSelfAfter()
{
    MonteCarloSimulation::setupSelfAfter();

    // properly size the array used to communicate between rundustXXX() and the corresponding parallel loop
    _Ncells = _pds ? _pds->numCells() : 0;
    if (_pds && _pds->hasDustEmission()) _Labsbolv.resize(_Ncells);
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::runSelf()
{
    runStellarEmission();
    if (_pds && _pds->hasDustEmission())
    {
        if (_pds && _pds->includeSelfAbsorption()) runDustSelfAbsorption();
        runDustEmission();
    }

    write();
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::runDustSelfAbsorption()
{
    TimeLogger logger(log(), "the dust self-absorption phase");

    // Initialize the total absorbed luminosity in the previous cycle
    double prevLabsdusttot = 0.;

    // Perform three "stages" of max 100 cycles each; the first stage uses 10 times less photon packages
    const int Nstages = 3;
    const char* stage_name[] = {"first-stage", "second-stage", "last-stage"};
    const double stage_factor[] = {1./10., 1./3., 1.};
    const double stage_epsmax[] = {0.010, 0.007, 0.005};
    for (int stage=0; stage<Nstages; stage++)
    {
        bool fixedNcycles = _pds->numCycles()!=0;
        const int Ncyclesmax = fixedNcycles ? _pds->numCycles() : 100;
        bool convergence = false;
        int cycle = 1;
        while (cycle<=Ncyclesmax && (!convergence || fixedNcycles))
        {
            TimeLogger logger(log(), "the " + string(stage_name[stage]) + " dust self-absorption cycle "
                              + std::to_string(cycle));

            // Construct the dust emission spectra
            log()->info("Calculating dust emission spectra...");
            _pds->calculateDustEmission();
            log()->info("Dust emission spectra calculated.");

            // Determine the bolometric luminosity that is absorbed in every cell (and that will hence be re-emitted).
            _Labsbolv = _pds->absorbedLuminosity();
            // Set the absorbed dust luminosity to zero in all cells
            _pds->resetDustAbsorption();

            // Perform dust self-absorption, using the appropriate number of packages for the current stage
            setChunkParams(numPackages()*stage_factor[stage]);
            initProgress(string(stage_name[stage]) + " dust self-absorption cycle " + std::to_string(cycle));

            Parallel* parallel = find<ParallelFactory>()->parallel();
            if (_lambdagrid->assigner())
                parallel->call(this, &PanMonteCarloSimulation::doDustSelfAbsorptionChunk,
                               _lambdagrid->assigner(), _Nchunks);
            else
                parallel->call(this, &PanMonteCarloSimulation::doDustSelfAbsorptionChunk, _Nlambda, _Nchunks);

            // Wait for the other processes to reach this point
            communicator()->wait("this self-absorption cycle");
            _pds->sumResults();

            // Determine and log the total absorbed luminosity in the vector Labstotv.
            double Labsdusttot = _pds->absorbedDustLuminosity();
            log()->info("The total absorbed stellar luminosity is "
                       + StringUtils::toString(units()->obolluminosity(_pds->absorbedStellarLuminosity())) + " "
                       + units()->ubolluminosity() );
            log()->info("The total absorbed dust luminosity is "
                       + StringUtils::toString(units()->obolluminosity(Labsdusttot)) + " "
                       + units()->ubolluminosity() );

            // Check the criteria to terminate the self-absorption cycle:
            // - the total absorbed dust luminosity should change by less than epsmax compared to the previous cycle;
            // - the last stage must perform at least 2 cycles (to make sure that the energy is properly distributed)
            double eps = fabs((Labsdusttot-prevLabsdusttot)/Labsdusttot);
            prevLabsdusttot = Labsdusttot;
            if ( (stage<Nstages-1 || cycle>1) && eps<stage_epsmax[stage])
            {
                log()->info("Convergence reached; the last increase in the absorbed dust luminosity was "
                           + StringUtils::toString(eps*100, 'f', 2) + "%");
                convergence = true;
            }
            else
            {
                log()->info("Convergence not yet reached; the increase in the absorbed dust luminosity was "
                           + StringUtils::toString(eps*100, 'f', 2) + "%");
            }
            cycle++;
        }
        if (!convergence && !emulationMode())
        {
            log()->error("Convergence not yet reached after " + std::to_string(Ncyclesmax) + " "
                        + string(stage_name[stage]) + " cycles!");
        }
    }
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::doDustSelfAbsorptionChunk(size_t index)
{
    // Determine the wavelength index for this chunk
    int ell = index % _Nlambda;

    // Determine the luminosity to be emitted at this wavelength index
    Array Lv(_Ncells);
    for (int m=0; m<_Ncells; m++)
    {
        double Labsbol = _Labsbolv[m];
        if (Labsbol>0.0) Lv[m] = Labsbol * _pds->emittedDustLuminosity(m,ell);
    }
    double Ltot = Lv.sum();

    // Emit photon packages
    if (Ltot > 0)
    {
        Array Xv;
        NR::cdf(Xv, Lv);

        PhotonPackage pp;
        double L = Ltot / _Npp;
        double Lthreshold = L / minWeightReduction();

        uint64_t remaining = _chunksize;
        while (remaining > 0)
        {
            uint64_t count = min(remaining, _logchunksize);
            for (uint64_t i=0; i<count; i++)
            {
                double X = random()->uniform();
                int m = NR::locateClip(Xv,X);
                Position bfr = _pds->randomPositionInCell(m);
                Direction bfk = random()->direction();
                pp.launch(L,ell,bfr,bfk);
                while (true)
                {
                    _pds->fillOpticalDepth(&pp);
                    simulateEscapeAndAbsorption(&pp,true);
                    double L = pp.luminosity();
                    if (L==0.0) break;
                    if (L<=Lthreshold && pp.numScatt()>=minScattEvents()) break;
                    simulatePropagation(&pp);
                    simulateScattering(&pp);
                }
            }
            logProgress(count);
            remaining -= count;
        }
    }
    else logProgress(_chunksize);
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::runDustEmission()
{
    TimeLogger logger(log(), "the dust emission phase");

    // Construct the dust emission spectra
    log()->info("Calculating dust emission spectra...");
    _pds->calculateDustEmission();
    log()->info("Dust emission spectra calculated.");

    // Determine the bolometric luminosity that is absorbed in every cell (and that will hence be re-emitted).
    _Labsbolv = _pds->absorbedLuminosity();

    // Perform the actual dust emission, possibly using more photon packages to obtain decent resolution
    setChunkParams(numPackages()*_pds->emissionBoost());
    initProgress("dust emission");
    Parallel* parallel = find<ParallelFactory>()->parallel();
    if (_lambdagrid->assigner())
        parallel->call(this, &PanMonteCarloSimulation::doDustEmissionChunk, _lambdagrid->assigner(), _Nchunks);
    else
        parallel->call(this, &PanMonteCarloSimulation::doDustEmissionChunk, _Nlambda, _Nchunks);

    // Wait for the other processes to reach this point
    communicator()->wait("the dust emission phase");
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::doDustEmissionChunk(size_t index)
{
    // Determine the wavelength index for this chunk
    int ell = index % _Nlambda;

    // Determine the luminosity to be emitted at this wavelength index
    Array Lv(_Ncells);
    for (int m=0; m<_Ncells; m++)
    {
        double Labsbol = _Labsbolv[m];
        if (Labsbol>0.0) Lv[m] = Labsbol * _pds->emittedDustLuminosity(m,ell);
    }
    double Ltot = Lv.sum();  // the total luminosity to be emitted at this wavelength index

    // Emit photon packages
    if (Ltot > 0)
    {
        // We consider biasing in the selection of the cell from which the photon packages are emitted.
        // A fraction of the cells is selected from the "natural" distribution, in which each cell is
        // weighted according to its total luminosity (Lv[em]). The other cells are selected from
        // a uniform distribution in which each cell has a equal probability.
        double xi = _pds->emissionBias();    // the fraction to be selected from a uniform distribution

        // the cumulative distribution of the natural pdf (Lv does not need to be normalized before)
        Array cumLv;
        NR::cdf(cumLv, Lv);

        PhotonPackage pp,ppp;
        double Lmean = Ltot/_Ncells;
        double Lem = Ltot / _Npp;
        double Lthreshold = Lem / minWeightReduction();

        uint64_t remaining = _chunksize;
        while (remaining > 0)
        {
            uint64_t count = min(remaining, _logchunksize);
            for (uint64_t i=0; i<count; i++)
            {
                int m;
                double X = random()->uniform();
                if (X<xi)
                {
                    // rescale the deviate from [0,xi[ to [0,Ncells[
                    m = max(0,min(_Ncells-1,static_cast<int>(_Ncells*X/xi)));
                }
                else
                {
                    // rescale the deviate from [xi,1[ to [0,1[
                    m = NR::locateClip(cumLv,(X-xi)/(1-xi));
                }
                double weight = 1.0/(1-xi+xi*Lmean/Lv[m]);
                Position bfr = _pds->randomPositionInCell(m);
                Direction bfk = random()->direction();
                pp.launch(Lem*weight,ell,bfr,bfk);
                peelOffEmission(&pp,&ppp);
                while (true)
                {
                    _pds->fillOpticalDepth(&pp);
                    if (continuousScattering()) continuousPeelOffScattering(&pp,&ppp);
                    simulateEscapeAndAbsorption(&pp,false);
                    double L = pp.luminosity();
                    if (L==0.0) break;
                    if (L<=Lthreshold && pp.numScatt()>=minScattEvents()) break;
                    simulatePropagation(&pp);
                    if (!continuousScattering()) peelOffScattering(&pp,&ppp);
                    simulateScattering(&pp);
                }
            }
            logProgress(count);
            remaining -= count;
        }
    }
    else logProgress(_chunksize);
}

////////////////////////////////////////////////////////////////////
