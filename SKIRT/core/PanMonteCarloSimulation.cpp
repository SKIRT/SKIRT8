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
    if (_pds && _pds->hasDustEmission() && _pds->dustDistribution()->mass() > 0)
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

    // Get the parameters controlling the self-absorption iteration
    int minIters = _pds->minIterations();
    int maxIters = _pds->maxIterations();
    double fractionOfStellar = _pds->maxFractionOfStellar();
    double fractionOfPrevious = _pds->maxFractionOfPrevious();

    // Initialize the total absorbed luminosity in the previous iteration
    double prevLabsdusttot = 0.;

    // Iterate over the maximum number of iterations; the loop body returns from the function
    // when convergence is reached after the minimum number of iterations have been completed
    for (int iter = 1; iter<=maxIters; iter++)
    {
        // Construct the dust emission spectra
        {
            TimeLogger logger(log(), "calculation of emission spectra for dust self-absorption iteration "
                                            + std::to_string(iter));
            _pds->calculateDustEmission();
        }

        // Shoot the photons
        {
            TimeLogger logger(log(), "photon shooting for dust self-absorption iteration " + std::to_string(iter));

            // Determine the bolometric luminosity that is absorbed in every cell (and that will hence be re-emitted).
            _Labsbolv = _pds->absorbedLuminosity();
            // Set the absorbed dust luminosity to zero in all cells
            _pds->resetDustAbsorption();

            // Perform dust self-absorption
            initProgress("dust self-absorption iteration " + std::to_string(iter));
            Parallel* parallel = find<ParallelFactory>()->parallel();
            if (_lambdagrid->assigner())
                parallel->call(this, &PanMonteCarloSimulation::doDustSelfAbsorptionChunk,
                               _lambdagrid->assigner(), _Nchunks);
            else
                parallel->call(this, &PanMonteCarloSimulation::doDustSelfAbsorptionChunk, _Nlambda, _Nchunks);

            // Wait for the other processes to reach this point
            communicator()->wait("this self-absorption iteration");
            _pds->sumResults();
        }

        // Determine and log the total absorbed luminosity
        double Labsstellartot = _pds->absorbedStellarLuminosity();
        double Labsdusttot = _pds->absorbedDustLuminosity();
        log()->info("The total absorbed stellar luminosity is "
                   + StringUtils::toString(units()->obolluminosity(Labsstellartot)) + " "
                   + units()->ubolluminosity() );
        log()->info("The total absorbed dust luminosity in iteration " + std::to_string(iter) + " is "
                   + StringUtils::toString(units()->obolluminosity(Labsdusttot)) + " "
                   + units()->ubolluminosity() );

        // Log the current performance and corresponding convergence criteria
        if (Labsstellartot > 0. && Labsdusttot > 0.)
        {
            if (iter == 1)
            {
                log()->info("--> absorbed dust luminosity is "
                            + StringUtils::toString(Labsdusttot/Labsstellartot*100., 'f', 2)
                            + "% of absorbed stellar luminosity (convergence criterion is "
                            + StringUtils::toString(fractionOfStellar*100., 'f', 2) + "%)");
            }
            else
            {
                log()->info("--> absorbed dust luminosity changed by "
                            + StringUtils::toString(abs((Labsdusttot-prevLabsdusttot)/Labsdusttot)*100., 'f', 2)
                            + "% compared to previous iteration (convergence criterion is "
                            + StringUtils::toString(fractionOfPrevious*100., 'f', 2) + "%)");
            }
        }

        // Force at least the minimum number of iterations
        if (iter < minIters)
        {
            log()->info("Continuing until " + std::to_string(minIters) + " iterations have been performed");
        }
        else
        {
            // The self-absorption iteration has reached convergence if one or more of the following conditions holds:
            // - the absorbed stellar luminosity is zero
            // - the absorbed dust luminosity is zero
            // - the absorbed dust luminosity is less than a given fraction of the absorbed stellar luminosity
            // - the absorbed dust luminosity has changed by less than a given fraction compared to the previous iter
            if (Labsstellartot <= 0. || Labsdusttot <= 0.
                || Labsdusttot/Labsstellartot < fractionOfStellar
                || abs((Labsdusttot-prevLabsdusttot)/Labsdusttot) < fractionOfPrevious)
            {
                log()->info("Convergence reached after " + std::to_string(iter) + " iterations");
                return; // end the iteration by returning from the function
            }
            else
            {
                log()->info("Convergence not yet reached after " + std::to_string(iter) + " iterations");
            }
        }
        prevLabsdusttot = Labsdusttot;
    }

    // If the loop terminates, convergence was not reached even after the maximum number of iterations
    log()->error("Convergence not yet reached after " + std::to_string(maxIters) + " iterations");
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
    {
        TimeLogger logger(log(), "calculation of emission spectra for dust emission phase");
        log()->info("Calculating dust emission spectra...");
        _pds->calculateDustEmission();
    }

    // Shoot the photons
    {
        TimeLogger logger(log(), "photon shooting for dust emission phase");

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
