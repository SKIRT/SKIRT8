/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LumFitN.hpp"
#include "FatalError.hpp"
#include "Image.hpp"
#include "GAPopulation.h"
#include "GARealGenome.h"
#include "GASStateGA.h"
#include <mutex>

////////////////////////////////////////////////////////////////////

void LumFitN::setMinLuminosities(const vector<double>& value)
{
    _minLum = value;
}

//////////////////////////////////////////////////////////////////////

void LumFitN::setMaxLuminosities(const vector<double>& value)
{
    _maxLum = value;
}

////////////////////////////////////////////////////////////////////

namespace
{
    double objectiveFunction(GAGenome& g, vector<Image>* sim)
    {
        // Read the array list and take out the reference image in the back
        Image ref = sim->back();
        sim->pop_back();
        GARealGenome& genome = (GARealGenome&)g;
        if (sim->size() != static_cast<size_t>(genome.size()))
            throw FATALERROR("Number of luminosities and components do not match");

        // Read the suggested luminosities for each component
        vector<double> lumis;
        for (int i = 0; i < genome.length(); i++)
        {
            lumis.push_back(pow(10,genome.gene(i)));
        }

        // Determine the chi2 value for this genome
        int arraysize = (*sim)[0].size();
        double chi=0;

        for (int m=0; m<arraysize; m++)
        {
            double total_sim = 0;

            // Create the correct summed image and take over masks from the reference image
            for (size_t n=0; n<sim->size(); n++)
            {
                if (ref[m]==0)
                {
                    ((*sim)[n])[m] = 0;
                }
                else
                {
                    total_sim += lumis[n] * (((*sim)[n])[m]);
                }
            }

            // Calculate the chi2 value if for non masked regions
            double sigma = sqrt(abs(ref[m]) + total_sim);
            if (ref[m]==0)
            {
                total_sim = 0;
                sigma = 0;
            }
            else
            {
                chi += pow( (ref[m] - total_sim) / sigma,2);
            }
        }

        // Returning the reference frame in the back of the list
        sim->push_back(ref);
        return chi;
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    void evaluatorFunction(GAPopulation& p)
    {
        vector<Image>* sim = (vector<Image>*)p.userData();

        // loop over all individuals and make replacement for all unevaluated individuals
        for (int i=0; i<p.size(); i++)
        {
            if (p.individual(i).isEvaluated()==gaFalse)
            {
                double value =  objectiveFunction(p.individual(i),sim);
                p.individual(i).score(value);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // Mutex to guard GALib operations
    std::mutex _mutex;
}

////////////////////////////////////////////////////////////////////

void LumFitN::optimize(const Image& refframe, vector<Image>& frames, vector<double>& luminosities, double& chi2)
{
    // Perform GALib operations in a critical section because GALib is not thread-safe
    std::unique_lock<std::mutex> lock(_mutex);

    // Create the boundaries, set to be uniform in logscale
    GARealAlleleSetArray allelesetarray;
    for (size_t i = 0; i < _minLum.size(); i++)
    {
        GARealAlleleSet RealAllele(log10(_minLum[i]), log10(_maxLum[i]));
        allelesetarray.add(RealAllele);
    }

    // Set the initializers, mutator and crossover scheme
    GARealGenome* genome = new GARealGenome(allelesetarray);
    genome->initializer(GARealGenome::UniformInitializer);
    genome->mutator(GARealGaussianMutator);
    genome->crossover(GARealUniformCrossover);
    GASteadyStateGA* ga = new GASteadyStateGA(*genome);
    GASigmaTruncationScaling scaling;
    ga->minimize();
    GAPopulation popu = ga->population();
    frames.push_back(refframe);
    popu.userData(&frames);
    popu.evaluator(evaluatorFunction);

    // Set the population size and generations to scale with the number of components
    ga->population(popu);
    ga->populationSize(frames.size()*30);
    ga->nGenerations(frames.size()*20);
    ga->pMutation(0.03);
    ga->pCrossover(0.65);
    ga->scaling(scaling);
    ga->scoreFrequency(0);  // was 1e-22 implicitly converted to integer 0
    ga->selectScores(GAStatistics::AllScores);
    ga->flushFrequency(0);  // was 1e-22 implicitly converted to integer 0
    ga->initialize();

    // Loop over the GA until it is done and print out the best fitting values
    GARealGenome& best_genome = (GARealGenome&) ga->statistics().bestIndividual();
    while (!ga->done()) ga->step();
    for (size_t i = 0; i < _minLum.size(); i++) luminosities.push_back(pow(10, best_genome.gene(i)));
    chi2 = best_genome.score();
}

////////////////////////////////////////////////////////////////////
