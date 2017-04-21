/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP

#include "SimulationItem.hpp"
#include "Image.hpp"
#include "GAPopulation.h"
#include "GARealGenome.h"
#include "GASStateGA.h"
#include "SerializedData.hpp"

////////////////////////////////////////////////////////////////////

/** The optimization class contains all information necessary to do the actual fitting. This class
    uses the genetic algorithm library, GAlib. The ParameterRanges object residing in this fit
    scheme is used to set the boundaries and to interpret the output values. The
    evaluatePopulation() function in this class is called-back by GAlib. It feeds the genome values
    to the adjustable SKIRT simulation residing in this fit scheme, and calculates the goal
    function from the simulation results. This is done in parallalel for all individuals over the
    available threads or processes. */
class Optimization: public SimulationItem
{
    ITEM_CONCRETE(Optimization, SimulationItem, "The optimization setup")

    PROPERTY_INT(numIndividuals, "the number of individuals in the population")
        ATTRIBUTE_MIN_VALUE(popsize, "0")
        ATTRIBUTE_MAX_VALUE(popsize, "1000")
        ATTRIBUTE_DEFAULT_VALUE(popsize, "100")

    PROPERTY_INT(numGenerations, "the number of generations to be evaluated ")
        ATTRIBUTE_MIN_VALUE(generations, "0")
        ATTRIBUTE_MAX_VALUE(generations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(generations, "100")

    PROPERTY_DOUBLE(mutationProbability, "the mutation probability")
        ATTRIBUTE_MIN_VALUE(pmut, "]0")
        ATTRIBUTE_MAX_VALUE(pmut, "1[")
        ATTRIBUTE_DEFAULT_VALUE(pmut, "0.03")

    PROPERTY_DOUBLE(crossoverProbability, "the crossover probability")
        ATTRIBUTE_MIN_VALUE(pcross, "]0")
        ATTRIBUTE_MAX_VALUE(pcross, "1[")
        ATTRIBUTE_DEFAULT_VALUE(pcross, "0.65")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor releases the GA resources created during setup. */
    ~Optimization();

protected:
    /** This function readies the GA library for use, without final initialization. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** Initializes the GA library. */
    void initialize();

    /** Proceed one step in the GA optimization process. */
    void step();

    /** Checks if the GA optimization process is done. */
    bool done();

    /** Evaluates all individuals of a certain population. This is done by creating a temporary
        folder to store all simulations. The individual evaluations are parallelised over the
        available number of threads or processes and the function values are stored. At the end of
        each generation the contents of the temporary folder is removed, the scores for each
        individual are set and the best solutions are stored. */
    void evaluatePopulation(GAPopulation& pop);

private:
    /** Translates input/output variables to/from SerializedData and performs the SKIRT simulations
        in parallel. */
    void performSimulations();

    /** Performs the SKIRT simulation corresponding to the serialized input data, calculates the
        \f$\chi^2\f$ values and luminosities for the simulation result, and returns them in a
        serialized data object. */
    SerializedData performSimulation(const SerializedData& input);

    /** Reads the simulation output frames for the specified individual into the given table. There
        is a frame for each luminosity component (inner index) and for each wavelength (outer
        index). */
    void readSimulationOutputFrames(int individualIndex, vector<vector<Image>>& outputFrames);

    /** Write information about the specified individual as the best result so far. */
    void writeBest(int individualIndex);

    /** This function writes a series of comments lines to the specified output file, describing
        each of the columns written by the writeLine() function. The specified identifier
        description is used as the description of the first column. */
    void writeHeader(std::ofstream& stream, string identifierDescription);

    /** This function writes a summary line of information about the simulation result for the
        specified individual to the specified output file. The identifier is output as the first
        item of the line. */
    void writeLine(std::ofstream& stream, int identifier, int individualIndex);

    //======================== Data Members ========================

private:
    double _bestChi{DBL_MAX};
    int _bestSerial{0};
    GARealAlleleSetArray _allelesetarray;
    GARealGenome* _genome{nullptr};
    GASteadyStateGA* _ga{nullptr};
    std::ofstream _allstream;
    std::ofstream _beststream;
    vector<int> _genIndices;
    vector<vector<double>> _genValues;
    vector<vector<double>> _genUnitsValues;
    vector<double> _genScores;
    vector<vector<double>> _genLuminosities;
    vector<vector<double>> _genChis;
};

////////////////////////////////////////////////////////////////////

#endif
