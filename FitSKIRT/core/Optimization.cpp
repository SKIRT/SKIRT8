/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Optimization.hpp"

#include "AdjustableSkirtSimulation.hpp"
#include "Console.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "MasterSlaveCommunicator.hpp"
#include "OligoFitScheme.hpp"
#include "ParameterRanges.hpp"
#include "ReferenceImages.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "Units.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    void evaluate(GAPopulation& p)
    {
        Optimization* opt = (Optimization*)p.userData();
        opt->evaluatePopulation(p);
    }
}

//////////////////////////////////////////////////////////////////////

Optimization::~Optimization()
{
    delete _genome;
    delete _ga;
}

//////////////////////////////////////////////////////////////////////

void Optimization::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    ParameterRanges* ranges = find<ParameterRanges>();
    for (ParameterRange* range : ranges->ranges())
    {
        GARealAlleleSet RealAllele(range->minValue(),range->maxValue());
        _allelesetarray.add(RealAllele);
    }

    _genome = new GARealGenome(_allelesetarray);
    _genome->initializer(GARealGenome::UniformInitializer);
    _genome->mutator(GARealGaussianMutator);
    _genome->crossover(GARealUniformCrossover);
    _genome->userData(this);
    _ga = new GASteadyStateGA(*_genome);
    GASigmaTruncationScaling scaling;
    _ga->minimize();
    GAPopulation popu = _ga->population();
    popu.userData(this);
    popu.evaluator(evaluate);
    _ga->population(popu);
    _ga->populationSize(_numIndividuals);
    _ga->nGenerations(_numGenerations);
    _ga->pMutation(_mutationProbability);
    _ga->pCrossover(_crossoverProbability);
    _ga->scaling(scaling);
    _ga->scoreFrequency(0);
    _ga->selectScores(GAStatistics::AllScores);
    _ga->flushFrequency(0);

    // must be done before initializing the GA
    MasterSlaveCommunicator* comm = find<MasterSlaveCommunicator>();
    comm->setLocalSlaveCount(find<FitScheme>()->parallelSimulationCount());
    comm->registerTask([this] (SerializedData input) { return performSimulation(input); });
    if (comm->isMaster())
    {
        FilePaths* path = find<FilePaths>();
        _allstream = System::ofstream(path->output("all_simulations.dat"));
        _beststream = System::ofstream(path->output("best_simulations.dat"));
        writeHeader(_allstream, "generation index");
        writeHeader(_beststream, "serial number of best result");
    }
}

//////////////////////////////////////////////////////////////////////

void Optimization::initialize()
{
    _ga->initialize();
}

//////////////////////////////////////////////////////////////////////

void Optimization::step()
{
    _ga->step();
}

//////////////////////////////////////////////////////////////////////

bool Optimization::done()
{
   return _ga->done() != gaFalse;
}

//////////////////////////////////////////////////////////////////////

void Optimization::evaluatePopulation(GAPopulation& pop)
{
    auto generationIndex = pop.geneticAlgorithm()->generation();
    find<Log>()->info("Evaluating generation " + std::to_string(generationIndex));
    Units* units = find<Units>();
    auto ranges = find<ParameterRanges>()->ranges();

    // loop over all individuals and create replacement info for all unevaluated individuals
    _genIndices.clear();
    _genValues.clear();
    _genUnitsValues.clear();
    for (int i=0; i<pop.size(); ++i)
    {
        if (pop.individual(i).isEvaluated()==gaFalse)
        {
            GARealGenome& genome = (GARealGenome&)pop.individual(i);

            // loop over all parameters using the genome values to create the replacement info
            vector<double> currentValues;
            vector<double> currentUnitsValues;
            for (size_t j=0; j < ranges.size(); ++j)
            {
                double value = genome.gene(j);
                currentValues.push_back(value);
                string qty = ranges[j]->quantityString();
                if (!qty.empty()) value = units->out(qty, value);
                currentUnitsValues.push_back(value);
            }
            _genIndices.push_back(i);
            _genValues.push_back(currentValues);
            _genUnitsValues.push_back(currentUnitsValues);
        }
    }

    // create a temporary directory to store the SKIRT simulation results
    string tmpdirpath = find<FilePaths>()->output("tmp");
    if (!System::makeDir(tmpdirpath))
        throw FATALERROR("Can't create temporary directory " + tmpdirpath);

    // perform the simulations and calculate the objective function values in parallel
    performSimulations();

    // set the individual's scores and write a summary line for each simulation
    find<Log>()->info("Setting Scores");
    for (size_t i=0; i<_genIndices.size(); ++i)
    {
        pop.individual(_genIndices[i]).score(_genScores[i]);
        writeLine(_allstream, generationIndex, i);
    }

    // find the best (newly evaluated) individual in this population
    size_t bestIndex = 0;
    double bestChi = _genScores[0];
    for (size_t i=1; i<_genIndices.size(); ++i)
    {
        if (_genScores[i] < bestChi)
        {
             bestIndex = i;
             bestChi = _genScores[i];
        }
    }

    // if the best individual is the best one so far, write information about it
    if (bestChi < _bestChi)
    {
        _bestChi = bestChi;
        writeBest(bestIndex);
    }

    // remove all files in the temporary folder
    for (string filename : System::filesInDirectory(tmpdirpath))
    {
        System::removeFile(StringUtils::joinPaths(tmpdirpath, filename));
    }
}

//////////////////////////////////////////////////////////////////////

void Optimization::performSimulations()
{
    // serialize input data for each of the simulations to perform
    size_t n = _genValues.size();
    vector<SerializedData> datav(n);
    for (size_t i = 0; i < n; ++i)
    {
        datav[i].push(static_cast<double>(i));
        datav[i].push(_genValues[i]);
    }

    // perform the simulations in parallel
    MasterSlaveCommunicator* comm = find<MasterSlaveCommunicator>();
    datav = comm->performTask(datav);

    // deserialize output data from each of the simulations
    _genScores.resize(n);
    _genLuminosities.resize(n);
    _genChis.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        datav[i].pop(_genChis[i]);      // pop in reverse order!
        datav[i].pop(_genLuminosities[i]);
        datav[i].pop(_genScores[i]);
    }
}

//////////////////////////////////////////////////////////////////////

SerializedData Optimization::performSimulation(const SerializedData& input)
{
    // deserialize input data
    SerializedData data = input;
    vector<double> genValues;
    data.pop(genValues);
    size_t individualIndex = static_cast<size_t>(data.pop());

    // create a replacement dictionary from the input genome values
    AdjustableSkirtSimulation::ReplacementDict replacement;
    auto ranges = find<ParameterRanges>()->ranges();
    for (size_t j=0; j < ranges.size(); ++j)
    {
        replacement[ranges[j]->label()] = std::make_pair(genValues[j], ranges[j]->quantityString());
    }

    // perform the adjusted SKIRT simulation
    // HACK: we issue messages directly to the console, bypassing the regular mechanism,
    // to ensure that these messages are always logged even if sent from a slave process
    string individualString = std::to_string(individualIndex);
    string slaveString =  std::to_string(find<MasterSlaveCommunicator>()->slave());
    Console::info("  Slave " + slaveString + " running SKIRT model for individual " + individualString);
    auto simulation = find<AdjustableSkirtSimulation>();
    simulation->performWith(replacement, "tmp/tmp_" + std::to_string(individualIndex));
    Console::info("  Slave " + slaveString + " fitting luminosities for individual " + individualString);

    // read the simulation output frames
    vector<vector<Image>> frames;
    readSimulationOutputFrames(individualIndex, frames);

    // determine the best fitting luminosities and lowest chi2 value
    vector<vector<double>> luminosities;
    vector<double> chis;
    auto refImages = find<ReferenceImages>();
    double chi_sum = refImages->optimizeLuminosities(frames, luminosities, chis);

    // flatten the luminosities into a single vector
    vector<double> flatluminosities;
    for (const auto& vect : luminosities) for (double value : vect) flatluminosities.push_back(value);

    // serialize the output data
    data.push(chi_sum);
    data.push(flatluminosities);
    data.push(chis);
    return data;
}

//////////////////////////////////////////////////////////////////////

void Optimization::readSimulationOutputFrames(int individualIndex, vector<vector<Image>>& outputFrames)
{
    string prefix = "tmp_" + std::to_string(individualIndex);
    string tmpdirpath = find<FilePaths>()->output("tmp");
    auto simulation = find<AdjustableSkirtSimulation>();
    string instrname = simulation->instrumentName();
    size_t numWavelengths = simulation->numWavelengths();

    outputFrames.resize(numWavelengths);
    for (size_t ell=0; ell < numWavelengths; ++ell)
    {
        outputFrames[ell].clear();
        for (size_t k = 0; k < simulation->numComponents(); k++)
        {
            string filename = prefix + "_" + instrname + "_stellar_" + std::to_string(k) + "_" + std::to_string(ell);
            outputFrames[ell].emplace_back(this, filename, tmpdirpath, false);
        }
    }
}

//////////////////////////////////////////////////////////////////////

void Optimization::writeBest(int individualIndex)
{
    find<Log>()->info("Found new best fit; serial nr " + std::to_string(_bestSerial));

    // read the simulation output frames for the specified individual
    vector<vector<Image>> frames;
    readSimulationOutputFrames(individualIndex, frames);

    // determine the best fitting luminosities
    vector<vector<double>> luminosities;
    vector<double> chis; // not used
    auto refImages = find<ReferenceImages>();
    refImages->optimizeLuminosities(frames, luminosities, chis);

    // calculate the total and residual images
    vector<Image> totals;
    vector<Image> residuals;
    refImages->getTotalAndResidual(frames, luminosities, totals, residuals);

    // get information for saving
    auto simulation = find<AdjustableSkirtSimulation>();
    auto units = find<Units>();
    string xyUnits = units->ulength();
    string sbUnits = units->usurfacebrightness();

    // save the best fitting total and residual images
    for (size_t ell=0; ell < totals.size() ; ++ell)
    {
        string filename1 = "best_" + std::to_string(_bestSerial) + "_" + std::to_string(ell);
        totals[ell].save(this, "best fitting frame", filename1,
                         units->olength(simulation->pixelSizeX(ell)),
                         units->olength(simulation->pixelSizeY(ell)),
                         sbUnits, xyUnits);

        string filename2 =  "residual_" + std::to_string(_bestSerial) + "_" + std::to_string(ell);
        residuals[ell].save(this, "residuals frame", filename2,
                            units->olength(simulation->pixelSizeX(ell)),
                            units->olength(simulation->pixelSizeY(ell)),
                            string(), xyUnits);
    }

    // write a summary line in the "best results" file
    writeLine(_beststream, _bestSerial, individualIndex);
    _bestSerial++;
}

//////////////////////////////////////////////////////////////////////

void Optimization::writeHeader(std::ofstream& stream, string identifierDescription)
{
    // identifier
    size_t ncol = 0;
    stream << "# column " << ++ncol << ": " << identifierDescription << '\n';

    // parameter values
    auto units = find<Units>();
    auto ranges = find<ParameterRanges>()->ranges();
    for (size_t j=0; j < ranges.size(); ++j)
    {
        stream << "# column " << ++ncol << ": " << ranges[j]->label();
        string qty = ranges[j]->quantityString();
        if (qty.empty()) stream << '\n';
        else stream << " (" << units->unit(qty) << ")\n";
    }

    // luminosities
    auto simulation = find<AdjustableSkirtSimulation>();
    size_t numWavelengths = simulation->numWavelengths();
    size_t numComponents = simulation->numComponents();
    for (size_t ell=0; ell < numWavelengths; ++ell)
    {
        for (size_t k=0; k < numComponents; ++k)
        {
            stream << "# column " << ++ncol << ": luminosity for component " << k
                   << " at wavelength " << ell
                   << " (" << StringUtils::toString(units->omonluminosity(simulation->luminosity(k, ell)))
                   << ' ' << units->umonluminosity() << ")\n";
        }
    }

    // total chi^2 and chi^2 per wavelength
    stream << "# column " << ++ncol << ": total chi^2\n";
    for (size_t ell=0; ell < numWavelengths; ++ell)
    {
        stream << "# column " << ++ncol << ": chi^2 for wavelength " << ell << '\n';
    }

    stream.flush();
}

//////////////////////////////////////////////////////////////////////

void Optimization::writeLine(std::ofstream& stream, int identifier, int individualIndex)
{
    stream << identifier << "   ";
    for (double value : _genUnitsValues[individualIndex]) stream << value << " ";
    for (double value : _genLuminosities[individualIndex]) stream << value << " ";
    stream << "  " << _genScores[individualIndex] << " ";
    for (double value : _genChis[individualIndex]) stream << value << " ";
    stream << std::endl;
}

//////////////////////////////////////////////////////////////////////
