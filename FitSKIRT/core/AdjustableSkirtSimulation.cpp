/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdjustableSkirtSimulation.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "FitScheme.hpp"
#include "InstrumentFrame.hpp"
#include "InstrumentSystem.hpp"
#include "Log.hpp"
#include "MonteCarloSimulation.hpp"
#include "MultiFrameInstrument.hpp"
#include "OligoStellarComp.hpp"
#include "ParallelFactory.hpp"
#include "SimulationItemRegistry.hpp"
#include "StellarSystem.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "Units.hpp"
#include "XmlHierarchyCreator.hpp"
#include <streambuf>

////////////////////////////////////////////////////////////////////

void AdjustableSkirtSimulation::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // if the file does not exist as specified, try adding the .ski extension,
    string filepath = find<FilePaths>()->input(_skiFilename);
    std::ifstream infile = System::ifstream(filepath);
    if (!infile)
    {
        filepath = StringUtils::addExtension(filepath, "ski");
        infile = System::ifstream(filepath);
        if (!infile) throw FATALERROR("This ski file does not exist: " + filepath);
    }

    // read the complete file into our string
    _skiContent.assign(std::istreambuf_iterator<char>(infile), std::istreambuf_iterator<char>());
    infile.close();
    if (_skiContent.empty()) throw FATALERROR("Could not read the ski file " + filepath);

    // construct the simulation from the default ski content
    find<Log>()->info("Constructing default simulation hierarchy from ski file " + filepath + "...");
    auto schema = SimulationItemRegistry::getSchemaDef();
    auto topitem = XmlHierarchyCreator::readString(schema, _skiContent, _skiFilename);  // unique pointer to Item
    auto simulation = dynamic_cast<MonteCarloSimulation*>(topitem.get());

    // limit the number of threads to avoid multiple setup messages from Random
    simulation->parallelFactory()->setMaxThreadCount(1);

    // verify that the unit systems have the same type
    if (find<Units>()->type() != simulation->find<Units>()->type() ||
        find<Units>()->fluxOutputStyle() != simulation->find<Units>()->fluxOutputStyle())
        throw FATALERROR("Fit scheme and ski file must have the same type of unit system and flux output style");

    // copy information about the stellar system from the default simulation
    StellarSystem* stelsys = simulation->find<StellarSystem>();
    _ncomponents = stelsys->numComponents();
    for (auto stelcomp : stelsys->components())
    {
        _luminosities.push_back(stelcomp->find<OligoStellarComp>()->luminosities());
    }

    // copy information about the mandatory MultiFrameInstrument instrument from the default simulation
    MultiFrameInstrument* multiframe = simulation->find<InstrumentSystem>()->find<MultiFrameInstrument>();
    _instrname = multiframe->instrumentName();
    _nwavelengths = multiframe->frames().size();
    for (InstrumentFrame* insFrame : multiframe->frames())
    {
        _pixelsX.push_back(insFrame->numPixelsX());
        _pixelsY.push_back(insFrame->numPixelsY());
        _pixelSizeX.push_back(insFrame->fieldOfViewX()/insFrame->numPixelsX());
        _pixelSizeY.push_back(insFrame->fieldOfViewY()/insFrame->numPixelsY());
    }

    // inform the user
    find<Log>()->info("Number of stellar components in this simulation: " + std::to_string(_ncomponents));
    find<Log>()->info("Number of wavelengths in this simulation: " + std::to_string(_nwavelengths));
    find<Log>()->info("The simulation's instrument name is : " + _instrname);
}

////////////////////////////////////////////////////////////////////

size_t AdjustableSkirtSimulation::numComponents() const
{
    return _ncomponents;
}

////////////////////////////////////////////////////////////////////

size_t AdjustableSkirtSimulation::numWavelengths() const
{
    return _nwavelengths;
}

////////////////////////////////////////////////////////////////////

string AdjustableSkirtSimulation::instrumentName() const
{
    return _instrname;
}

////////////////////////////////////////////////////////////////////

int AdjustableSkirtSimulation::pixelsX(int ell) const
{
    return _pixelsX[ell];
}

////////////////////////////////////////////////////////////////////

int AdjustableSkirtSimulation::pixelsY(int ell) const
{
    return _pixelsY[ell];
}

////////////////////////////////////////////////////////////////////

double AdjustableSkirtSimulation::pixelSizeX(int ell) const
{
    return _pixelSizeX[ell];
}

////////////////////////////////////////////////////////////////////

double AdjustableSkirtSimulation::pixelSizeY(int ell) const
{
    return _pixelSizeY[ell];
}

////////////////////////////////////////////////////////////////////

double AdjustableSkirtSimulation::luminosity(int k, int ell) const
{
    return _luminosities[k][ell];
}

////////////////////////////////////////////////////////////////////

void AdjustableSkirtSimulation::performWith(const ReplacementDict& replacements, string prefix)
{
    // construct the simulation from the ski content after adjustment as requested
    auto schema = SimulationItemRegistry::getSchemaDef();
    auto topitem = XmlHierarchyCreator::readString(schema, adjustedSkiContent(replacements), _skiFilename);
    auto simulation = dynamic_cast<MonteCarloSimulation*>(topitem.get());

    // setup any simulation attributes that are not loaded from the ski content
    // -> copy file paths
    FilePaths* myfilepaths = find<FilePaths>();
    FilePaths* itsfilepaths = simulation->filePaths();
    itsfilepaths->setOutputPrefix(myfilepaths->outputPrefix() + "_" + prefix);
    itsfilepaths->setInputPath(myfilepaths->inputPath());
    itsfilepaths->setOutputPath(myfilepaths->outputPath());
    // -> copy number of threads
    simulation->parallelFactory()->setMaxThreadCount(find<FitScheme>()->parallelThreadCount());
    // -> suppress log messages
    simulation->log()->setLowestLevel(Log::Level::Error);

    // run the simulation
    simulation->setupAndRun();
}

////////////////////////////////////////////////////////////////////

string AdjustableSkirtSimulation::adjustedSkiContent(const ReplacementDict& replacements)
{
    Units* units = find<Units>();

    // declare shorthand for input string and container for output string
    const string& in{_skiContent};
    string out;

    // loop over input; index points at the next byte to be processed
    size_t index = 0;
    while (true)
    {
        // look for left bracket
        auto leftindex = in.find('[', index);
        if (leftindex == string::npos) break;

        // look for right bracket
        auto rightindex = in.find(']', leftindex+1);
        if (rightindex == string::npos) throw FATALERROR("Square brackets not balanced in ski file");

        // copy everything up to the left bracket and put the input index beyond the right bracket
        out += in.substr(index, leftindex-index);
        index = rightindex+1;

        // get the segment containing the attribute value, without the brackets
        string segment = in.substr(leftindex+1, rightindex-leftindex-1);
        if (segment.find('[') != string::npos) throw FATALERROR("Square brackets not balanced in ski file");

        // look for colon
        auto colonindex = segment.find(':');
        if (colonindex == string::npos) throw FATALERROR("Square brackets don't enclose colon in ski file");

        // get the label
        string label = segment.substr(0, colonindex);

        // if the label is in the replacements dict, insert the replacement value, otherwise fail
        if (replacements.count(label))
        {
            auto pair = replacements.at(label);
            double value = pair.first;
            string qty = pair.second;
            if (qty.empty()) out += StringUtils::toString(value);
            else out += StringUtils::toString(units->out(qty, value)) + " " + units->unit(qty);
        }
        else throw FATALERROR("No replacement found for label '" + label + "' in ski file");
    }

    // no more brackets -> copy the rest of the input
    out += in.substr(index);
    if (out.find(']') != string::npos) throw FATALERROR("Square brackets not balanced in ski file");

    return out;
}

////////////////////////////////////////////////////////////////////
