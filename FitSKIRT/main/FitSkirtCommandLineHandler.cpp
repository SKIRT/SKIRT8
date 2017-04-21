/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FitSkirtCommandLineHandler.hpp"
#include "BuildInfo.hpp"
#include "Console.hpp"
#include "ConsoleHierarchyCreator.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "FitScheme.hpp"
#include "FitSchemeItemRegistry.hpp"
#include "ProcessManager.hpp"
#include "SchemaDef.hpp"
#include "StopWatch.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "TimeLogger.hpp"
#include "XmlHierarchyCreator.hpp"
#include "XmlHierarchyWriter.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // the allowed options list, in the format consumed by the CommandLineArguments constructor
    static const char* allowedOptions = "-t* -s* -i* -o* -k -x";
}

////////////////////////////////////////////////////////////////////

FitSkirtCommandLineHandler::FitSkirtCommandLineHandler()
    : _args(System::arguments(), allowedOptions)
{
    // issue welcome message
    _producerInfo = "FitSKIRT " + BuildInfo::projectVersion()
                    + " (" + BuildInfo::codeVersion() + " " + BuildInfo::timestamp() + ")";
    _hostUserInfo = "Running on " + System::hostname() + " for " + System::username();
    _console.info("Welcome to " + _producerInfo);
    _console.info(_hostUserInfo);
}

////////////////////////////////////////////////////////////////////

int FitSkirtCommandLineHandler::perform()
{
    // catch and properly report any exceptions
    try
    {
        // if there are no arguments at all --> interactive mode
        // if there is at least one file path argument --> batch mode
        // if the -x option is present --> export smile schema (undocumented option)
        // otherwise --> error
        if (_args.isValid() && !_args.hasOptions() && !_args.hasFilepaths()) return doInteractive();
        if (_args.hasFilepaths()) return doBatch();
        if (_args.isPresent("-x")) return doSmileSchema();
        _console.error("Invalid command line arguments");
        printHelp();
    }
    catch (FatalError& error)
    {
        for (string line : error.message()) _console.error(line);
    }
    catch (const std::exception& except)
    {
        _console.error("Standard Library Exception: " + string(except.what()));
    }
    return EXIT_FAILURE;
}

////////////////////////////////////////////////////////////////////

int FitSkirtCommandLineHandler::doInteractive()
{
    if (ProcessManager::isMultiProc()) throw FATALERROR("Interactive mode cannot be run with multiple processes");
    Console::info("Interactively constructing a fit scheme...");

    // ask for the name of the fski file in which to save the result
    string filename;
    while (true)
    {
        // get a file name, adding the .ski extension if needed
        filename = Console::promptForString("Enter the name of the fski file to be created", false, string());
        filename = StringUtils::addExtension(filename, "fski");

        // verify that the file does not yet exist
        // (we test whether the file can be opened, which is the best we can do in standard C++14)
        if (System::ifstream(filename)) Console::error("This file already exists; enter another name");
        else break;
    }

    // interactively construct the fit scheme
    auto schema = FitSchemeItemRegistry::getSchemaDef();
    auto fitscheme = ConsoleHierarchyCreator::create(schema);  // unique pointer

    // create the fski file reflecting this fit scheme
    XmlHierarchyWriter::write(fitscheme.get(), schema, filename, _producerInfo);

    // report success
    Console::info("Successfully created fski file '" + filename + "'.");
    Console::info("To run the fit use the command: fitskirt " + filename.substr(0, filename.length()-5));
    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////

int FitSkirtCommandLineHandler::doBatch()
{
    // build a list of paths to existing fski files
    vector<string> filepaths;
    for (string filepath : _args.filepaths())
    {
        // if the file does not exist as specified, try adding the .fski extension
        if (!System::ifstream(filepath)) filepath = StringUtils::addExtension(filepath, "fski");
        if (!System::ifstream(filepath))
        {
            _console.error("This ski file does not exist: " + filepath);
            return EXIT_FAILURE;
        }
        filepaths.push_back(filepath);
    }

    // a single filepath -> just run it
    if (filepaths.size() == 1)
    {
        doFitScheme(filepaths[0]);
    }

    // multiple filepaths -> add a time logger around the complete set
    else
    {
        TimeLogger logger(&_console, "a set of " + std::to_string(filepaths.size()) + " fit schemes");
        for (size_t index = 0; index < filepaths.size(); index++)
        {
            _console.warning("Performing fit scheme #" + std::to_string(index+1) +
                             " of " + std::to_string(filepaths.size()));
            doFitScheme(filepaths[index]);
        }
    }

    // report stopwatch results, if any
    for (string line : StopWatch::report()) _console.warning(line);
    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////

int FitSkirtCommandLineHandler::doSmileSchema()
{
    auto schema = FitSchemeItemRegistry::getSchemaDef();
    schema->save("fitskirt.smile", _producerInfo);
    _console.info("Successfully created SMILE schema file 'fitskirt.smile'.");
    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////

void FitSkirtCommandLineHandler::doFitScheme(string filepath)
{
    _console.info("Constructing a fit scheme from fski file '" + filepath + "'...");

    // construct the fit scheme from the fski file
    auto schema = FitSchemeItemRegistry::getSchemaDef();
    auto topitem = XmlHierarchyCreator::readFile(schema, filepath);  // unique pointer to Item
    auto fitscheme = dynamic_cast<FitScheme*>(topitem.get());

    // setup fit scheme attributes that are not loaded from the fski file:
    //  - file paths
    fitscheme->filePaths()->setOutputPrefix(StringUtils::filenameBase(filepath));
    string base = _args.isPresent("-k") ? StringUtils::dirPath(filepath) : "";
    string inpath = _args.value("-i");
    string outpath = _args.value("-o");
    if (!StringUtils::isAbsolutePath(inpath)) inpath = StringUtils::joinPaths(base, inpath);
    if (!StringUtils::isAbsolutePath(outpath)) outpath = StringUtils::joinPaths(base, outpath);
    fitscheme->filePaths()->setInputPath(inpath);
    fitscheme->filePaths()->setOutputPath(outpath);
    //  - parallel simulations and threads
    if (_args.intValue("-s") > 0) fitscheme->setParallelSimulationCount(_args.intValue("-s"));
    if (_args.intValue("-t") > 0) fitscheme->setParallelThreadCount(_args.intValue("-t"));

    // run the fit scheme
    fitscheme->setupAndRun();
}

////////////////////////////////////////////////////////////////////

void FitSkirtCommandLineHandler::printHelp()
{
    _console.warning("");
    _console.warning("To create a new fski file interactively:   fitskirt");
    _console.warning("To run a fit scheme with default options:  fitskirt <fski-filename>");
    _console.warning("");
    _console.warning(" fitskirt [-k] [-i <dirpath>] [-o <dirpath>]");
    _console.warning("          [-s <simulations>] [-t <threads>] ");
    _console.warning("          {<filepath>}*");
    _console.warning("");
    _console.warning("  -k : makes the input/output paths relative to the fski file being processed");
    _console.warning("  -i <dirpath> : the relative or absolute path for input files");
    _console.warning("  -o <dirpath> : the relative or absolute path for output files");
    _console.warning("  -s <simulations> : the number of parallel SKIRT simulations in single-process mode");
    _console.warning("  -t <threads> : the number of parallel threads for each SKIRT simulation");
    _console.warning("  <filepath> : the relative or absolute file path for an fski file");
    _console.warning("");
}

//////////////////////////////////////////////////////////////////////
