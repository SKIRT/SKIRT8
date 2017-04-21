/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FITSKIRTCOMMANDLINEHANDLER_HPP
#define FITSKIRTCOMMANDLINEHANDLER_HPP

#include "CommandLineArguments.hpp"
#include "ConsoleLog.hpp"

////////////////////////////////////////////////////////////////////

/**
This class processes the command line arguments for FitSKIRT and invokes the appropriate
high-level functions to perform the actions requested by the user. When invoked with invalid
command line arguments, it prints a brief help message. When invoked without any arguments, it
enters interactive mode, constructing a FitSKIRT fit scheme from the user's responses and
saving the result in an fski file, without actually performing the fit. Otherwise, it runs the
fit schemes in the fski files specified on the command line according to the following syntax:

\verbatim
    fitskirt [-k] [-i <dirpath>] [-o <dirpath>]
             [-s <simulations>] [-t <threads>]
             {<filepath>}*
\endverbatim

The first three options affect where FitSKIRT and the SKIRT simulations executed by FitSKIRT
look for input files and place output files. The -k option causes the input/output paths to
be relative to the fski file being processed, rather than to the current directory. The -i
option specifies the absolute or relative path for input files. The -o option specifies the
absolute or relative path for output files.

The following two options affect the parallization of SKIRT simulations executed by FitSKIRT.
The -s option specifies the number of SKIRT simulations to be executed in parallel (within a
single process). The -t option specifies the number of parallel threads for each SKIRT
simulation. For both options, the default value is one.

Finally, each \<filepath\> argument specifies the relative or absolute file path for a single
fski file, with or without the .fski extension. If there are multiple file paths, the specified
fski files are processed one after the other (i.e. serialized).
*/
class FitSkirtCommandLineHandler
{
public:
    /** This constructor obtains the program's command line arguments and issues a welcome message
        to the console log. */
    FitSkirtCommandLineHandler();

    /** This function processes the command line arguments and invokes the appropriate high-level
        functions to perform the actions requested by the user. The function returns an appropriate
        application exit value. */
    int perform();

private:
    /** This function conducts an interactive session to construct a fit scheme and save the result
        in an fski file. The function returns an appropriate application exit value. */
    int doInteractive();

    /** This function scans the filepaths specified on the command line for fski files and performs
        the corresponding fits according to the specified command line options. The function
        returns an appropriate application exit value. */
    int doBatch();

    /** This function exports a smile schema. This is an undocumented option. */
    int doSmileSchema();

    /** This function actually runs a single fit scheme constructed from the specified fski file. */
    void doFitScheme(string filepath);

    /** This function prints a brief help message to the console. */
    void printHelp();

private:
    // data members
    CommandLineArguments _args;
    ConsoleLog _console;
    string _producerInfo;
    string _hostUserInfo;
};

////////////////////////////////////////////////////////////////////

#endif
