# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library that supports genetic algorithm mechanisms.
# The source code is adapted from Matthew Wall's GAlib (c) 1995 MIT
# which can be downloaded at http://lancet.mit.edu/ga/GAlib.html
#---------------------------------------------------------------------

# set the target name
set(TARGET galib)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.h")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# adjust C++ compiler flags to our needs
set(NO_EXTRA_WARNINGS true)  # to avoid warnings in the legacy GALib code
include("../../SMILE/build/CompilerFlags.cmake")
