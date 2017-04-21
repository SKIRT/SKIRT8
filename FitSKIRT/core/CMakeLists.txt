# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library for the FitSKIRT core layer
# ------------------------------------------------------------------

# set the target name
set(TARGET fitskirtcore)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} serialize schema fundamentals)
include_directories(../../SMILE/serialize ../../SMILE/schema ../../SMILE/fundamentals)

# add SKIRT library dependencies
target_link_libraries(${TARGET} skirtcore mpi utils)
include_directories(../../SKIRT/core ../../SKIRT/mpi ../../SKIRT/utils)

# add FitSKIRT library dependencies
target_link_libraries(${TARGET} convolution galib)
include_directories(../convolution ../galib)

# adjust C++ compiler flags to our needs
set(NO_EXTRA_WARNINGS true)  # to avoid warnings in the legacy GALib code
include("../../SMILE/build/CompilerFlags.cmake")
