# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library encapsulating the FFTW library for Fast Fourier
# Transforms, specifically providing an interface for convolving an
# array of data with a given kernel. The code in this library is
# inspired by a package called FFTConvolution developed by Jeremy Fix
# (see https://github.com/jeremyfix/FFTConvolution).
# ------------------------------------------------------------------

# set the target name
set(TARGET convolution)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} fundamentals)
include_directories(../../SMILE/fundamentals)

# define a user-configurable option to build with FFT convolution,
# which requires the FFTW3 library to be installed on the system
option(BUILD_WITH_FFT "build with FFT convolution - requires the FFTW3 library")

# if requested, find FFTW3 library and configure accordingly
if (BUILD_WITH_FFT)
    # find and configure fftw3
    find_library(FFT_LIBRARIES fftw3 PATHS ~/FFTW/lib)
    target_link_libraries(${TARGET} ${FFT_LIBRARIES})
    find_path(FFT_INCLUDES fftw3.h PATHS ~/FFTW/include)
    include_directories(${FFT_INCLUDES})
    # add the C++ preprocessor symbol that triggers compilation with MPI
    add_definitions(-DBUILD_WITH_FFT)
endif()

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
