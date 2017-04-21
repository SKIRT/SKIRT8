/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTGRIDDENSITYINTERFACE_HPP
#define DUSTGRIDDENSITYINTERFACE_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** DustGridDensityInterface is a pure interface. It is implemented by dust grids that offer a fast
    way to obtain the density of a grid cell, given the cell index \em m, for each dust component
    \em h. If the DustSystem class detects a dust grid that implements this interface, it will call
    the density() function in this interface rather than sampling the dust distribution in a number
    of random points. In some special cases, for example when the grid is lined up with some cell
    structure in the dust geometry, this can dramatically enhance performance. */
class DustGridDensityInterface
{
protected:
    /** The empty constructor for the interface. */
    DustGridDensityInterface() { }

public:
    /** The empty destructor for the interface. */
    virtual ~DustGridDensityInterface() { }

    /** This function returns the density for the dust component \em h in the dust grid cell with
        index \em m. */
    virtual double density(int h, int m) const = 0;
};

/////////////////////////////////////////////////////////////////////////////

#endif
