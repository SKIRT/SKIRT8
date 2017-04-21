/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STELLARCOMP_HPP
#define STELLARCOMP_HPP

#include "SimulationItem.hpp"
class PhotonPackage;

//////////////////////////////////////////////////////////////////////

/** StellarComp is an abstract class used to represent a stellar component. A stellar
    component is characterized by a geometrical distribution of stars, the spectral energy
    distribution of each star, and the total luminosity of the system. */
class StellarComp : public SimulationItem
{
    ITEM_ABSTRACT(StellarComp, SimulationItem, "a stellar component")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the stellar component, which depends on the (lack of)
        symmetry of its geometry. A value of 1 means spherical symmetry, 2 means axial symmetry and
        3 means none of these symmetries. */
    virtual int dimension() const = 0;

    /** This function returns the luminosity \f$L_\ell\f$ of the stellar component in the
        wavelength bin at index \f$\ell\f$. */
    virtual double luminosity(int ell) const = 0;

    /** This function simulates the emission of a monochromatic photon package with a luminosity
        \f$L\f$ at wavelength index \f$\ell\f$ from the stellar component. */
    virtual void launch(PhotonPackage* pp, int ell, double L) const = 0;
};

//////////////////////////////////////////////////////////////////////

#endif
