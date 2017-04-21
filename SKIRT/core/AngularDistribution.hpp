/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ANGULARDISTRIBUTION_HPP
#define ANGULARDISTRIBUTION_HPP

#include "Position.hpp"

////////////////////////////////////////////////////////////////////

/** AngularDistribution is a pure interface representing a possibly wavelength- and
    position-dependent angular probability distribution. For each wavelength, and at each position
    in space, the probability distribution \f$P(\Omega)\,{\mathrm{d}}\Omega\f$ is normalized on the
    unit sphere as follows: \f[ \int P(\Omega) \,{\mathrm{d}}\Omega = \int_{\phi=0}^{2\pi}
    \int_{\theta=0}^{\pi} P(\theta,\phi)\sin\theta \,{\mathrm{d}}\theta \,{\mathrm{d}}\phi =
    4\pi\f] An object that implements the interface, when given a particular wavelength and
    position, must provide a means to obtain the corresponding probability \f$P(\Omega)\f$ for any
    direction \f$(\theta,\phi)\f$, and to generate a random direction \f$(\theta,\phi)\f$ drawn
    from the probability distribution \f$P(\Omega)\,{\mathrm{d}}\Omega\f$. */
class AngularDistribution
{
protected:
    /** The empty constructor for the interface. */
    AngularDistribution() { }

public:
    /** The empty destructor for the interface. */
    virtual ~AngularDistribution() { }

    /** This function returns the probability \f$P(\Omega)\f$ for a given direction
        \f$(\theta,\phi)\f$ at the specified wavelength and position. For an isotropic
        distribution, this function would return 1 for any direction. */
    virtual double probabilityForDirection(int ell, Position bfr, Direction bfk) const = 0;

    /** This function generates a random direction \f$(\theta,\phi)\f$ drawn from the probability
        distribution \f$P(\Omega)\,{\mathrm{d}}\Omega\f$ at the specified wavelength and position.
        */
    virtual Direction generateDirection(int ell, Position bfr) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
