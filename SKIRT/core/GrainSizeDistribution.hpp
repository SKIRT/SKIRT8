/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GRAINSIZEDISTRIBUTION_HPP
#define GRAINSIZEDISTRIBUTION_HPP

#include "SimulationItem.hpp"
#include "GrainSizeDistributionInterface.hpp"

////////////////////////////////////////////////////////////////////

/** GrainSizeDistribution is an abstract class that represents a size distribution for the dust
    grains in a particular dust population. Specifically, it represents a function \f[
    \Omega(a)=(\frac{\text{d}n_\text{D}}{\text{d}a})/n_\text{H} \qquad \text{for}\quad a_\text{min}
    \leq a \leq a_\text{max}, \f] that specifies the number of dust grains with size \f$a\f$ per
    hydrogen atom.

    The GrainSizeDistribution class publishes the GrainSizeDistributionInterface that provides
    access to the size distribution range and function. It expects each subclass to implement the
    functions declared in this interface, i.e. the functions amin() and amax() to specify the grain
    size range, and the function dnda() to specify the grain size distribution function within
    that range. For historical reasons, the latter function is named dnda() while it in fact
    returns the value of \f$\Omega(a)\f$ defined above.

    This base class manages the attribute \f$C\f$, a proportionality factor that should be used by
    subclasses as front factor in the function \f$\Omega(a)\f$.
*/
class GrainSizeDistribution: public SimulationItem, public GrainSizeDistributionInterface
{
    ITEM_ABSTRACT(GrainSizeDistribution, SimulationItem, "a dust grain size distribution")

    PROPERTY_DOUBLE(proportionalityFactor, "the proportionality factor in the size distribution function")
        ATTRIBUTE_MIN_VALUE(proportionalityFactor, "]0")
        ATTRIBUTE_DEFAULT_VALUE(proportionalityFactor, "1")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** Sets the proportionality factor \f$C\f$ in the size distribution function. */
    void setProportionalityFactor(double value);
};

////////////////////////////////////////////////////////////////////

#endif
