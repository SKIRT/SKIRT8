/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef POWERLAWGRAINSIZEDISTRIBUTION_HPP
#define POWERLAWGRAINSIZEDISTRIBUTION_HPP

#include "RangeGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** PowerLawGrainSizeDistribution is a GrainSizeDistribution subclass that represents a dust grain
    size distribution of the form \f[ \Omega(a)=(\frac{\text{d}n_\text{D}}{\text{d}a})/n_\text{H}
    = C \, a^{-\gamma} \qquad \text{for}\quad a_\text{min} \leq a \leq a_\text{max}, \f] where the
    the exponent \f$\gamma>0\f$ can be configured as an attribute. The size range and the
    proportionality factor \f$C\f$ can be configured in the GrainSizeDistribution base class. */
class PowerLawGrainSizeDistribution: public RangeGrainSizeDistribution
{
    ITEM_CONCRETE(PowerLawGrainSizeDistribution, RangeGrainSizeDistribution, "a power-law dust grain size distribution")

    PROPERTY_DOUBLE(exponent, "the (absolute value of the) exponent in the power-law distribution function")
        ATTRIBUTE_MIN_VALUE(exponent, "]0")
        ATTRIBUTE_MAX_VALUE(exponent, "10]")
        ATTRIBUTE_DEFAULT_VALUE(exponent, "3.5")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function implements the missing part of the GrainSizeDistributionInterface. It returns
        the value of \f$\Omega(a) = C\, a^{-\gamma}\f$. */
    double dnda(double a) const override;
};

////////////////////////////////////////////////////////////////////

#endif
