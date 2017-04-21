/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ZUBKOGRAPHITEGRAINSIZEDISTRIBUTION_HPP
#define ZUBKOGRAPHITEGRAINSIZEDISTRIBUTION_HPP

#include "GrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** ZubkoGraphiteGrainSizeDistribution represents the dust grain size distribution and grain size
    range for the graphite population in model BARE_GR_S of Zubko, Dwek & Arendt (2004, ApJS, 152,
    211). The proportionality factor \f$C\f$ configured in the GrainSizeDistribution base class is
    combined with the built-in front-factor. It should usually be set to the default value of
    \f$C=1\f$. */
class ZubkoGraphiteGrainSizeDistribution: public GrainSizeDistribution
{
    ITEM_CONCRETE(ZubkoGraphiteGrainSizeDistribution, GrainSizeDistribution,
                  "a Zubko, Dwek & Arendt size distribution for graphite dust grains")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by dust mix classes that wish to hard-code the creation of
        a new grain size distribution object of this type (as opposed to creation through the ski
        file). Before the constructor returns, the newly created object is hooked up as a child to
        the specified parent in the simulation hierarchy (so it will automatically be deleted), and
        it's setup() function has been called. The optional second argument specifies the
        proportionality factor configured in the GrainSizeDistribution base class, with a default
        value of 1. */
    explicit ZubkoGraphiteGrainSizeDistribution(SimulationItem* parent, double C = 1.);

    //======================== Other Functions =======================

public:
    /** This function implements part of the GrainSizeDistributionInterface. It returns the
        built-in minimum grain size \f$a_\text{min}\f$ as described in the header for this class.
        */
    double amin() const override;

    /** This function implements part of the GrainSizeDistributionInterface. It returns the
        built-in maximum grain size \f$a_\text{max}\f$ as described in the header for this class.
        */
    double amax() const override;

    /** This function implements part of the GrainSizeDistributionInterface. It returns the
        built-in value of \f$\Omega(a)\f$ as described in the header for this class. */
    double dnda(double a) const override;
};

////////////////////////////////////////////////////////////////////

#endif
