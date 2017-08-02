/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef AMHYDROCARBONGRAINCOMPOSITION_HPP
#define AMHYDROCARBONGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The AmHydrocarbonGrainComposition class represents amorphous carbonaceous dust grains described
    by Jones et al. 2017 (A&A, 602, A46) and the references therein. The optical and calorimetric
    properties are loaded from data files calculated with DustEM. */
class AmHydrocarbonGrainComposition: public GrainComposition
{
    ITEM_CONCRETE(AmHydrocarbonGrainComposition, GrainComposition,
                  "a Jones et al. 2017 amorphous hydrocarbon dust grain composition")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by dust mix classes that wish to hard-code the creation of
        a new grain composition object of this type (as opposed to creation through the ski file).
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), and its
        setup() function has been called. The bulk mass density is set to the value specified in
        the constructor rather than the value hardcoded in setupSelfBefore(). */
    explicit AmHydrocarbonGrainComposition(SimulationItem* parent, double bulkDensity);

protected:
    /** This function reads the optical and calorimetric data from resource files, and sets the
        bulk mass density to the value of 1.6 g/cm\f$^3\f$ specified by Jones. */
    void setupSelfBefore() override;

    //====================== Identifying =====================

public:
    /** This function returns a brief human-readable identifier for the type of grain composition
        represented by the instance. */
    string name() const override;
};

////////////////////////////////////////////////////////////////////

#endif
