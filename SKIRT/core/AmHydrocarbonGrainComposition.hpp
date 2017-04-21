/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef AMHYDROCARBONGRAINCOMPOSITION_HPP
#define AMHYDROCARBONGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The AmHydrocarbonGrainComposition class represents the optical properties of a-C(:H) dust grains
    (carbonaceous grains) from Jones et al. 2013 (A&A, 558, A62). The calorimetric properties are
    calculated in DustEM. */
class AmHydrocarbonGrainComposition: public GrainComposition
{
    ITEM_CONCRETE(AmHydrocarbonGrainComposition, GrainComposition,
                  "Jones et al. 2013 amorphous hydrocarbon dust grain composition")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by dust mix classes that wish to hard-code the creation of
        a new grain composition object of this type (as opposed to creation through the ski file).
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), and
        its setup() function has been called. */
    explicit AmHydrocarbonGrainComposition(SimulationItem* parent);

protected:
    /** This function reads the raw optical and calorimetric data from resource files, and sets the
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
