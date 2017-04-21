/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DRAINESILICATEGRAINCOMPOSITION_HPP
#define DRAINESILICATEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The DraineSilicateGrainComposition class represents the optical and calorimetric properties of
    silicate dust grains. The optical properties are taken from Bruce Draine's website, the
    calorimetric properties follow the prescription of Draine & Li (2001). */
class DraineSilicateGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(DraineSilicateGrainComposition, GrainComposition, "a Draine silicate dust grains composition")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by dust mix classes that wish to hard-code the creation of
        a new grain composition object of this type (as opposed to creation through the ski file).
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), and
        its setup() function has been called. */
    explicit DraineSilicateGrainComposition(SimulationItem* parent);

    /** This static function calculates and returns the specific enthalpy of silicate according to
        equation (11) of Draine & Li (2001), properly integrated to obtain internal energy rather
        than heat capacity. It serves as a call-back for the calculateEnthalpyGrid() function
        defined in the GrainComposition class. It is declared public so that it can also be used by
        other grain composition classes. */
    static double enthalpyFunction(double T);

protected:
    /** This function reads the raw optical property data from a resource file, calculates the
        enthalpy data using the analytical function for silicate derived in Draine & Li 2001, and
        sets the bulk mass density to the standard value of 3000 kg m\f$^{-3}\f$ for silicate grains.
        */
    void setupSelfBefore() override;

    //====================== Identifying =====================

public:
    /** This function returns a brief human-readable identifier for the type of grain composition
        represented by the instance. */
    string name() const override;
};

////////////////////////////////////////////////////////////////////

#endif
