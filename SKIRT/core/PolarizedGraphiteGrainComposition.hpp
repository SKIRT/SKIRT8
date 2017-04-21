/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef POLARIZEDGRAPHITEGRAINCOMPOSITION_HPP
#define POLARIZEDGRAPHITEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The PolarizedGraphiteGrainComposition class represents the optical and calorimetric properties
    of graphite dust grains with support for polarization. The optical data, including scattering and
    absorption efficiency coefficients and Mueller matrix coefficients, are read from an
    appropriate resource file in the format as used by the STOKES code. The calorimetric properties
    follow the prescription of Draine & Li (2001). */
class PolarizedGraphiteGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(PolarizedGraphiteGrainComposition, GrainComposition,
                  "a graphite dust grain composition with support for polarization")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the optical and calorimetric properties from the appropriate resource files
        and sets the appropriate bulk mass density value. */
    void setupSelfBefore() override;

    //====================== Identifying =====================

public:
    /** This function returns a brief human-readable identifier for the type of grain composition
        represented by the instance. */
    string name() const override;
};

////////////////////////////////////////////////////////////////////

#endif
