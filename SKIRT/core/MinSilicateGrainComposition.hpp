/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MINSILICATEGRAINCOMPOSITION_HPP
#define MINSILICATEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The MinSilicateGrainComposition class represents the optical properties of amorphous silicate
    dust grains taken from Min et al. (2007, A&A, 462, 667). This model was designed to match the
    observations of interstellar dust towards the galactic center. The calorimetric properties are
    taken from the DustEM data included with SKIRT (see the DustEmGrainComposition class). */
class MinSilicateGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(MinSilicateGrainComposition, GrainComposition,
                  "a Min 2007 amorphous silicate dust grains composition")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the raw optical and calorimetric data from resource files, and sets the
        bulk mass density to the value of 3090 kg/m\f$^3\f$ specified by Min for silicate grains. */
    void setupSelfBefore() override;

    //====================== Identifying =====================

public:
    /** This function returns a brief human-readable identifier for the type of grain composition
        represented by the instance. */
    string name() const override;
};

////////////////////////////////////////////////////////////////////

#endif
