/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ENSTATITEGRAINCOMPOSITION_HPP
#define ENSTATITEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The EnstatiteGrainComposition class represents the optical properties of Enstatite dust grains,
    in two different forms:

    - \em Crystalline silicate MgSiO3 grains, for which Michiel Min <michiel.min.science@gmail.com>
      prepared the data. The refractive index data was taken from Jaeger et al. 1998. Further data was
      obtained from the Jena group (Fabian 2001, Zeidler 2011) for UV to near-IR. Below 0.2 micron the
      results are extrapolated using theoretical formulas. The calculations were performed with DHS
      using \f$f_\mathrm{max}=0.8\f$ (see Min et al. 2005). The calorimetric properties are taken from
      the DustEM data included with <tt>SKIRT</tt> (see the DustEmGrainComposition class).

    - \em Amorphous silicates with enstatite-normative composition from Köhler et al. 2014 (A&A,
      565, L9). Together with the amorphous silicates with forsterite-normative composition, they
      represent the silicate grains of Jones et al. 2017 (A&A, 602, A46). The optical and calorimetric
      properties are loaded from data files calculated for DustEM.  */
class EnstatiteGrainComposition : public GrainComposition
{
    /** The enumeration type indicating the type of Enstatite grains. */
    ENUM_DEF(GrainType, Crystalline, Amorphous)
    ENUM_VAL(GrainType, Crystalline, "crystalline (Min et al.)")
    ENUM_VAL(GrainType, Amorphous, "amorphous (Jones et al. 2017)")
    ENUM_END()

    ITEM_CONCRETE(EnstatiteGrainComposition, GrainComposition, "an Enstatite dust grain composition")

    PROPERTY_ENUM(grainType, GrainType, "the type of Enstatite grains")
        ATTRIBUTE_DEFAULT_VALUE(grainType, "Crystalline")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by dust mix classes that wish to hard-code the creation of
        a new grain composition object of this type (as opposed to creation through the ski file).
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), and
        its setup() function has been called. */
    explicit EnstatiteGrainComposition(SimulationItem* parent, GrainType type);

protected:
    /** This function reads the raw optical and calorimetric data from resource files, and sets the
        bulk mass density to the value of 2800 kg m\f$^{-3}\f$ specified by Min for crystalline enstatite
        and 2190 kg m\f$^{-3}\f$ specified by Köhler for amorphous enstatite. */
    void setupSelfBefore() override;

    //====================== Identifying =====================

public:
    /** This function returns a brief human-readable identifier for the type of grain composition
        represented by the instance. */
    string name() const override;
};

////////////////////////////////////////////////////////////////////

#endif
