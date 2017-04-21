/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FORSTERITEGRAINCOMPOSITION_HPP
#define FORSTERITEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The ForsteriteGrainComposition class represents the optical properties of Forsterite dust
    grains, in two different forms:

    - \em Crystalline silicate Mg2SiO4 grains, for which Michiel Min <michiel.min.science@gmail.com>
      prepared the data. The refractive index data was taken from Suto et al. 2006, using the lowest
      temperature 50K. Further data was obtained from the Jena group (Fabian 2001, Zeidler 2011) for
      UV to near-IR. Below 0.2 micron the results are extrapolated using theoretical formulas. The
      calculations were performed with DHS using \f$f_\mathrm{max}=0.8\f$ (see Min et al. 2005). The
      calorimetric properties are taken from the DustEM data included with <tt>SKIRT</tt> (see the
      DustEmGrainComposition class).

    - \em Amorphous silicates with forsterite-normative composition from Köhler et al. 2014 (A&A,
      565, L9). Together with the amorphous silicates with enstatite-normative composition, they
      replace the silicate grains of Jones et al. 2013 (A&A, 558, A62). The calorimetric properties
      are calculated in DustEM. */
class ForsteriteGrainComposition : public GrainComposition
{
    /** The enumeration type indicating the type of Forsterite grains. */
    ENUM_DEF(GrainType, Crystalline, Amorphous)
    ENUM_VAL(GrainType, Crystalline, "crystalline")
    ENUM_VAL(GrainType, Amorphous, "amorphous")
    ENUM_END()

    ITEM_CONCRETE(ForsteriteGrainComposition, GrainComposition, "a Forsterite dust grains composition")

    PROPERTY_ENUM(grainType, GrainType, "the type of Forsterite grains")
        ATTRIBUTE_DEFAULT_VALUE(grainType, "Crystalline")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by dust mix classes that wish to hard-code the creation of
        a new grain composition object of this type (as opposed to creation through the ski file).
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), and
        its setup() function has been called. */
    explicit ForsteriteGrainComposition(SimulationItem* parent, GrainType type);

protected:
    /** This function reads the raw optical and calorimetric data from resource files, and sets the
        bulk mass density to the value of 3330 kg m\f$^{-3}\f$ specified by Min for crystalline Forsterite
        and 2190 kg m\f$^{-3}\f$ specified by Köhler for amorphous Forsterite. */
    void setupSelfBefore() override;

    //====================== Identifying =====================

public:
    /** This function returns a brief human-readable identifier for the type of grain composition
        represented by the instance. */
    string name() const override;
};

////////////////////////////////////////////////////////////////////

#endif
