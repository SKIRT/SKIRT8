/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEMGRAINCOMPOSITION_HPP
#define DUSTEMGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The DustEmGrainComposition class represents the optical and calorimetric properties of dust
    grains, obtained from the DustEM data included with SKIRT. The user must provide the name of a
    particlar DustEM grain type as an attribute of this class. Example grain types include
    "PAH0_DL07", "PAH1_DL07", "Gra", and "aSil".

    DustEM is described in Compiègne et al. 2011 (AA, 525, A103) and the data was downloaded from
    http://www.ias.u-psud.fr/DUSTEM/. */
class DustEmGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(DustEmGrainComposition, GrainComposition, "a dust grain composition based on DustEM data")

    PROPERTY_STRING(grainType, "the DustEM grain type name")

    PROPERTY_DOUBLE(bulkMassDensity, "the bulk mass density for this grain material")
        ATTRIBUTE_QUANTITY(bulkMassDensity, "bulkmassdensity")
        ATTRIBUTE_MIN_VALUE(bulkMassDensity, "[100 kg/m3")
        ATTRIBUTE_MAX_VALUE(bulkMassDensity, "10000 kg/m3]")
        ATTRIBUTE_DEFAULT_VALUE(bulkMassDensity, "2240 kg/m3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by dust mix classes that wish to hard-code the creation of
        a new grain composition object of this type (as opposed to creation through the ski file).
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), the
        grain type attribute has been set to the value specified in the second argument, the bulk
        mass density has been set to the value specified in the third argument, and the setup()
        function has been called. */
    explicit DustEmGrainComposition(SimulationItem* parent, string graintype, double rhobulk);

protected:
    /** This function reads the optical and calorimetric properties from the DustEM resource files
        corresponding to the value of the grain type attribute, and sets the bulk mass density to
        the value of the corresponding attribute. */
    void setupSelfBefore() override;

    //====================== Identifying =====================

public:
    /** This function returns a brief human-readable identifier for the type of grain composition
        represented by the instance. */
    string name() const override;
};

////////////////////////////////////////////////////////////////////

#endif
