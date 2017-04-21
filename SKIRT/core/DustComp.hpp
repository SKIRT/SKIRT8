/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTCOMP_HPP
#define DUSTCOMP_HPP

#include "SimulationItem.hpp"
#include "DustCompNormalization.hpp"
#include "DustMix.hpp"
#include "Geometry.hpp"

//////////////////////////////////////////////////////////////////////

/** DustComp is a class that can be used to simulate a dust component. A dust component is
    characterized by a geometrical distribution of the dust, the optical properties of the dust,
    and the total amount of dust. A DustComp class object contains as data members a pointer to a
    Geometry object, a pointer to a dust mixture object of the DustMix class, and a
    normalization factor. */
class DustComp : public SimulationItem
{
    ITEM_CONCRETE(DustComp, SimulationItem, "a dust component")

    PROPERTY_ITEM(geometry, Geometry, "the geometry of the dust component")
        ATTRIBUTE_DEFAULT_VALUE(geometry, "ExpDiskGeometry")

    PROPERTY_ITEM(mix, DustMix, "the dust mixture of the dust component")
        ATTRIBUTE_DEFAULT_VALUE(mix, "InterstellarDustMix")

    PROPERTY_ITEM(normalization, DustCompNormalization, "the type of normalization for the dust component")
        ATTRIBUTE_DEFAULT_VALUE(normalization, "DustMassDustCompNormalization")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates the normalization factor based on the chosen normalization. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the dust component, which depends on the (lack of)
        symmetry of its geometry. A value of 1 means spherical symmetry, 2 means axial symmetry and
        3 means none of these symmetries. */
    int dimension();

    /** This function returns the mass density \f$\rho({\bf{r}})\f$ of the dust component at the
        position \f${\bf{r}}\f$. It is just the mass density of the Geometry object multiplied
        by the normalization factor. */
    double density(Position bfr) const;

    /** This function returns the total dust mass of the dust component at the position
        \f${\bf{r}}\f$. It is just the total dust mass of the Geometry object multiplied by the
        normalization factor. */
    double mass() const;

    /** This function returns the X-axis surface density of the dust component. It is just the
        X-axis surface density of the Geometry object multiplied by the normalization factor. */
    double SigmaX() const;

    /** This function returns the Y-axis surface density of the dust component. It is just the
        Y-axis surface density of the Geometry object multiplied by the normalization factor. */
    double SigmaY() const;

    /** This function returns the Z-axis surface density of the dust component. It is just the
        Z-axis surface density of the Geometry object multiplied by the normalization factor. */
    double SigmaZ() const;

    //======================== Data Members ========================

private:
    // initialized during setup
    double _nf{0.};
};

//////////////////////////////////////////////////////////////////////

#endif
