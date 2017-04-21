/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef COMPDUSTDISTRIBUTION_HPP
#define COMPDUSTDISTRIBUTION_HPP

#include "DustDistribution.hpp"
#include "Array.hpp"
#include "DustComp.hpp"

//////////////////////////////////////////////////////////////////////

/** The CompDustDistribution class is a subclass of the DustDistribution class and represents dust
    distributions consisting of different dust components. The class is basically just a vector of
    pointers to objects of the DustComp class. */
class CompDustDistribution : public DustDistribution
{
    ITEM_CONCRETE(CompDustDistribution, DustDistribution, "a dust distribution composed of various dust components")

    PROPERTY_ITEM_LIST(components, DustComp, "the dust components")
        ATTRIBUTE_DEFAULT_VALUE(components, "DustComp")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function contructs a cumulative mass distribution over the dust components. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the dust distribution, which depends on the (lack
        of) symmetry in the geometries of its dust components. A value of 1 means spherical
        symmetry, 2 means axial symmetry and 3 means none of these symmetries. The dust component
        with the least symmetry (i.e. the highest dimension) determines the result for the whole
        distribution. */
    int dimension() const override;

    /** This function returns the number of dust components that are involved in the dust
        distribution. */
    int numComponents() const override;

    /** This function returns a pointer to the dust mixture corresponding to the \f$h\f$'th dust
        component. */
    DustMix* mix(int h) const override;

    /** This function returns the mass density \f$\rho_h({\bf{r}})\f$ of the \f$h\f$'th component
        of the dust distribution at the position \f${\bf{r}}\f$. */
    double density(int h, Position bfr) const override;

    /** This function returns the total mass density \f$\rho({\bf{r}})\f$ of the dust distribution
        at the position \f${\bf{r}}\f$. For a component-based dust distribution, it just sums the
        contribution of the different components. */
    double density(Position bfr) const override;

    /** This function generates a random position from the dust distribution. It randomly chooses a
        dust component from the normalized cumulative density distribution constructed during the setup
        phase. Then a position is generated for the selected component. */
    Position generatePosition() const override;

    /** This function returns the total mass of the \f$h\f$'th component of the dust distribution. */
    double mass(int h) const override;

    /** This function returns the total dust mass of the dust distribution. For a component-based
        dust distribution, it just sums the contribution of the different components. */
    double mass() const override;

    /** This function returns the X-axis surface density of the dust distribution. For a
        component-based dust distribution, it just sums the contribution of the different
        components. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density of the dust distribution. For a
        component-based dust distribution, it just sums the contribution of the different
        components. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density of the dust distribution. For a
        component-based dust distribution, it just sums the contribution of the different
        components. */
    double SigmaZ() const override;

    /** This function is used by the interface() template function in the SimulationItem class. It
        returns a list of simulation items that should be considered in the search for an item that
        implements the requested interface. The implementation in this class returns the default
        list (i.e. the receiving CompDustDistribution instance) as the first item. If there is
        exactly one dust component, the geometry held by that dust component is added to
        the list. */
    virtual vector<SimulationItem*> interfaceCandidates(const std::type_info& interfaceTypeInfo) override;

    //======================== Data Members ========================

private:
    // initialized during setup
    Array _cumrhov;
};

//////////////////////////////////////////////////////////////////////

#endif
