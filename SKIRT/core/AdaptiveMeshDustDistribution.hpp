/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHDUSTDISTRIBUTION_HPP
#define ADAPTIVEMESHDUSTDISTRIBUTION_HPP

#include "BoxDustDistribution.hpp"
#include "AdaptiveMeshInterface.hpp"
#include "AdaptiveMeshFile.hpp"
#include "Array.hpp"
#include "MeshDustComponent.hpp"

////////////////////////////////////////////////////////////////////

/** The AdaptiveMeshDustDistribution class represents a dust distribution imported from an adaptive
    mesh data file. The data file must have one of the supported formats; refer to the
    AdaptiveMeshFile class and its subclasses. Since the adaptive mesh data format does not specify
    the size of the domain, this information must be provided as properties of this class. This
    class supports multiple dust components, as long as the dust density distributions for all
    components are defined on the same mesh in the same adaptive mesh data file. Each dust
    component is represented by an instance of the MeshDustComponent class, which specifies the
    data column index defining the dust density distribution for the component and the
    corresponding dust mix. */
class AdaptiveMeshDustDistribution : public BoxDustDistribution, public AdaptiveMeshInterface
{
    ITEM_CONCRETE(AdaptiveMeshDustDistribution, BoxDustDistribution,
                  "a dust distribution imported from an adaptive mesh data file")

    PROPERTY_ITEM(adaptiveMeshFile, AdaptiveMeshFile, "the adaptive mesh data file")
        ATTRIBUTE_DEFAULT_VALUE(adaptiveMeshFile, "AdaptiveMeshAsciiFile")

    PROPERTY_DOUBLE(densityUnits, "the units in which the file specifies density values")
        ATTRIBUTE_QUANTITY(densityUnits, "massvolumedensity")
        ATTRIBUTE_MIN_VALUE(densityUnits, "]0")
        ATTRIBUTE_DEFAULT_VALUE(densityUnits, "1 Msun/pc3")

    PROPERTY_ITEM_LIST(components, MeshDustComponent, "the dust components")
        ATTRIBUTE_DEFAULT_VALUE(components, "MeshDustComponent")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

    /** The destructor deletes the data structures allocated during setup. */
    ~AdaptiveMeshDustDistribution();

protected:
    /** This function imports the adaptive mesh data (we need to know the number of required data
        fields, so our dust components must already have been setup). */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the dust distribution, which for this class is
        always 3 since there are no symmetries in the geometry. */
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
        at the position \f${\bf{r}}\f$. For this type of dust distribution, it just sums the
        contribution of the different components. */
    double density(Position bfr) const override;

    /** This function generates a random position from the dust distribution. It randomly chooses a
        mesh cell from the normalized cumulative density distribution constructed during the setup
        phase. Then a position is determined randomly within the cell boundaries. */
    Position generatePosition() const override;

    /** This function returns the total mass of the \f$h\f$'th component of the dust distribution. */
    double mass(int h) const override;

    /** This function returns the total dust mass of the dust distribution. For this type of dust
        distribution, it just sums the contribution of the different components. */
    double mass() const override;

    /** This function returns the X-axis surface density of the dust distribution, defined as the
        mass density integrated along the entire X-axis, \f[ \Sigma_X =
        \int_{-\infty}^\infty \rho(x,0,0)\, {\text{d}}x.\f] For this type of dust distribution, it
        just sums the contribution of the different components. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density of the dust distribution, defined as the
        mass density integrated along the entire X-axis, \f[ \Sigma_Y =
        \int_{-\infty}^\infty \rho(0,y,0)\, {\text{d}}y.\f] For this type of dust distribution, it
        just sums the contribution of the different components. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density of the dust distribution, defined as the
        mass density integrated along the entire Z-axis, \f[ \Sigma_Z =
        \int_{-\infty}^\infty \rho(0,0,z)\, {\text{d}}z.\f] For this type of dust distribution, it
        just sums the contribution of the different components. */
    double SigmaZ() const override;

    /** This function implements the AdaptiveMeshInterface interface. It returns a pointer to the
        adaptive mesh maintained by this dust distribution. */
    AdaptiveMesh* mesh() const override;

    //======================== Data Members ========================

private:
    // initialized during setup
    AdaptiveMesh* _mesh{nullptr};
    Array _cumrhov;
};

////////////////////////////////////////////////////////////////////

#endif
