/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VORONOIDUSTDISTRIBUTION_HPP
#define VORONOIDUSTDISTRIBUTION_HPP

#include "BoxDustDistribution.hpp"
#include "DustParticleInterface.hpp"
#include "VoronoiMeshInterface.hpp"
#include "Array.hpp"
#include "MeshDustComponent.hpp"
#include "VoronoiMeshFile.hpp"

////////////////////////////////////////////////////////////////////

/** The VoronoiDustDistribution class represents a dust distribution imported from a Voronoi mesh
    data file. The data file must have one of the supported formats; refer to the VoronoiMeshFile
    class and its subclasses. Since the Voronoi mesh data format does not specify the size of the
    domain, this information must be provided as properties of this class. This class supports
    multiple dust components, as long as the dust density distributions for all components are
    defined on the same mesh in the same Voronoi mesh data file. Each dust component is represented
    by an instance of the MeshDustComponent class, which specifies the data column index defining
    the dust density distribution for the component and the corresponding dust mix. */
class VoronoiDustDistribution : public BoxDustDistribution, public VoronoiMeshInterface, public DustParticleInterface
{
    ITEM_CONCRETE(VoronoiDustDistribution, BoxDustDistribution,
                  "a dust distribution imported from a Voronoi mesh data file")

    PROPERTY_ITEM(voronoiMeshFile, VoronoiMeshFile, "the Voronoi mesh data file")
        ATTRIBUTE_DEFAULT_VALUE(voronoiMeshFile, "VoronoiMeshAsciiFile")

    PROPERTY_DOUBLE(densityUnits, "the units in which the file specifies density values")
        ATTRIBUTE_QUANTITY(densityUnits, "massvolumedensity")
        ATTRIBUTE_MIN_VALUE(densityUnits, "]0")
        ATTRIBUTE_DEFAULT_VALUE(densityUnits, "1 Msun/pc3")

    PROPERTY_ITEM_LIST(components, MeshDustComponent, "the dust components")
        ATTRIBUTE_DEFAULT_VALUE(components, "MeshDustComponent")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor deletes the data structures allocated during setup. */
    ~VoronoiDustDistribution();

protected:
    /** This function imports the Voronoi mesh data (we need to know the number of required data
        fields, so our dust components must already have been setup). */
    virtual void setupSelfAfter() override;

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

    /** This function returns the total mass of the \f$h\f$'th component. */
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
        mass density integrated along the entire Y-axis, \f[ \Sigma_Y =
        \int_{-\infty}^\infty \rho(0,y,0)\, {\text{d}}y.\f] For this type of dust distribution, it
        just sums the contribution of the different components. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density of the dust distribution, defined as the
        mass density integrated along the entire Z-axis, \f[ \Sigma_Z =
        \int_{-\infty}^\infty \rho(0,0,z)\, {\text{d}}z.\f] For this type of dust distribution, it
        just sums the contribution of the different components. */
    double SigmaZ() const override;

    /** This function implements the VoronoiMeshInterface interface. It returns a pointer to the
        Voronoi mesh maintained by this dust distribution. */
    VoronoiMesh* mesh() const override;

    /** This function implements (part of) the DustParticleInterface interface. It returns the
        number of particles defining the dust distribution. */
    int numParticles() const override;

    /** This function implements (part of) the DustParticleInterface interface. It returns the
        coordinates of the dust-distribution-defining particle with the specified zero-based index.
        If the index is out of range, a fatal error is thrown. */
    Vec particleCenter(int index) const override;

    //======================== Data Members ========================

private:
    // initialized during setup
    VoronoiMesh* _mesh{nullptr};
    Array _cumrhov;
};

////////////////////////////////////////////////////////////////////

#endif
