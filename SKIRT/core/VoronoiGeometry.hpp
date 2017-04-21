/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VORONOIGEOMETRY_HPP
#define VORONOIGEOMETRY_HPP

#include "BoxGeometry.hpp"
#include "DustParticleInterface.hpp"
#include "VoronoiMeshInterface.hpp"
#include "VoronoiMeshFile.hpp"
#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** The VoronoiGeometry class describes an arbitrary 3D geometry defined by the
    probability distribution imported from a Voronoi mesh data file. The data file must have a
    format supported by one of the VoronoiMeshFile subclasses. The VoronoiGeometry class
    allows selecting the data column containing the probability distribution, and it offers the
    option to use a second column as a multiplication factor (i.e. the probability distribution is
    then defined by the product of the two column values). The geometry will be normalized anyway
    after importing the probability distribution, so the probability distribution in the data file
    does not have to be normalized, and the units of the values in the data file are irrelevant.
    Since the Voronoi mesh data format does not specify the size of the domain, this information
    must be provided as properties of this class as well. */
class VoronoiGeometry : public BoxGeometry, public VoronoiMeshInterface, public DustParticleInterface
{
    ITEM_CONCRETE(VoronoiGeometry, BoxGeometry, "a geometry imported from a Voronoi mesh data file")

    PROPERTY_ITEM(voronoiMeshFile, VoronoiMeshFile, "the Voronoi mesh data file")
        ATTRIBUTE_DEFAULT_VALUE(voronoiMeshFile, "VoronoiMeshAsciiFile")

    PROPERTY_INT(densityIndex, "the index of the column defining the density distribution")
        ATTRIBUTE_MIN_VALUE(densityIndex, "0")
        ATTRIBUTE_MAX_VALUE(densityIndex, "99")
        ATTRIBUTE_DEFAULT_VALUE(densityIndex, "0")

    PROPERTY_INT(multiplierIndex, "the index of the column defining an extra multiplication factor, or -1")
        ATTRIBUTE_MIN_VALUE(multiplierIndex, "-1")
        ATTRIBUTE_MAX_VALUE(multiplierIndex, "99")
        ATTRIBUTE_DEFAULT_VALUE(multiplierIndex, "-1")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor releases the data structure allocated during setup. */
    ~VoronoiGeometry();

protected:
    /** This function verifies the property values, imports the Voronoi mesh data, and calculates
        the total density integrated over the complete domain. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho({\bf{r}})\f$ for this geometry at the
        position \f${\bf{r}}\f$. It forwards the call to the mesh, and applies the necessary
        normalization (the density on an adaptive mesh is not necessarily normalized to one). */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry, by
        drawing a random point from the three-dimensional probability density \f$p({\bf{r}})\,
        {\text{d}}{\bf{r}} = \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density of the geometry, defined as the integration
        of the density along the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,
        {\text{d}}x.\f] It forwards the call to the mesh, and applies the necessary
        normalization (the density on an adaptive mesh is not necessarily normalized to one). */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density of the geometry, defined as the integration
        of the density along the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\,
        {\text{d}}y.\f] It forwards the call to the mesh, and applies the necessary
        normalization (the density on an adaptive mesh is not necessarily normalized to one). */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density of the geometry, defined as the integration
        of the density along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,
        {\text{d}}z.\f] It forwards the call to the mesh, and applies the necessary
        normalization (the density on an adaptive mesh is not necessarily normalized to one). */
    double SigmaZ() const override;

    /** This function implements the VoronoiMeshInterface interface. It returns a
        pointer to the Voronoi mesh maintained by this geometry. */
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
    // data members initialized during setup
    VoronoiMesh* _mesh{nullptr};
    Array _cumrhov;
};

////////////////////////////////////////////////////////////////////

#endif
