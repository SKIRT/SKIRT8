/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHGEOMETRY_HPP
#define ADAPTIVEMESHGEOMETRY_HPP

#include "BoxGeometry.hpp"
#include "AdaptiveMeshInterface.hpp"
#include "AdaptiveMeshFile.hpp"
#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** The AdaptiveMeshGeometry class describes an arbitrary 3D geometry defined by the
    probability distribution imported from an adaptive mesh data file. The data file must have a
    format supported by one of the AdaptiveMeshFile subclasses. The AdaptiveMeshGeometry class
    allows selecting the data column containing the probability distribution, and it offers the
    option to use a second column as a multiplication factor (i.e. the probability distribution is
    then defined by the product of the two column values). The geometry will be normalized anyway
    after importing the probability distribution, so the probability distribution in the data file
    does not have to be normalized, and the units of the values in the data file are irrelevant.
    Since the adaptive mesh data format does not specify the size of the domain, this information
    must be provided as properties of this class as well. */
class AdaptiveMeshGeometry : public BoxGeometry, public AdaptiveMeshInterface
{
    ITEM_CONCRETE(AdaptiveMeshGeometry, BoxGeometry, "a geometry imported from an adaptive mesh data file")

    PROPERTY_ITEM(adaptiveMeshFile, AdaptiveMeshFile, "the adaptive mesh data file")
        ATTRIBUTE_DEFAULT_VALUE(adaptiveMeshFile, "AdaptiveMeshAsciiFile")

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
    ~AdaptiveMeshGeometry();

protected:
    /** This function imports the adaptive mesh data and calculates the total density integrated
        over the complete domain. */
    void setupSelfBefore() override;

    //======== Setters & Getters for Discoverable Attributes =======

    /** \fn densityIndex
        The geometry will be normalized after importing the probability distribution, so the
        probability distribution in the data file does not have to be normalized, and the units of
        the values in the data file are irrelevant. */

    /** \fn multiplierIndex
        If this column index is nonnegative, the probability distribution is effectively defined by
        the product of the two column values specified by \em densityIndex and \em multiplierIndex.
        */

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

    /** This function implements the AdaptiveMeshInterface interface. It returns a
        pointer to the adaptive mesh maintained by this geometry. */
    AdaptiveMesh* mesh() const override;

    //======================== Data Members ========================

private:
    // other data members
    AdaptiveMesh* _mesh{nullptr};
    Array _cumrhov;
};

////////////////////////////////////////////////////////////////////

#endif
