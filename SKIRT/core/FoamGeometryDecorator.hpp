/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FOAMGEOMETRYDECORATOR_HPP
#define FOAMGEOMETRYDECORATOR_HPP

#include "BoxGeometry.hpp"
#include "FoamDensity.hpp"
class Foam;

////////////////////////////////////////////////////////////////////

/** The FoamGeometryDecorator class is a decorator for the Geometry class that provides an
    alternative method to generate random positions, using a three-dimensional cell structure,
    called the foam. A foam is based on the three-dimensional unit cube \f$[0,1]^3\f$, subdivided
    into a large number of small cuboidal cells. The distribution of the grid cells is performed
    completely automatically, based on the density distribution of the geometry that is being
    decorated. The foam implementation, characterized by the Foam class, allows to efficiently
    generate random numbers drawn from this probability distribution on this unit cube. One problem
    is that the stellar density \f$\rho({\bf{r}})\f$ is typically defined on the entire 3D space,
    whereas the foam requires a density distribution on the unit cube. We solve this problem using
    a simple linear transformation, where we map the volume from which we sample (assumed to be a
    box) to the unit cube. */
class FoamGeometryDecorator : public BoxGeometry, FoamDensity
{
    ITEM_CONCRETE(FoamGeometryDecorator, BoxGeometry,
                  "a decorator that provides an alternative random position generator")

    PROPERTY_ITEM(geometry, Geometry, "the geometry to be alternatively sampled")

    PROPERTY_INT(numCells, "the number of cells in the foam")
        ATTRIBUTE_MIN_VALUE(numCells, "1000")
        ATTRIBUTE_MAX_VALUE(numCells, "1000000")
        ATTRIBUTE_DEFAULT_VALUE(numCells, "10000")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor; it deallocates the foam object created during setup. */
    ~FoamGeometryDecorator();

protected:
    /** This function sets up the foam. The Foam constructor needs a pointer to an instance of the
        FoamDensity interface, which is implemented for this purpose by the FoamGeometryDecorator
        class. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho({\bf{r}})\f$ at the position
        \f${\bf{r}}\f$. It simply calls the density() function for the geometry being
        decorated. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry, by drawing a random
        point from the three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}
        {\bf{r}} = \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. This task is accomplished by
        drawing a random point \f$(\bar{\bf{r}})\f$ from the foam and converting this to
        a position \f${\bf{r}}\f$ using the appropriate transformation. */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density, i.e. the integration of the
        density along the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)
        \,{\text{d}}x. \f] It just calls the corresponding function for the geometry being
        decorated. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the
        density along the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)
        \,{\text{d}}y. \f] It just calls the corresponding function for the geometry being
        decorated. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the
        density along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)
        \,{\text{d}}z. \f] It just calls the corresponding function for the geometry being
        decorated. */
    double SigmaZ() const override;

    /** This function implements the FoamDensity interface to convert the stellar density
        distribution \f$\rho({\bf{r}})\f$ to a new probability density distribution
        \f$\bar\rho(\bar{\bf{r}})\f$ on the three-dimensional unit cube \f$[0,1]^3\f$ such
        that it can be used as input for the construction of a foam. The function returns the
        density \f$\bar\rho(\bar{\bf{r}})\f$ in a point \f$\bar{\bf{r}}\f$ in the unit cube,
        defined by the input parameters (\f$n_{\text{par}}\f$ is always equal to 3 and the
        pointer contains the three cartesian parameters \f$(\bar{x},\bar{y},\bar{z})\f$ that
        define the position \f$\bar{\bf{r}}\f$). The required density \f$\bar\rho(
        \bar{\bf{r}})\f$ is determined as \f[ \bar\rho(\bar{x},\bar{y},\bar{z}) = \rho
        \Bigl( x(\bar{x}), y(\bar{y}), z(\bar{z}) \Bigr)\, \frac{{\text{d}}x}{{\text{d}}
        \bar{x}}\, \frac{{\text{d}}y}{{\text{d}}\bar{y}}\, \frac{{\text{d}}z}{{\text{d}}
        \bar{z}}, \f] with \f$\rho(x,y,z)\f$ the stellar density. The transformation between
        the coordinates \f$(x,y,z)\f$ and \f$(\bar{x},\bar{y},\bar{z})\f$ is a simple linear
        transformation from the cuboid defined by \f$x_{\text{min}} < x < x_{\text{max}}\f$,
        \f$y_{\text{min}} < y < y_{\text{max}}\f$ and \f$z_{\text{min}} < z <
        z_{\text{max}}\f$ to the unit cube. */
    double foamDensity(int ndim, double* par) const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    double _jacobian{0.};
    Foam* _foam{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
