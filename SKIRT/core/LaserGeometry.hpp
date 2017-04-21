/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LASERGEOMETRY_HPP
#define LASERGEOMETRY_HPP

#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

/** The LaserGeometry class is a subclass of the Geometry class. It represents a
    point source that emits all its radiation towards the positive Z-axis. As a
    SKIRT Geometry class, this is a fairly simple geometry, with the main peculiarity
    that the emission is anisotropic. The emission pattern is axisymmetric, so this geometry
    has a dimension of 2. */
class LaserGeometry : public Geometry
{
    ITEM_CONCRETE(LaserGeometry, Geometry, "a laser geometry")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the geometry, which is 2 in this
        case. */
    int dimension() const override;

    /** This function returns the density \f$\rho({\bf{r}})\f$ at the position
        \f${\bf{r}}\f$. For this geometry, the density takes the form of a
        Dirac delta function, \f$\rho({\bf{r}}) = \delta({\bf{r}})\f$. The function
        returns infinity if \f${\bf{r}} = {\bf{0}}\f$ and zero in all other cases. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry, by
        drawing a random point from the three-dimensional probability density
        \f$p({\bf{r}})\,{\text{d}}{\bf{r}} = \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. For
        this geometry, it always returns the origin of the coordinate system. */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density, i.e. the integration of the
        density along the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,
        {\text{d}}x. \f] For this geometry, this integral is infinity. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the
        density along the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\,
        {\text{d}}y. \f] For this geometry, this integral is infinity. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the
        density along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,
        {\text{d}}z. \f] For this geometry, this integral is infinity. */
    double SigmaZ() const override;

    // - - - - - - AngularDistribution interface - - - - - -

    /** This function returns the normalized probability for a given direction \f${\bf{k}}
        = (\theta,\phi)\f$ for radiation emitted at the position \f${\bf{r}}\f$. When
        \f${\bf{r}}\f$ is not the origin, this probability function is not defined and a fatal
        error is returned. When \f${\bf{r}}\f$ is the origin, the appropriate probability is a
        Dirac delta function, \f$4\pi\delta({\bf{k}}-{\bf{e}}_z)\f$. The function returns
        infinity if \f${\bf{k}} = {\bf{e}}_z\f$, or equivalently if \f$\theta=0\f$, and
        zero in all other cases. */
    double probabilityForDirection(int ell, Position bfr, Direction bfk) const override;

    /** This function generates a random direction drawn from the appropriate angular
        probability distribution at the specified position. When \f${\bf{r}}\f$ is not the
        origin, a fatal error is returned. When \f${\bf{r}}\f$ is the origin, this function
        simply returns \f${\bf{k}} = {\bf{e}}_z\f$. */
    Direction generateDirection(int ell, Position bfr) const override;
};

////////////////////////////////////////////////////////////////////

#endif
