/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONICALPOINTGEOMETRY_HPP
#define CONICALPOINTGEOMETRY_HPP

#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

/** The ConicalPointGeometry class is a subclass of the Geometry class. It represents a single
    point positioned in the centre of the coordinate system, with an anisotropic emission that can
    represent a simple model for an AGN: the emission is isotropic within a cone around the
    symmetry axis (on both sides), and completely blocked outside this cone. The (half) opening
    angle \f$\Delta\f$ of the cone is the only free parameter. The emission pattern is
    axisymmetric, so this geometry has a dimension of 2. */
class ConicalPointGeometry : public Geometry
{
    ITEM_CONCRETE(ConicalPointGeometry, Geometry,
                  "a point-like geometry with anisotropic conical emission")

    PROPERTY_DOUBLE(openingAngle, "the opening angle of the cone")
    ATTRIBUTE_QUANTITY(openingAngle, "posangle")
    ATTRIBUTE_MIN_VALUE(openingAngle, "[0 deg")
    ATTRIBUTE_MAX_VALUE(openingAngle, "90 deg]")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates some frequently used values. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the geometry, which is 2 in this case. */
    int dimension() const override;

    /** This function returns the density \f$\rho({\bf{r}})\f$ at the position \f${\bf{r}}\f$. For
        this geometry, the density takes the form of a Dirac delta function, \f$\rho({\bf{r}}) =
        \delta({\bf{r}})\f$. The function returns infinity if \f${\bf{r}} = {\bf{0}}\f$ and zero in
        all other cases. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry, by drawing a random point from
        the three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}{\bf{r}} =
        \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. For this geometry, it always returns the origin of
        the coordinate system. */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density, i.e. the integration of the density along
        the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\, {\text{d}}x. \f] For
        this geometry, this integral is infinity. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the density along
        the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\, {\text{d}}y. \f] For
        this geometry, this integral is infinity. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\, {\text{d}}z. \f] For
        this geometry, this integral is infinity. */
    double SigmaZ() const override;

    // - - - - - - AngularDistribution interface - - - - - -

    /** This function returns the normalized probability \f$p({\bf{k}})\f$ for a given direction
        \f${\bf{k}} = (\theta,\phi)\f$ for radiation emitted at the position \f${\bf{r}}\f$. When
        \f${\bf{r}}\f$ is not the origin, this probability function is not defined and a fatal
        error is returned. When \f${\bf{r}}\f$ is the origin, this function returns the normalized
        probability for a given direction \f${\bf{k}} = (\theta,\phi)\f$, subject to the
        normalization \f[ \frac{1}{4\pi} \iint p({\bf{k}})\, {\text{d}}\Omega = 1 \f] For an
        isotropic emission within a cone with (half) opening angle \f$\Delta\f$, we have \f[
        p({\bf{k}}) = \frac{1}{1-\cos\Delta}\times\begin{cases}\; 1 & \qquad\text{if
        $0\le\theta<\Delta$} \\ \;0 & \qquad\text{if $\Delta\le\theta<\pi-\Delta$} \\ \;1 &
        \qquad\text{if $\pi-\Delta\le\theta\le\pi$} \end{cases} \f] */
    double probabilityForDirection(int ell, Position bfr, Direction bfk) const override;

    /** This function generates a random direction drawn from the appropriate angular probability
        distribution at the specified position. When \f${\bf{r}}\f$ is not the origin, a fatal
        error is returned. When \f${\bf{r}}\f$ is the origin, this function generates a random
        direction \f${\bf{k}} = (\phi,\theta)\f$ according to an isotropic distribution within a
        cone with (half) opening angle \f$Delta\f$. This can easily be calculated: \f$\phi\f$
        distributed uniformly over the interval \f$0\le\phi\le 2\pi\f$, and \f$\theta\f$ determined
        by picking a uniform deviate \f${\cal{X}}\f$ and setting \f[ \theta = \begin{cases}
        \;\arccos\left[ 1-2{\cal{X}}(1-\cos\Delta)\right] & \qquad\text{if $0\le {\cal{X}} <
        \dfrac{\pi}{2}$} \\ \;\arccos\left[ 1-2\cos\Delta - 2{\cal{X}}(1-\cos\Delta) \right] &
        \qquad\text{if $\dfrac12 \le {\cal{X}} < 1$} \end{cases} \f] */
    Direction generateDirection(int ell, Position bfr) const override;

    //======================== Data Members ========================

private:
    // alias to discoverable data member for ease of notation and backwards compatibility
    const double& _Delta{_openingAngle};

    // data members initialized during setup
    double _cosDelta;
};

////////////////////////////////////////////////////////////////////

#endif
