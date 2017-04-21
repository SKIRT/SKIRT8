/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef POINTGEOMETRY_HPP
#define POINTGEOMETRY_HPP

#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

/** The PointGeometry class is a subclass of the Geometry class, and
    describes the geometry of a single point positioned in the centre of the coordinate
    system. It's rather trivial as it has no data members. */
class PointGeometry : public Geometry
{
    ITEM_CONCRETE(PointGeometry, Geometry, "a point source geometry")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the geometry, which is 0 in this case. */
    int dimension() const override;

    /** This function returns the stellar density \f$\rho({\bf{r}})\f$ at the position
        \f${\bf{r}}\f$. In the case of a point geometry, the density takes the form of a
        Dirac delta function, \f$\rho({\bf{r}}) = \delta({\bf{r}})\f$. The function returns
        infinity if \f${\bf{r}} = {\bf{0}}\f$ and zero in all other cases. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry, by
        drawing a random point from the three-dimensional probability density \f$p({\bf{r}})\,
        {\text{d}}{\bf{r}} = \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. In the case of a point
        geometry, it always returns the origin of the coordinate system. */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density, i.e. the integration of
        the density along the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,
        {\text{d}}x. \f] For the point geometry, this integral is infinity. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of
        the density along the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\,
        {\text{d}}y. \f] For the point geometry, this integral is infinity. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of
        the density along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,
        {\text{d}}z. \f] For the point geometry, this integral is infinity. */
    double SigmaZ() const override;
};

////////////////////////////////////////////////////////////////////

#endif
