/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GEOMETRICSTELLARCOMP_HPP
#define GEOMETRICSTELLARCOMP_HPP

#include "StellarComp.hpp"
#include "Array.hpp"
#include "Geometry.hpp"
class PhotonPackage;

//////////////////////////////////////////////////////////////////////

/** GeometricStellarComp represents a stellar component in which the spatial distribution of stars
    is characterized by a built-in geometry. This abstract class handles an instance of the
    Geometry class to define the spatial distribution. Each subclasses is expected to define the
    spectral energy distribution of the stars (which is constant across the spatial distribution)
    and some form of normalization to specify the total luminosity of the component. */
class GeometricStellarComp : public StellarComp
{
    ITEM_ABSTRACT(GeometricStellarComp, StellarComp, "a stellar component with a built-in geometry")

    PROPERTY_ITEM(geometry, Geometry, "the geometry of the spatial stellar distribution")
        ATTRIBUTE_DEFAULT_VALUE(geometry, "PointGeometry")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the stellar component, which depends on the (lack of)
        symmetry of its geometry. A value of 1 means spherical symmetry, 2 means axial symmetry and
        3 means none of these symmetries. */
    int dimension() const override;

    /** This function simulates the emission of a monochromatic photon package with a luminosity
        \f$L\f$ at wavelength index \f$\ell\f$ from the stellar component. The position and
        propagation direction of the emission are determined randomly from the geometry of the
        stellar component. */
    void launch(PhotonPackage* pp, int ell, double L) const override;
};

//////////////////////////////////////////////////////////////////////

#endif
