/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHSTELLARCOMP_HPP
#define ADAPTIVEMESHSTELLARCOMP_HPP

#include "BoxStellarComp.hpp"
#include "ArrayTable.hpp"
#include "AdaptiveMeshFile.hpp"
class AdaptiveMesh;
class Random;

//////////////////////////////////////////////////////////////////////

/** The AdaptiveMeshStellarComp class represents a stellar component defined by the stellar density
    and properties imported from an adaptive mesh data file. The data file must have one of the
    supported formats; refer to the AdaptiveMeshFile class and its subclasses. The
    AdaptiveMeshStellarComp class allows selecting the data columns respectively containing the
    initial stellar density \f$\rho\f$ (in \f$M_\odot\,\text{pc}^{-3}\f$ at \f$t=0\f$), the
    metallicity \f$Z\f$ of the stellar population (dimensionless fraction), and the age of the
    stellar population (in yr). Since the adaptive mesh data format does not specify the size of
    the domain, this information must be provided as properties of this class as well. */
class AdaptiveMeshStellarComp : public BoxStellarComp
{
    ITEM_CONCRETE(AdaptiveMeshStellarComp, BoxStellarComp,
                  "a stellar component imported from an adaptive mesh data file")

    PROPERTY_ITEM(adaptiveMeshFile, AdaptiveMeshFile, "the adaptive mesh data file")
        ATTRIBUTE_DEFAULT_VALUE(adaptiveMeshFile, "AdaptiveMeshAsciiFile")

    PROPERTY_INT(densityIndex, "the index of the column defining the stellar density distribution")
        ATTRIBUTE_MIN_VALUE(densityIndex, "0")
        ATTRIBUTE_MAX_VALUE(densityIndex, "99")
        ATTRIBUTE_DEFAULT_VALUE(densityIndex, "0")

    PROPERTY_INT(metallicityIndex, "the index of the column defining the metallicity of the stellar population")
        ATTRIBUTE_MIN_VALUE(metallicityIndex, "0")
        ATTRIBUTE_MAX_VALUE(metallicityIndex, "99")
        ATTRIBUTE_DEFAULT_VALUE(metallicityIndex, "1")

    PROPERTY_INT(ageIndex, "the index of the column defining the age of the stellar population")
        ATTRIBUTE_MIN_VALUE(ageIndex, "0")
        ATTRIBUTE_MAX_VALUE(ageIndex, "99")
        ATTRIBUTE_DEFAULT_VALUE(ageIndex, "2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor releases the data structure allocated during setup. */
    ~AdaptiveMeshStellarComp();

protected:
    /** This function imports the adaptive mesh data and calculates the luminosity \f$L_{\ell,m}\f$
        for each mesh cell \f$m\f$ and at each wavelength grid point \f$\ell\f$. The luminosity
        distribution in each cell is determined using bilinear interpolation on age and metallicity
        in the family of Bruzual & Charlot SEDs (read from the appropriate resource files). These
        cell-specific SEDs are resampled on the simulation's global wavelength grid. Summing over
        all cells, a vector with the total luminosity for each wavelength bin is constructed.
        Finally, a matrix \f$X_{\ell,m}\f$ is filled that contains the normalized cumulative
        luminosity, \f[ X_{\ell,m} = \frac{ \sum_{i=0}^{m-1} L_{\ell,i} }{ \sum_{i=0}^{N-1}
        L_{\ell,i} } \f] where \f$N\f$ is the total number of mesh cells. This matrix will be used
        for the efficient generation of random photon packages from the stellar component. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the luminosity \f$L_\ell\f$ of the stellar component in the
        wavelength bin with index \f$\ell\f$. It just reads the appropriate number from the
        internally stored vector. */
    double luminosity(int ell) const override;

    /** This function simulates the emission of a monochromatic photon package with a luminosity
        \f$L\f$ at wavelength index \f$\ell\f$ from the stellar component. It randomly chooses a
        mesh cell from the normalized cumulative luminosity matrix defined in the setup phase and
        stored internally. Once the cell has been determined, a position is determined randomly
        within the cell boundaries. */
    void launch(PhotonPackage* pp, int ell, double L) const override;

    //======================== Data Members ========================

private:
    // other data members
    Random* _random{nullptr};
    AdaptiveMesh* _mesh{nullptr};
    Array _Ltotv;
    ArrayTable<2> _Xvv;
};

////////////////////////////////////////////////////////////////////

#endif
