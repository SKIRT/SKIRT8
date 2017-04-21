/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STELLARSYSTEM_HPP
#define STELLARSYSTEM_HPP

#include "SimulationItem.hpp"
#include "ArrayTable.hpp"
#include "StellarComp.hpp"
class PhotonPackage;
class Random;

//////////////////////////////////////////////////////////////////////

/** An instance of the StellarSystem class represents a complete stellar system, which is the
    superposition of one or more stellar components. Each stellar component provides a complete
    description of the spatial and spectral distribution of the stars (or any other primary source
    of radiation, such as an AGN). The main function of this class is to manage a list of objects
    of type StellarComp, and to implement the superposition of the distributions defined in these
    objects. */
class StellarSystem : public SimulationItem
{
    ITEM_CONCRETE(StellarSystem, SimulationItem, "a stellar system")

    PROPERTY_ITEM_LIST(components, StellarComp, "the stellar components")
        ATTRIBUTE_DEFAULT_VALUE(components, "GeometricStellarComp")

    PROPERTY_DOUBLE(emissionBias, "the stellar emission bias")
        ATTRIBUTE_MIN_VALUE(emissionBias, "[0")
        ATTRIBUTE_MAX_VALUE(emissionBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(emissionBias, "0.5")
        ATTRIBUTE_SILENT(emissionBias)

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches some pointers. */
    void setupSelfBefore() override;

    /** This function calculates and cashes luminosity information about the components for later
        use. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the stellar system, which depends on the (lack of)
        symmetry in the geometries of its stellar components. A value of 1 means spherical
        symmetry, 2 means axial symmetry and 3 means none of these symmetries. The stellar
        component with the least symmetry (i.e. the highest dimension) determines the result for
        the whole system. */
    int dimension() const;

    /** This function returns the number of components in the stellar system. */
    int numComponents() const;

    /** This function returns the monochromatic luminosity \f$L_\ell\f$ of the stellar system at
        the wavelength index \f$\ell\f$, which is the sum of the luminosities of the stellar
        components in the system. */
    double luminosity(int ell) const;

    /** This function simulates the emission of a monochromatic photon package with a monochromatic
        luminosity \f$L\f$ at wavelength index \f$\ell\f$ from the stellar system. It randomly
        chooses a stellar component from which to emit the photon and then simulates the emission
        through the corresponding StellarComp::launch() function. */
    void launch(PhotonPackage* pp, int ell, double L) const;

    //======================== Data Members ========================

private:
    Array _Lv;
    ArrayTable<2> _Xvv;
    Random* _random{nullptr};
};

////////////////////////////////////////////////////////////////

#endif
