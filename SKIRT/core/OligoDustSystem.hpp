/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OLIGODUSTSYSTEM_HPP
#define OLIGODUSTSYSTEM_HPP

#include "DustSystem.hpp"

//////////////////////////////////////////////////////////////////////

/** An OligoDustSystem class object represents a complete dust system for use with oligochromatic
    simulations. All functionality is implemented in the DustSystem base class. */
class OligoDustSystem : public DustSystem
{
    ITEM_CONCRETE(OligoDustSystem, DustSystem, "a dust system for use with oligochromatic simulations")

    PROPERTY_BOOL(writeMeanIntensity, "output FITS files displaying the mean intensity of the radiation field")
        ATTRIBUTE_DEFAULT_VALUE(writeMeanIntensity, "false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function  resizes the absorption rate matrix, if required. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function always returns false since oligochromatic simulations do no support dust
        emission. */
    bool hasDustEmission() const override;

    /** This function returns a flag that indicate whether the absorption rates in each cell need
        to be stored for this dust system. For a panchromatic simulation, absorption rates are only
        calculated if dust emission is turned on. */
    bool hasDustAbsorption() const override;

    /** The function simulates the absorption of a monochromatic luminosity package in the
        specified dust cell. For an oligochromatic simulation, it should only be invoked if the
        absorption rates are explicitly requested, and the ynstellar flag should always be true
        (there is only stellar emission in oligochromatic simulations). If these conditions are
        met, the routine just adds the absorbed luminosity \f$\Delta L\f$ to the appropriate item
        in the absorption rate table. The addition is performed in a thread-safe manner so this
        function may be concurrently called from multiple threads. */
    void absorb(int m, int ell, double DeltaL, bool ynstellar) override;

    /** This function returns the absorbed luminosity \f$L_{\ell,m}\f$ at wavelength index
        \f$\ell\f$ in the dust cell with cell number \f$m\f$. For an oligochromatic dust system, it
        simply reads the corresponding absorption rate counter. */
    double absorbedLuminosity(int m, int ell) const override;

    /** If the writeMeanIntensity attribute is true, this function writes out FITS files (named
        <tt>prefix_ds_JXX.fits</tt>) with the mean radiation field in the coordinate planes. Each
        of these maps contains 1024 x 1024 pixels, and covers as a field of view the total
        extension of the grid. The number of data files written depends on the dimension of the
        dust system's geometry: for spherical symmetry only the intersection with the xy plane is
        written, for axial symmetry the intersections with the xy and xz planes are written, and
        for general geometries all three intersections are written. Each FITS file is a data cube,
        where a map is created for each wavelength in the global wavelength grid. */
    void write() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Table<2> _Labsvv;   // absorbed emission for each cell and each wavelength (indexed on m,ell)
};

//////////////////////////////////////////////////////////////////////

#endif
