/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FSPSVARIABLEIMFSEDFAMILY_HPP
#define FSPSVARIABLEIMFSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "ArrayTable.hpp"
#include <unordered_map>
class WavelengthGrid;

//////////////////////////////////////////////////////////////////////

/** An instance of the FSPSVariableIMFSEDFamily class represents one of two SED families for single
    stellar populations with a variable initial mass function (IMF), parameterized on metallicity
    and age in addition to the IMF shape. For the first family ("VaryLowMass"), the IMF slope
    below 0.5 Msun is varied while the high mass slope is fixed at the Kroupa value. For the second
    family ("VaryHighMass"), the slope below 0.5 Msun is fixed at the Kroupa value and the
    high mass slope is varied instead.

    The SED templates were prepared by Chris Barber (Leiden Obervatory, The Netherlands) using the
    FSPS code (Conroy, Gunn, & White 2009, ApJ, 699, 486; Conroy & Gunn 2010, ApJ, 712, 833)

    The templates are read from the appropriate resource files during setup, and are subsequently
    interpolated to the desired parameters and wavelength grid points by calling the luminosities()
    function as often as needed. */
class FSPSVariableIMFSEDFamily : public SEDFamily
{
    /** The enumeration type indicating the type of initial mass function (IMF). */
    ENUM_DEF(InitialMassFunction, VaryLowMass, VaryHighMass)
    ENUM_VAL(InitialMassFunction, VaryLowMass, "IMF with variable slope below 0.5 Msun")
    ENUM_VAL(InitialMassFunction, VaryHighMass, "IMF with variable slope above 0.5 Msun")
    ENUM_END()

    ITEM_CONCRETE(FSPSVariableIMFSEDFamily, SEDFamily, "a variable-IMF SED family for single stellar populations")

    PROPERTY_ENUM(initialMassFunction, InitialMassFunction, "the type of initial mass function")
        ATTRIBUTE_DEFAULT_VALUE(initialMassFunction, "VaryLowMass")

     ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit FSPSVariableIMFSEDFamily(SimulationItem* parent);

protected:
    /** This function reads the grid points in each dimension of the parameter space from the
        appropriate resource files and initializes the data structure where the SED templates will
        be cached after being loaded by the luminosities() function. */
    void setupSelfBefore() override;

    //====================== Retrieving an SED =====================

public:
    /** This function returns the luminosity \f$L_\ell\f$ at each wavelength in the simulation's
        wavelength grid for a stellar population with given initial mass \em M (in \f$M_\odot\f$ at
        \f$t=0\f$), metallicity \em Z (as a dimensionless fraction), age \em t (in years), and IMF
        slope \f$\alpha\f$ (a dimensionless exponent). The luminosity is defined as the emissivity
        multiplied by the width of the wavelength bin.

        If the redshift argument \f$z\f$ is present and nonzero, the spectrum is redshifted
        according to the specified value. In function of the rest wavelength \f$\lambda_0\f$, the
        redshifted wavelength \f$\lambda_z\f$ is given by \f$\lambda_z=(1+z)\,\lambda_0\f$. For
        \f$z<<1\f$, the observed luminosity at each wavelength in the simulation's wavelength grid
        can then be written as \f[L_z[\lambda_\ell] = L_0[(1-z)\,\lambda_\ell]\f] */
    Array luminosities(double M, double Z, double t, double alpha, double z=0) const;

    /** This function returns the number of parameters used by this particular %SED family, in
        other words the number of arguments of the luminosities() function. */
    int numParams() const override;

    /** This function returns the luminosity \f$L_\ell\f$ at each wavelength in the simulation's
        wavelength grid for the specified parameter values, which should be in the same order and
        using the same units as the arguments described for the luminosities() function. The first
        \em skipvals values in the \em params array are ignored. */
    Array luminositiesGeneric(const Array& params, int skipvals=0, double z=0) const override;

    /** This function returns the mass (in \f$M_\odot\f$) of the source represented by the
        specified set of parameter values. The \em params array must contain the appropriate number
        of parameter values in the order specified by the particular %SED family subclass. The
        first \em skipvals values in the array are ignored. */
    double massGeneric(const Array& params, int skipvals=0) const override;

    /** This function returns a short name for the type of sources typically assigned to this
        particular %SED family. */
     string sourceName() const override;

     /** This function returns a description for the type of sources typically assigned to this
         particular %SED family. */
      string sourceDescription() const override;

private:
      /** This function returns (a reference to) the emissivity table corresponding to the
          specified indices for IMF slope and metallicity. The function loads the table from the
          appropriate resource file if needed; if the table has already been loaded, it is returned
          from the cache. */
      const ArrayTable<2>& getEmissivityTable(size_t s, size_t m) const;

    //====================== Data members =====================

private:
    WavelengthGrid* _lambdagrid{nullptr};

    // coordinate grid points along each dimension, loaded by constructor
    Array _alphav;      // index s
    Array _Zv;          // index m
    Array _tv;          // index p
    Array _lambdav;     // index k

    // emissivity tables, loaded on demand during operation
    mutable std::unordered_map<size_t, ArrayTable<2>> _jmap;  // key(s,m); indices p,k
};

////////////////////////////////////////////////////////////////////

#endif
