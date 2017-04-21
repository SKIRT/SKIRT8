/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MAPPINGSSEDFAMILY_HPP
#define MAPPINGSSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "ArrayTable.hpp"
class WavelengthGrid;

//////////////////////////////////////////////////////////////////////

/** An instance of the MappingsSEDFamily class represents the family of MAPPINGS III
    starburst template SEDs, parameterized on metallicity, compactness, ISM pressure and PDR
    covering factor, as described in Groves et al. (2008) ApJS,176,438. The data was downloaded
    from http://www.mpia-hd.mpg.de/~brent/starburst.html -> Cparam_models.save and converted to
    plain text files using a simple IDL script. The MAPPINGS III library data is read from the
    appropriate resource files during setup, and it is subsequently interpolated to the
    desired parameters and wavelength grid points by calling the luminosities() function as often
    as needed. */
class MappingsSEDFamily : public SEDFamily
{
    ITEM_CONCRETE(MappingsSEDFamily, SEDFamily, "a MAPPINGS III SED family for starbursting regions")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit MappingsSEDFamily(SimulationItem* parent);

protected:
    /** This function reads the MAPPINGS III library data from the appropriate resource files and
        stores all relevant information internally. */
    void setupSelfBefore() override;

    //====================== Retrieving an SED =====================

public:
    /** This function returns the luminosity \f$L_\ell\f$ at each wavelength in the simulation's
        wavelength grid for a starbursting population, given the star formation rate \f$\dot{M}\f$
        (assumed to be constant over the past 10 Myr, in \f$M_\odot\,{\text{yr}}^{-1}\f$),
        metallicity \f$Z\f$, the logarithm of the compactness \f$\log C\f$, the ISM pressure
        \f$p\f$ (in Pa), and the dimensionless PDR covering factor \f$f_{\text{PDR}}\f$.

        If the redshift argument \f$z\f$ is present and nonzero, the spectrum is redshifted
        according to the specified value. In function of the rest wavelength \f$\lambda_0\f$, the
        redshifted wavelength \f$\lambda_z\f$ is given by \f$\lambda_z=(1+z)\,\lambda_0\f$. For
        \f$z<<1\f$, the observed luminosity at each wavelength in the simulation's wavelength grid
        can then be written as \f[L_z[\lambda_\ell] = L_0[(1-z)\,\lambda_\ell]\f] */
    Array luminosities(double SFR, double Z, double logC, double pressure, double fPDR, double z=0) const;

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

    //====================== Data members =====================

private:
    WavelengthGrid* _lambdagrid{nullptr};

    // contents of the library, read by constructor
    Array _lambdav;
    Array _Zrelv;
    Array _logCv;
    Array _logpv;
    ArrayTable<4> _j0vv;
    ArrayTable<4> _j1vv;
};

////////////////////////////////////////////////////////////////////

#endif
