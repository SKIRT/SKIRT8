/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHSTELLARCOMP_HPP
#define SPHSTELLARCOMP_HPP

#include "StellarComp.hpp"
#include "ArrayTable.hpp"
#include "SEDFamily.hpp"
#include "Vec.hpp"
class Random;
namespace SPHStellarComp_Private { class VelocityAnisotropy; }

//////////////////////////////////////////////////////////////////////

/** The SPHStellarComp class represents a stellar component defined from a set of SPH source
    particles, such as for example resulting from a cosmological simulation. The information on the
    SPH source particles is read from a file formatted as described below.

    The input text file should contain columns of numbers separated by whitespace; lines starting
    with # are ignored. The first three columns are the \f$x\f$, \f$y\f$ and \f$z\f$ coordinates of
    the particles (in pc) and the fourth column is the SPH smoothing length \f$h\f$ (in pc). If the
    \em velocity flag is set to true, the subsequent three columns (i.e. the fifth through the
    seventh columns) contain the \f$v_\mathrm{x}\f$, \f$v_\mathrm{y}\f$ and \f$v_\mathrm{z}\f$
    components of the particle's velocity (in km/s). If the \em velocity flag is set to false,
    these columns are not present.

    The number and interpretation of the subsequent columns depends on the %SED family specified as
    a property to this particular SPHStellarComp instance.

    For the Bruzual-Charlot %SED family, the remaining three columns provide the properties of the
    stellar population. The fifth column is the initial mass of the stellar population (in
    \f$M_\odot\f$ at \f$t=0\f$), the sixth column is the metallicity \f$Z\f$ of the stellar
    population (dimensionless fraction), and the seventh column is the age of the stellar
    population (in yr).

    For the Mappings III %SED family, the remaining five columns provide the properties of the
    corresponding HII region in the following order: the star formation rate \f$\dot{M}\f$, assumed
    to be constant over the past 10 Myr (in \f$M_\odot\,{\text{yr}}^{-1}\f$), metallicity \f$Z\f$
    (as a dimensionless fraction), the logarithm of the compactness \f$\log C\f$ (as a
    dimensionless fraction), the ISM pressure \f$p\f$ (in Pa), and the dimensionless PDR covering
    factor \f$f_{\text{PDR}}\f$.
*/
class SPHStellarComp : public StellarComp
{
    ITEM_CONCRETE(SPHStellarComp, StellarComp, "a stellar component derived from an SPH output file")

    PROPERTY_STRING(filename, "the name of the file with the SPH source particles")

    PROPERTY_BOOL(importVelocity, "import velocity information from the SPH particle file")
        ATTRIBUTE_DEFAULT_VALUE(importVelocity, "false")

    PROPERTY_ITEM(sedFamily, SEDFamily, "the SED family used to assign a spectral energy distribution to each particle")
        ATTRIBUTE_DEFAULT_VALUE(sedFamily, "BruzualCharlotSEDFamily")

    PROPERTY_BOOL(writeLuminosities, "output a data file with the luminosities per wavelength bin")
        ATTRIBUTE_DEFAULT_VALUE(writeLuminosities, "false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor. */
    ~SPHStellarComp();

protected:
    /** This function caches some pointers for later use. */
    void setupSelfBefore() override;

    /** This function performs setup for the SPH stellar component. The first step is to load the
        properties for each of the SPH source particles from the specified file. The second step is
        to calculate the luminosity \f$L_{\ell,i}\f$ for each particle \f$i\f$ and at each
        wavelength grid point \f$\ell\f$ through the %SED family specified as a property. Summing
        over all particles, a vector with the total luminosity for each wavelength bin is
        constructed. Finally, a matrix \f$X_{\ell,i}\f$ is filled that contains the normalized
        cumulative luminosity, \f[ X_{\ell,i} = \frac{ \sum_{j=0}^{i-1} L_{\ell,j} }{
        \sum_{j=0}^{N-1} L_{\ell,j} }. \f] This matrix will be used for the efficient generation
        of random photon packages from the stellar component. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the stellar component, which for this class is
        always 3 since there are no symmetries in the geometry. */
    int dimension() const override;

    /** This function returns the luminosity \f$L_\ell\f$ of the stellar component in the
        wavelength bin with index \f$\ell\f$. It just reads the appropriate number from the
        internally stored vector. */
    double luminosity(int ell) const override;

    /** This function simulates the emission of a monochromatic photon package with a luminosity
        \f$L\f$ at wavelength index \f$\ell\f$ from the stellar component. It randomly chooses an
        SPH particle from the \f$N\f$ possible particles by generating a random number
        \f${\cal{X}}\f$ and determining the particle number \f$i\f$ for which
        \f$X_{\ell,i}\leq{\cal{X}}<X_{\ell,i+1}\f$, with \f$X{\ell,i}\f$ the normalized cumulative
        luminosity matrix defined in the setup phase and stored internally. Once the SPH particle
        has been determined, a position is determined randomly from the smoothed distribution
        around the particle centre, a random propagation direction is determined, and a photon
        package with these properties is constructed and returned. The function assumes the scaled
        Gaussian smoothing kernel \f[ W(h,r) = \frac{a^3}{\pi^{3/2}\,h^3} \,\exp({-\frac{a^2
        r^2}{h^2}}) \f] with the empirically determined value of \f$a=2.42\f$, which approximates
        the standard cubic spline kernel to within two percent accuracy. */
    void launch(PhotonPackage* pp, int ell, double L) const override;

    //======================== Data Members ========================

private:
    // particle position and size
    std::vector<Vec> _rv;
    std::vector<double> _hv;

    // luminosity info
    Array _Ltotv;         // total luminosity for each wavelength bin -- [ell]
    ArrayTable<2> _Xvv;   // cumulative luminosity over particles, for each wavelength bin -- [ell, i]

    // anisotropy information for each particle (only if _velocity is true)
    std::vector<SPHStellarComp_Private::VelocityAnisotropy*> _av;  // [i]

    // cached
    Random* _random{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
