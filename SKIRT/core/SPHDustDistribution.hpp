/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHDUSTDISTRIBUTION_HPP
#define SPHDUSTDISTRIBUTION_HPP

#include "DustDistribution.hpp"
#include "DustMassInBoxInterface.hpp"
#include "DustParticleInterface.hpp"
#include "Array.hpp"
#include "DustMix.hpp"
#include "SPHGasParticle.hpp"
class SPHGasParticleGrid;

////////////////////////////////////////////////////////////////////

/** The SPHDustDistribution represents dust distributions defined from a set of SPH gas particles,
    such as for example resulting from a cosmological simulation. The information on the SPH gas
    particles is read from a text file formatted as follows. The text file should contain 6 or 7
    columns of numbers separated by whitespace; lines starting with # are ignored. The first three
    columns are the \f$x\f$, \f$y\f$ and \f$z\f$ coordinates of the particle (in pc), the fourth
    column is the SPH smoothing length \f$h\f$ (in pc), the fifth column is the mass \f$M\f$ of the
    particle (in \f$M_\odot\f$), and the sixth column is the metallicity \f$Z\f$ of the gas
    (dimensionless fraction). The optional seventh column is the temperature of the gas (in K). If
    this value is provided and it is higher than the maximum temperature the particle is ignored.
    If the temperature value is missing, the particle is never ignored. */
class SPHDustDistribution : public DustDistribution, public DustMassInBoxInterface, public DustParticleInterface
{
    ITEM_CONCRETE(SPHDustDistribution, DustDistribution, "a dust distribution derived from an SPH output file")

    PROPERTY_STRING(filename, "the name of the file with the SPH gas particles")

    PROPERTY_DOUBLE(dustFraction, "the fraction of the metal content locked up in dust grains")
        ATTRIBUTE_MIN_VALUE(dustFraction, "]0")
        ATTRIBUTE_MAX_VALUE(dustFraction, "1]")
        ATTRIBUTE_DEFAULT_VALUE(dustFraction, "0.3")

    PROPERTY_DOUBLE(maxTemperature, "the maximum temperature for a gas particle to contain dust")
        ATTRIBUTE_QUANTITY(maxTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(maxTemperature, "[0 K")
        ATTRIBUTE_MAX_VALUE(maxTemperature, "1000000 K]")
        ATTRIBUTE_DEFAULT_VALUE(maxTemperature, "75000 K")

    PROPERTY_ITEM(dustMix, DustMix, "the dust mix describing the attributes of the dust")
        ATTRIBUTE_DEFAULT_VALUE(dustMix, "InterstellarDustMix")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

    /** The destructor deletes the data structures allocated during setup. */
    ~SPHDustDistribution();

protected:
    /** This function performs setup for the SPH dust distribution. It reads the properties for each
        of the SPH gas particles from the specified file, converting them to program units and
        storing them in the internal vectors. */
    virtual void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the dust distribution, which for this class is always 3
        since there are no symmetries in the geometry. */
    int dimension() const override;

    /** This function returns the number of dust components that are involved in the dust
        distribution. For an SPH dust distribution this is equal to one. */
    int numComponents() const override;

    /** This function returns a pointer to the dust mixture corresponding to the \f$h\f$'th dust
        component. If \f$h\f$ is not equal to zero, an FatalError error is thrown. */
    DustMix* mix(int h) const override;

    /** This function returns the mass density \f$\rho_h({\bf{r}})\f$ of the \f$h\f$'th component
        of the dust distribution at the position \f${\bf{r}}\f$. If \f$h\f$ is not equal to zero,
        a FatalError error is thrown. In the other case, the call is passed to the total density
        function. */
    double density(int h, Position bfr) const override;

    /** This function returns the total mass density \f$\rho({\bf{r}})\f$ of the dust distribution
        at the position \f${\bf{r}}\f$. For an SPH dust distribution, the dust mass density is
        calculated by summing over all the particles \f[ \rho({\bf{r}}) = f_{\text{dust}} \sum_i
        Z_i\, M_i\, W(h_i,|{\bf{r}}-{\bf{r}}_i|) \f] with \f$f_{\text{dust}}\f$ the fraction of
        metals locked up in dust grains, \f$Z_i\f$ the metallicity and \f$M_i\f$ the
        (gas) mass of the \f$i\f$'th particle, \f$h_i\f$ the SPH smoothing length of the
        \f$i\f$'th particle, and \f$W(h,r)\f$ the SPH smoothing kernel. We assume a standard spline
        kernel, \f[ W(h,r) = \frac{8}{\pi\,h^3} \times \begin{cases} 1 - 6\,u^2\,(1-u) & \text{for
        }0<u<\tfrac12, \\ 2\,(1-u)^3 & \text{for }\tfrac12<u<1, \\ 0 & \text{else}. \end{cases} \f]
        with \f$u=r/h\f$. */
    double density(Position bfr) const override;

    /** This function generates a random position from the dust distribution. It randomly chooses a
        particle using the normalized cumulative density distribution constructed during the setup
        phase. Then a position is determined randomly from the smoothed distribution around the
        particle center. The function assumes the scaled Gaussian smoothing kernel \f[ W(h,r) =
        \frac{a^3}{\pi^{3/2}\,h^3} \,\exp({-\frac{a^2 r^2}{h^2}}) \f] with the empirically
        determined value of \f$a=2.42\f$, which approximates the standard cubic spline kernel to
        within two percent accuracy. */
    Position generatePosition() const override;

    /** This function returns the portion of the dust mass inside a given box (i.e. a cuboid lined
        up with the coordinate axes). If \f$h\f$ is not equal to zero, a FatalError error is
        thrown. In the other case, the call is passed to the total mass-in-box function. */
    double massInBox(int h, const Box& box) const override;

    /** This function returns the portion of the total dust mass (i.e. for all dust components)
        inside a given box (i.e. a cuboid lined up with the coordinate axes). For an SPH dust
        distribution, the dust mass inside the box is calculated by summing over all the particles
        and integrating over the box \f[ M_{\text{box}} = f_{\text{dust}} \sum_i Z_i\, M_i
        \int_{x_\text{min}}^{x_\text{max}} \int_{y_\text{min}}^{y_\text{max}}
        \int_{z_\text{min}}^{z_\text{max}} W(h_i,|{\bf{r}}-{\bf{r}}_i|) \,\text{d}x \,\text{d}y
        \,\text{d}z\f] with \f$f_{\text{dust}}\f$ the fraction of metals locked up in dust grains,
        \f$Z_i\f$ the metallicity and \f$M_i\f$ the (gas) mass of the \f$i\f$'th particle,
        \f$h_i\f$ the SPH smoothing length of the \f$i\f$'th particle, and \f$W(h,r)\f$ the SPH
        smoothing kernel. To speed up the calculations, this function uses the scaled Gaussian
        kernel \f[ W(h,r) = \frac{a^3}{\pi^{3/2}\,h^3} \,\exp({-\frac{a^2 r^2}{h^2}}) \f] with the
        empirically determined value of \f$a=2.42\f$ to make this kernel approximate the standard
        spline kernel to within two percent accuracy. The advantage of this kernel is that the
        integration over a box can be written in terms of the error function \f[
        \int_{x_\text{min}}^{x_\text{max}} \int_{y_\text{min}}^{y_\text{max}}
        \int_{z_\text{min}}^{z_\text{max}} W(h,\sqrt{x^2+y^2+z^2}) \,\text{d}x \,\text{d}y
        \,\text{d}z = \tfrac18 \left(\text{erf}(\frac{a\,x_\text{max}}{h}) -
        (\text{erf}(\frac{a\,x_\text{min}}{h})\right) \left(\text{erf}(\frac{a\,y_\text{max}}{h}) -
        (\text{erf}(\frac{a\,y_\text{min}}{h})\right) \left(\text{erf}(\frac{a\,z_\text{max}}{h}) -
        (\text{erf}(\frac{a\,y_\text{min}}{h})\right) \f] */
    double massInBox(const Box& box) const override;

    /** This function returns the total mass of the \f$h\f$'th component of the dust distribution. If
        \f$h\f$ is not equal to zero, a FatalError error is thrown. In the other case, the call is passed
        to the total mass. */
    double mass(int h) const override;

    /** This function returns the total dust mass of the dust distribution. For a SPH dust
        distribution, the total dust mass is calculated as \f[ M = f_{\text{dust}} \sum_i Z_i\, M_i
        \f] with \f$f_{\text{dust}}\f$ the fraction of metals locked up in dust grains, and
        \f$Z_i\f$ the metallicity and \f$M_i\f$ the (gas) mass of the \f$i\f$'th particle. */
    double mass() const override;

    /** This function returns the X-axis surface density of the dust distribution. For an SPH
        dust distribution, this integral is calculated numerically using 10000 samples along
        the X-axis. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density of the dust distribution. For an SPH
        dust distribution, this integral is calculated numerically using 10000 samples along
        the Y-axis. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density of the dust distribution. For an SPH
        dust distribution, this integral is calculated numerically using 10000 samples along
        the Z-axis. */
    double SigmaZ() const override;

    /** This function implements the DustParticleInterface. It returns the number of SPH particles
        defining this dust distribution. */
    int numParticles() const override;

    /** This function implements the DustParticleInterface. It returns the coordinates of the SPH
        particle with the specified zero-based index. If the index is out of range, a fatal error
        is thrown. */
    Vec particleCenter(int index) const override;

    /** This function is used by the interface() template function in the SimulationItem class. It
        returns a list of simulation items that should be considered in the search for an item that
        implements the requested interface. The implementation in this class returns the default
        list (i.e. the receiving object) except in the following case. If the requested interface
        is DustMassInBoxInterface (which is implemented by this class) and one or more of the
        imported SPH particles have a negative mass, the returned list is empty. We avoid using the
        (faster) mass-in-box calculation when there are negative masses because the result can be
        inconsistent with the densities sampled over random points throughout the volume. This
        inconsistency is caused by the fact that the sampled densities are individually clipped to
        zero, while the mass-in-box is only clipped to zero after the integration/summation has
        been performed. */
    vector<SimulationItem*> interfaceCandidates(const std::type_info& interfaceTypeInfo) override;

    //======================== Data Members ========================

private:
    // the SPH particles
    vector<SPHGasParticle> _pv;     // the particles in the order read from the file
    const SPHGasParticleGrid* _grid{nullptr};  // a list of particles overlapping each grid cell
    Array _cumrhov;                 // cumulative density distribution for particles in pv
    bool _negativeMasses{false};    // true if at least one of the imported particles has a negative mass
};

////////////////////////////////////////////////////////////////////

#endif
