/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PANDUSTSYSTEM_HPP
#define PANDUSTSYSTEM_HPP

#include "DustSystem.hpp"
#include "DustEmissivity.hpp"
#include "DustLib.hpp"
#include "ParallelTable.hpp"
class ProcessAssigner;

//////////////////////////////////////////////////////////////////////

/** A PanDustSystem class object represents a complete dust system for use with panchromatic
    simulations. This class relies on the functionality implemented in the DustSystem base class,
    and additionaly supports dust emission. It maintains information on the absorbed energy for
    each cell at each wavelength in a (potentially very large) table. It also holds a
    DustEmissivity object and a DustLib object used to calculate the dust emission spectrum for
    dust cells. */
class PanDustSystem : public DustSystem
{
    ITEM_CONCRETE(PanDustSystem, DustSystem, "a dust system for use with panchromatic simulations")

    PROPERTY_ITEM(dustEmissivity, DustEmissivity, "the dust emissivity type")
        ATTRIBUTE_DEFAULT_VALUE(dustEmissivity, "GreyBodyDustEmissivity")
        ATTRIBUTE_OPTIONAL(dustEmissivity)

    PROPERTY_ITEM(dustLib, DustLib, "the dust library mechanism")
        ATTRIBUTE_RELEVANT_IF(dustLib, "dustEmissivity")
        ATTRIBUTE_DEFAULT_VALUE(dustLib, "AllCellsDustLib")

    PROPERTY_BOOL(includeSelfAbsorption, "include dust self-absorption")
        ATTRIBUTE_RELEVANT_IF(includeSelfAbsorption, "dustEmissivity")
        ATTRIBUTE_DEFAULT_VALUE(includeSelfAbsorption, "false")

    PROPERTY_BOOL(writeTemperature, "output FITS files displaying the dust temperature distribution")
        ATTRIBUTE_RELEVANT_IF(writeTemperature, "dustEmissivity")
        ATTRIBUTE_DEFAULT_VALUE(writeTemperature, "true")

    PROPERTY_DOUBLE(emissionBias, "the dust emission bias")
        ATTRIBUTE_RELEVANT_IF(emissionBias, "dustEmissivity")
        ATTRIBUTE_MIN_VALUE(emissionBias, "[0")
        ATTRIBUTE_MAX_VALUE(emissionBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(emissionBias, "0.5")
        ATTRIBUTE_SILENT(emissionBias)

    PROPERTY_DOUBLE(emissionBoost, "the factor by which to boost the number of dust emission photon packages")
        ATTRIBUTE_RELEVANT_IF(emissionBoost, "dustEmissivity")
        ATTRIBUTE_MIN_VALUE(emissionBoost, "]0")
        ATTRIBUTE_MAX_VALUE(emissionBoost, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(emissionBoost, "1")
        ATTRIBUTE_SILENT(emissionBoost)

    PROPERTY_INT(minIterations, "the minimum number of dust self-absorption iterations")
        ATTRIBUTE_RELEVANT_IF(minIterations, "includeSelfAbsorption")
        ATTRIBUTE_MIN_VALUE(minIterations, "1")
        ATTRIBUTE_MAX_VALUE(minIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minIterations, "1")
        ATTRIBUTE_SILENT(minIterations)

    PROPERTY_INT(maxIterations, "the maximum number of dust self-absorption iterations")
        ATTRIBUTE_RELEVANT_IF(maxIterations, "includeSelfAbsorption")
        ATTRIBUTE_MIN_VALUE(maxIterations, "1")
        ATTRIBUTE_MAX_VALUE(maxIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxIterations, "10")
        ATTRIBUTE_SILENT(maxIterations)

    PROPERTY_DOUBLE(maxFractionOfStellar, "convergence is reached when the total absorbed dust luminosity "
                                          "is less than this fraction of the total absorbed stellar luminosity")
        ATTRIBUTE_RELEVANT_IF(maxFractionOfStellar, "includeSelfAbsorption")
        ATTRIBUTE_MIN_VALUE(maxFractionOfStellar, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfStellar, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfStellar, "0.01")
        ATTRIBUTE_SILENT(maxFractionOfStellar)

    PROPERTY_DOUBLE(maxFractionOfPrevious, "convergence is reached when the total absorbed dust luminosity "
                                           "has changed by less than this fraction compared to the previous iteration")
        ATTRIBUTE_RELEVANT_IF(maxFractionOfPrevious, "includeSelfAbsorption")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPrevious, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPrevious, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPrevious, "0.03")
        ATTRIBUTE_SILENT(maxFractionOfPrevious)

    PROPERTY_BOOL(writeEmissivity, "output a data file with the dust mix emissivities in the local ISRF")
        ATTRIBUTE_RELEVANT_IF(writeEmissivity, "dustEmissivity")
        ATTRIBUTE_DEFAULT_VALUE(writeEmissivity, "false")
        ATTRIBUTE_SILENT(writeEmissivity)

    PROPERTY_BOOL(writeISRF, "output a data file describing the interstellar radiation field")
        ATTRIBUTE_RELEVANT_IF(writeISRF, "dustEmissivity")
        ATTRIBUTE_DEFAULT_VALUE(writeISRF, "false")
        ATTRIBUTE_SILENT(writeISRF)

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function does some basic initialization and verifies that the property values are
        valid. */
    void setupSelfBefore() override;

    /** This function properly sizes the tables for storing absorption information, corresponding
        to the property settings for this simulation.

        If the writeEmissivity flag is enabled, this function outputs text files tabulating the
        emissivity for each dust component's dust mix, assuming the dust would be embedded in a a
        range of scaled local (i.e. solar neighborhood) interstellar radiation fields as defined by
        Mathis et al. (1983, A&A, 128, 212), and a range of diluted black body input fields. */
    void setupSelfAfter() override;

    //======== Setters & Getters for Discoverable Attributes =======

    /** \fn emissionBias
        The emission bias is the fraction of dust cells selected for emission from a uniform
        distribution rather than a distribution weighted according to the total dust luminosity of
        the cells. */

    /** \fn emissionBoost
        The emission boost is the multiplication factor by which to increase the number of photon
        packages sent during the dust emission phase. The default value is 1, i.e. the same number
        of photon packages are used as during the stellar emission phase. A higher value increases
        the resolution of infrared images at the cost of extra run-time. */

public:
    /** This function puts the dust system in emulation mode. Specifically, it sets an internal
        flag that can be queried by other classes, and if dust self-absorption is enabled, it
        forces the number of self-absorption iterations to one. */
    void setEmulationMode();

    /** This function returns true if the dust system has been put in emulation mode. */
    bool emulationMode();

    //======================== Other Functions =======================

public:
    /** This function returns the data parallelization process assigner for this dust system, which
        is created during setup. The process assigner is the object that assigns different dust
        cells to different processes, to parallelize the storage of the dust absorption and
        emission. If the simulation is not in data parallelization mode, the function returns the
        null pointer. */
    const ProcessAssigner* assigner() const;

    /** This function returns true if dust emission is turned on for this dust system, and false
        otherwise. */
    bool hasDustEmission() const override;

    /** This function returns a flag that indicate whether the absorption rates in each cell need
        to be stored for this dust system. For a panchromatic simulation, absorption rates are only
        calculated if dust emission is turned on. */
    bool hasDustAbsorption() const override;

    /** The function simulates the absorption of a monochromatic luminosity package in the dust
        cell with cell number \f$m\f$, i.e. it adds a fraction \f$\Delta L\f$ to the already
        absorbed luminosity at wavelength index \f$\ell\f$. The function adds the absorbed
        luminosity \f$\Delta L\f$ to the appropriate item in the absorption rate table for stellar
        or dust emission as indicated by the flag. The addition is performed in a thread-safe
        manner so this function may be concurrently called from multiple threads. */
    void absorb(int m, int ell, double DeltaL, bool ynstellar) override;

    /** This function resets the absorbed dust luminosity to zero in all cells of the dust system.
        */
    void resetDustAbsorption();

    /** This function returns the absorbed luminosity \f$L_{\ell,m}\f$ at wavelength index
        \f$\ell\f$ in the dust cell with cell number \f$m\f$. For a panchromatic dust system, it sums
        the individual absorption rate counters corresponding to the stellar and dust emission. */
    double absorbedLuminosity(int m, int ell) const override;

    /** This function returns the total (bolometric) absorbed luminosity in the dust cell with cell
        number \f$m\f$. It is calculated by summing the absorbed luminosity at all the wavelength
        indices. */
    double absorbedLuminosity(int m) const;

    /** This function returns a vector with the total (bolometric) absorbed luminosity in each dust cell. */
    Array absorbedLuminosity() const;

    /** This function returns the total (bolometric) absorbed dust luminosity in the entire dust system.
        It is calculated by summing the absorbed stellar luminosity of all the cells. */
    double absorbedStellarLuminosity() const;

    /** This function returns the total (bolometric) absorbed luminosity in the entire dust system.
        It is calculated by summing the absorbed dust luminosity of all the cells. */
    double absorbedDustLuminosity() const;

    /** This function (re-)calculates the relevant dust emission spectra for the dust system, based
        on the absorption data currently stored in the dust cells, and internally caches the
        results. If dust emission is turned off, this function does nothing. */
    void calculateDustEmission();

    /** This function synchronizes the results of the absorption by calling the sync() function on the
        absorption tables. **/
    void sumResults();

    /** This function returns the luminosity \f$L_\ell\f$ at the wavelength index \f$\ell\f$ in the
        normalized dust emission SED corresponding to the dust cell with dust cell number \f$m\f$.
        It just looks up the appropriate value in the cached results produced by calculate(). If
        dust emission is turned off, the function returns zero. */
    double emittedDustLuminosity(int m, int ell) const;

    /** If the writeTemperature flag is enabled, this function writes out FITS files (named
        <tt>prefix_ds_tempXX.fits</tt>) with the mean dust temperatures in the coordinate planes.
        Each of these maps contains 1024 x 1024 pixels, and covers as a field of view the total
        extension of the grid. The number of data files written depends on the dimension of the
        dust system's geometry: for spherical symmetry only the intersection with the xy plane is
        written, for axial symmetry the intersections with the xy and xz planes are written, and
        for general geometries all three intersections are written. Each FITS file is a data cube,
        where a map is created for each dust population of each dust component. The mean dust
        temperature \f${\bar{T}}_c({\boldsymbol{r}})\f$ corresponding to the \f$c\f$'th population
        at the position \f${\boldsymbol{r}}\f$ is calculated from the mean intensity
        \f$J_\lambda({\boldsymbol{r}})\f$ of the radiation field using the balance equation \f[
        \int_0^\infty \kappa_{\lambda,c}^{\text{abs}}\, J_\lambda({\boldsymbol{r}})\,
        {\text{d}}\lambda = \int_0^\infty \kappa_{\lambda,c}^{\text{abs}}\,
        B_\lambda({\bar{T}}_c({\boldsymbol{r}}))\, {\text{d}}\lambda, \f] where
        \f$\kappa_{\lambda,c}^{\text{abs}}\f$ is the absorption coefficient of the \f$c\f$'th dust
        population. Note that this mean dust temperature is only illustrative: it does not imply
        that the dust emits as a modified blackbody at an equibrium temperature.

        Furthermore, if the writeTemperature flag is enabled, this function writes out a text
        file (named <tt>prefix_ds_celltemps.dat</tt>) containing an indicative dust temperature for
        each cell in the dust grid. Each line contains the cell index \f$m\f$ and the indicative
        dust temperature, defined as the mean temperature (calculated as described above) averaged
        over all dust populations (weighed by mass in the dust mix) and all dust components
        (weighed by density in the dust cell).

        If the writeISRF flag is enabled, this function writes out a text file (named
        <tt>prefix_ds_isrf.dat</tt>) describing the interstellar radiation field for each cell in
        the dust grid. Each line has the following columns: the cell index \f$m\f$, the bolometric
        absorbed luminosity \f$L_\text{abs}\f$, and the mean radiation field \f$J_{\ell,m}\f$ for
        each wavelength \f$\lambda_\ell\f$ in the simulation's wavelength grid. */
    void write() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    int _Nlambda{0};
    ParallelTable _LabsStelvv;
    ParallelTable _LabsDustvv;
    bool _haveLabsStel{false};     // true if absorbed stellar emission is relevant for this simulation
    bool _haveLabsDust{false};     // true if absorbed dust emission is relevant for this simulation
    const ProcessAssigner* _assigner{nullptr}; // determines which cells will be given to the DustLib

    // data member to remember whether emulation mode is enabled
    bool _emulationMode{false};
};

//////////////////////////////////////////////////////////////////////

#endif
