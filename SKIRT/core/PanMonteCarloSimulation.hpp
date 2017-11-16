/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PANMONTECARLOSIMULATION_HPP
#define PANMONTECARLOSIMULATION_HPP

#include "MonteCarloSimulation.hpp"
#include "Array.hpp"
#include "PanDustSystem.hpp"
#include "PanWavelengthGrid.hpp"
#include "StellarSystem.hpp"

//////////////////////////////////////////////////////////////////////

/** This is a subclass of the general MonteCarloSimulation class representing a panchromatic Monte
    Carlo simulation, i.e. operating at a range of wavelengths. In such simulations, there can be
    absorption, scattering and thermal emission by dust grains. */
class PanMonteCarloSimulation : public MonteCarloSimulation
{
    ITEM_CONCRETE(PanMonteCarloSimulation, MonteCarloSimulation, "a panchromatic Monte Carlo simulation")

    PROPERTY_ITEM(wavelengthGrid, PanWavelengthGrid, "the wavelength grid")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthGrid, "LogWavelengthGrid")

    PROPERTY_ITEM(stellarSystem, StellarSystem, "the stellar system")
        ATTRIBUTE_DEFAULT_VALUE(stellarSystem, "StellarSystem")

    PROPERTY_ITEM(dustSystem, PanDustSystem, "the dust system")
        ATTRIBUTE_DEFAULT_VALUE(dustSystem, "PanDustSystem")
        ATTRIBUTE_OPTIONAL(dustSystem)

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function performs some basic initialization. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

protected:
    /** This function actually runs the simulation. For a panchromatic simulation, this includes
        the stellar emission phase, the dust self-absorption phase, and the dust emission phase
        (plus writing the results). */
    void runSelf() override;

private:
    /** This function drives the dust self-absorption phase in a panchromatic Monte Carlo
        simulation. The function implements a big loop that iterates over consecutive
        self-absorption calculations until a self-consistent state is reached.

        The minimum and maximum number of iterations can be specified as silent options of the
        PanDustSystem object in the simulation. Within these limits, the actual number of
        iterations performed is determined by convergence criteria which can also be specified as
        silent options of the PanDustSystem object. Convergence is reached (and the function exits)
        when (a) the absorbed dust luminosity is less than a given fraction of the absorbed stellar
        luminosity, \em OR (b) the absorbed dust luminosity has changed by less than a given
        fraction compared to the previous iteration.

        Within every iteration, the first task is to construct the dust emission library that
        describes the spectral properties of the dust emission. Subsequently the dust
        self-absorption phase implements a parallelized loop that iterates over
        \f$N_{\text{pp}}\times N_\lambda\f$ monochromatic photon packages. Within this loop, the
        function simulates the life cycle of a single dust photon package. Before we can start the
        life cycle of the photon packages at a given wavelength index \f$\ell\f$, we first have to
        determine the total luminosity that is emitted from every dust cell at that wavelength
        index. This is just the product of the total luminosity absorbed in the \f$m\f$'th dust
        cell, \f$L^{\text{abs}}_m\f$, and the normalized SED at wavelength index \f$\ell\f$
        corresponding to that cell, as obtained from the dust emission library. Once we know the
        luminosity \f$L_{\ell,m}\f$ emitted by each dust cell, we calculate the total dust
        luminosity, \f[ L_\ell = \sum_{m=0}^{N_{\text{cells}}-1} L_{\ell,m}, \f] and create a
        vector \f$X_m\f$ that describes the normalized cumulative luminosity distribution as a
        function of the cell number \f$m\f$, \f[ X_m = \frac{ \sum_{m'=0}^m L_{\ell,m'} }{ L_\ell
        }. \f] This vector is used to generate random dust cells from which photon packages can be
        launched. Now the actual dust self-absorption can start, i.e. we launch \f$N_{\text{pp}}\f$
        different photon packages at wavelength index \f$\ell\f$, with the original position chosen
        as a random position in the cell \f$m\f$ chosen randomly from the cumulative luminosity
        distribution \f$X_m\f$. The remaining life cycle of a photon package in the dust emission
        phase is very similar to the life cycle described in
        MonteCarloSimulation::runstellaremission(). */
    void runDustSelfAbsorption();

    /** This function implements the loop body for rundustselfabsorption(). */
    void doDustSelfAbsorptionChunk(size_t index);

    /** This function drives the dust emission phase in a panchromatic Monte Carlo simulation. The
        first task is to construct the dust emission library that describes the spectral properties
        of the dust emission. Subsequently the dust emission phase implements a parallelized loop
        that iterates over \f$N_{\text{pp}}\times N_\lambda\f$ monochromatic photon packages.
        Within this loop, the function simulates the life cycle of a single dust photon
        package. Before we can start the life cycle of the photon packages at a given wavelength
        index \f$\ell\f$, we first have to determine the total luminosity that is emitted from
        every dust cell at that wavelength index. This is just the product of the total luminosity
        absorbed in the \f$m\f$'th dust cell, \f$L^{\text{abs}}_m\f$, and the normalized SED at
        wavelength index \f$\ell\f$ corresponding to that cell, as obtained from the dust emission
        library. Once we know the luminosity \f$L_{\ell,m}\f$ emitted by each dust cell, we
        calculate the total dust luminosity, \f[ L_\ell = \sum_{m=0}^{N_{\text{cells}}-1}
        L_{\ell,m}, \f] and create a vector \f$X_m\f$ that describes the normalized cumulative
        luminosity distribution as a function of the cell number \f$m\f$, \f[ X_m = \frac{
        \sum_{m'=0}^m L_{\ell,m'} }{ L_\ell }. \f] This vector is used to generate random dust
        cells from which photon packages can be launched. Now the actual dust emission can start,
        i.e. we launch \f$N_{\text{pp}}\f$ different photon packages at wavelength index
        \f$\ell\f$, with the original position chosen as a random position in the cell \f$m\f$
        chosen randomly from the cumulative luminosity distribution \f$X_m\f$. The remaining life
        cycle of a photon package in the dust emission phase is very similar to the life cycle
        described in MonteCarloSimulation::runstellaremission(). */
    void runDustEmission();

    /** This function implements the loop body for rundustemission(). */
    void doDustEmissionChunk(size_t index);

    //======================== Data Members ========================

private:
    PanDustSystem*& _pds{_dustSystem};   // copy of _ds (pointer to the dust system) with the "pan" subtype

    // data members used to communicate between rundustXXX() and the corresponding parallel loop
    int _Ncells{0};        // number of dust cells
    Array _Labsbolv;       // vector that contains the bolometric absorbed luminosity in each cell
};

////////////////////////////////////////////////////////////////////

#endif
