/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WAVELENGTHGRID_HPP
#define WAVELENGTHGRID_HPP

#include "SimulationItem.hpp"
#include "Array.hpp"
class ProcessAssigner;

//////////////////////////////////////////////////////////////////////

/** WavelengthGrid is an abstract class that represents grids of wavelengths on which Monte Carlo
    simulations are defined. Each WavelengthGrid class object basically consists of a set of
    wavelength grid points \f$\lambda_\ell\f$ and the corresponding wavelength bin widths
    \f$\Delta_\ell\f$. A particular Monte Carlo simulation uses a single instance of this class. */
class WavelengthGrid : public SimulationItem
{
    ITEM_ABSTRACT(WavelengthGrid, SimulationItem, "a wavelength grid")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the wavelength bins have been initialized by a subclass calling
        the setWavelengthBins() function of this class in their setupSelfBefore() function. */
    void setupSelfAfter() override;

    /** This function sets the wavelength grid to the central wavelengths and corresponding bin
        widths given as arguments. It should be called from the setupSelfBefore() function in each
        WavelengthGrid subclass. The function verifies that both arrays have the same nonzero
        length and that the wavelengths are sorted in ascending order. If these requirements are
        not met, the function throws a fatal error. */
    void setWavelengthBins(const Array& lambdav, const Array& dlambdav);

    //======================== Other Functions =======================

public:
    /** This function returns true if the wavelength grid represents a sampled wavelength range (as
        required for panchromatic simulations), and false if it contains individual distinct
        wavelengths (as used by oligochromatic simulations). Subclasses must implement this
        function appropriately. */
    virtual bool isSampledRange() const = 0;

    /** This function returns the wavelength assigner for the simulation, or the null pointer if
        the simulation is not in data parallel multi-processing mode. The process assigner is the
        object that assigns different wavelengths to different processes, to parallelize the photon
        shooting algorithm, and the data storage. */
    const ProcessAssigner* assigner() const;

    /** This function returns the number of wavelength grid points in the grid. */
    int numWavelengths() const;

    /** This function returns the wavelength \f$\lambda_\ell\f$ corresponding to the index
        \f$\ell\f$. */
    double lambda(int ell) const;

    /** This function returns the width of the \f$\ell\f$'th wavelength bin. */
    double dlambda(int ell) const;

    /** This function returns the minimum border of the wavelength bin corresponding to the index
        \f$\ell\f$. */
    double lambdamin(int ell) const;

    /** This function returns the maximum border of the wavelength bin corresponding to the index
        \f$\ell\f$. */
    double lambdamax(int ell) const;

    /** This function returns the index \f$\ell\f$ of the grid point to which the wavelength
        \f$\lambda\f$ naturally belongs. It is the natural number \f$\ell\f$ that corresponds to
        the grid point \f$\lambda_\ell\f$ that lies nearest to \f$\lambda\f$ in logarithmic space.
        If \f$\lambda\f$ lies outside the wavelength interval, the value \f$\ell=-1\f$ is returned.
        */
    int nearest(double lambda) const;

    /** This function returns the entire table with the wavelength grid points. */
    const Array& lambdav() const;

    /** This function returns the entire table with the wavelength bin widths. */
    const Array& dlambdav() const;

    //======================== Data Members ========================

private:
    // subclasses should fill these tables by calling setWavelengthBins() in their setupSelfBefore()
    Array _lambdav;     // central wavelength for each bin
    Array _dlambdav;    // width of each wavelength bin

    // determines which wavelengths are assigned to each process
    ProcessAssigner* _assigner{nullptr};
};

//////////////////////////////////////////////////////////////////////

#endif
