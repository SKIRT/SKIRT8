/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DIM2DUSTLIB_HPP
#define DIM2DUSTLIB_HPP

#include "DustLib.hpp"

//////////////////////////////////////////////////////////////////////

/** The Dim2DustLib class calculates a relevant subset of the dust emission spectra for the
    simulation, and maps each dust cell to one of these library entries. This avoids performing the
    calculation explicitly for every dust cell in the simulation, in return for lack of accuracy.
    The library is built using two parameters that describe the interstellar radiation fields
    \f$J_\lambda\f$ in the \f$N_{\text{cells}}\f$ cells of the dust system: a mean temperature and
    a mean wavelength. For a dust cell \f$m\f$ with dust components \f$h\f$, the mean temperature
    is defined as \f[{\bar{T}}_m = \frac{\sum_h \rho_{m,h}\,{\bar{T}}_{m,h}} {\sum_h \rho_{m,h}}
    \f] where \f${\bar{T}}_{m,h}\f$ is the mean dust temperature for dust component \f$h\f$,
    obtained through the balance equation \f[ \int_0^\infty \varsigma_{h,\lambda}^{\text{abs}}\,
    J_{m,\lambda}\, {\text{d}}\lambda = \int_0^\infty \varsigma_{h,\lambda}^{\text{abs}}\,
    B_\lambda({\bar{T}}_{m,h})\, {\text{d}}\lambda. \f] The mean wavelength is similarly defined as
    \f[{\bar{\lambda}}_m = \frac{\sum_h \rho_{m,h}\,{\bar{\lambda}}_{m,h}} {\sum_h \rho_{m,h}} \f]
    where \f${\bar{\lambda}}_{m,h}\f$ is the mean wavelength for dust component \f$h\f$, given by
    \f[ {\bar{\lambda}}_{m,h} = \frac{\int_0^\infty \varsigma_{h,\lambda}^{\text{abs}}\,
    J_{m,\lambda}\, \lambda\, {\text{d}}\lambda }{ \int_0^\infty
    \varsigma_{h,\lambda}^{\text{abs}}\, J_{m,\lambda}\, {\text{d}}\lambda}. \f] The library is
    built from binning the \f$N_{\text{cells}}\f$ values of \f${\bar{T}}_m\f$ and
    \f${\bar{\lambda}}_m\f$ onto a two-dimensional grid. Each grid point now corresponds to an
    entry in the library, and the corresponding dust emissivity is calculated for each library
    entry. */
class Dim2DustLib : public DustLib
{
    ITEM_CONCRETE(Dim2DustLib, DustLib, "a dust library with a two-dimensional grid of emissivity entries")

    PROPERTY_INT(numTemperatures, "the number of mean temperature grid points")
        ATTRIBUTE_MIN_VALUE(numTemperatures, "3")
        ATTRIBUTE_MAX_VALUE(numTemperatures, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numTemperatures, "25")

    PROPERTY_INT(numWavelengths, "the number of mean wavelength grid points")
        ATTRIBUTE_MIN_VALUE(numWavelengths, "3")
        ATTRIBUTE_MAX_VALUE(numWavelengths, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numWavelengths, "10")

    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function returns the number of entries in the library. In this class the function
        returns the product of the number of mean temperature grid points and the number of mean
        wavelength grid points. */
    int numEntries() const override;

    /** This function returns a vector \em nv with length \f$N_{\text{cells}}\f$ that maps each
        cell \f$m\f$ to the corresponding library entry \f$n_m\f$. In this class the function loops
        over all dust cells and calculates \f${\bar{T}}_m\f$ and \f${\bar{\lambda}}_m\f$ for each
        dust cell \f$m\f$. A two-dimensional grid is established such that it fits all the measured
        values. The mean temperature grid points \f${\bar{T}}_{(i)}\f$ are distributed linearly,
        i.e. \f[ {\bar{T}}_{(i)} = {\bar{T}}_{\text{min}} + \frac{i}{N_{\bar{T}}}\,
        ({\bar{T}}_{\text{max}} - {\bar{T}}_{\text{min}}) \qquad i=0,\ldots,N_{\bar{T}} \f] where
        \f${\bar{T}}_{\text{min}}\f$ and \f${\bar{T}}_{\text{max}}\f$ represent the smallest and
        largest values of the mean temperature found among all dust cells. The mean wavelength grid
        points have logarithmic distribution, \f[ {\bar{\lambda}}_{(j)} =
        {\bar{\lambda}}_{\text{min}} \left( \frac{ {\bar{\lambda}}_{\text{max}} }{
        {\bar{\lambda}}_{\text{min}} } \right)^{j/N_{\bar{\lambda}}} \qquad
        j=0,\ldots,N_{\bar{\lambda}}. \f] The function then calculates for each cell \f$m\f$ its
        library entry \f$n \equiv (i,j)\f$. */
    vector<int> mapping() const override;
};

////////////////////////////////////////////////////////////////////

#endif
