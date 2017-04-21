/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ALLCELLSDUSTLIB_HPP
#define ALLCELLSDUSTLIB_HPP

#include "DustLib.hpp"

//////////////////////////////////////////////////////////////////////

/** The AllCellsDustLib class determines the relevant dust emission spectra for the simulation by
    simply performing the calculation for each and every dust cell. This straightforward mechanism
    is meaningful when the calculation is fast (for example when assuming LTE conditions), or as a
    reference for evaluating more complicated library mechanisms. */
class AllCellsDustLib : public DustLib
{
    ITEM_CONCRETE(AllCellsDustLib, DustLib, "a dust library that calculates the emissivity for every dust cell")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function returns the number of entries in the library. In this class the function
        returns the number of dust cells. */
    int numEntries() const override;

    /** This function returns a vector \em nv with length \f$N_{\text{cells}}\f$ that maps each
        cell \f$m\f$ to the corresponding library entry \f$n_m\f$. In this class the function
        returns the identity mapping. */
    vector<int> mapping() const override;
};

////////////////////////////////////////////////////////////////////

#endif
