/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SMOOTHINGKERNEL_HPP
#define SMOOTHINGKERNEL_HPP

#include "SimulationItem.hpp"

//////////////////////////////////////////////////////////////////////

/** SmoothingKernel is an abstract class that represents smoothing kernels that can be used to
    smear out point sources to smoothed particles. Each smoothing kernel is completely defined by
    the kernel density \f$W(u)\f$ which is just a function of the normalized radius \f$u=r/h\f$
    with \f$h\f$ the smoothing length. All smoothing kernels must be spherically symmetric and
    normalized to one, i.e. \f[ 4\pi \int W(u)\, u^2\, {\text{d}}u = 1. \f] */
class SmoothingKernel : public SimulationItem
{
    ITEM_ABSTRACT(SmoothingKernel, SimulationItem, "a smoothing kernel")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$W(u)\f$ of the smoothing kernel as a function of the
        normalized radius \f$u\f$. Subclasses must implement this function appropriately. */
    virtual double density(double u) const = 0;

    /** This pure virtual function generates a random normalized radius \f$u\f$ from the smoothing
        kernel, by drawing a number from the one-dimensional probability density \f$
        p(u)\,{\text{d}}u = 4\pi\,W(u)\,u^2\, {\text{d}}u \f$. Subclasses must implement this
        function appropriately. */
    virtual double generateRadius() const = 0;

};

//////////////////////////////////////////////////////////////////////

#endif
