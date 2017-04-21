/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SEDINSTRUMENT_HPP
#define SEDINSTRUMENT_HPP

#include "DistantInstrument.hpp"

////////////////////////////////////////////////////////////////////

/** An SEDInstrument object represents a basic instrument that just records the total integrated
    flux, and outputs this as an SED file. Internally, the class contains a single 1D vector (the
    F-vector) that stores the total integrated flux at every wavelength index. */
class SEDInstrument : public DistantInstrument
{
    ITEM_CONCRETE(SEDInstrument, DistantInstrument,
                  "a basic instrument that outputs the total integrated flux as an SED")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function completes setup for this instrument. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

protected:
    /** This function simulates the detection of a photon package by the instrument.
        See SimpleInstrument::detect() for more information. */
    void detect(PhotonPackage* pp) override;

    /** This function calibrates and outputs the instrument data.
        See SimpleInstrument::write() for more information. */
    void write() override;

    //======================== Data Members ========================

private:
    Array _Ftotv;
};

////////////////////////////////////////////////////////////////////

#endif
