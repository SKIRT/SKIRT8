/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FRAMEINSTRUMENT_HPP
#define FRAMEINSTRUMENT_HPP

#include "SingleFrameInstrument.hpp"
#include "ParallelDataCube.hpp"

////////////////////////////////////////////////////////////////////

/** A FrameInstrument object represents a basic instrument that just records the total flux in
    every pixel, and outputs this as a data cube in a FITS file. Internally, the class contains a
    single 3D vector (the f-vector) corresponding to the surface brightness in every pixel, at
    every wavelength index. */
class FrameInstrument : public SingleFrameInstrument
{
    ITEM_CONCRETE(FrameInstrument, SingleFrameInstrument,
                  "a basic instrument that outputs the total flux in every pixel as a data cube")
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
    ParallelDataCube _distftotv;
};

////////////////////////////////////////////////////////////////////

#endif
