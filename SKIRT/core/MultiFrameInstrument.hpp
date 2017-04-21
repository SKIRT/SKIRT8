/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MULTIFRAMEINSTRUMENT_HPP
#define MULTIFRAMEINSTRUMENT_HPP

#include "DistantInstrument.hpp"
#include "InstrumentFrame.hpp"

////////////////////////////////////////////////////////////////////

/** MultiFrameInstrument is a specialty instrument class for use with oligochromatic simulations in
    combination with an external tool such as FitSKIRT. It is similar to the FrameInstrument class
    in the sense that each pixel stores the incoming total flux per wavelength. However,
    MultiFrameInstrument allows a different frame extent and/or pixel resolution at each
    wavelength. All frames do share the direction and position angles determined by the properties
    of the DistantInstrument base class. It is assumed that the distance to the system is
    sufficiently large so that parallel projection can be used. */

class MultiFrameInstrument : public DistantInstrument
{
    ITEM_CONCRETE(MultiFrameInstrument, DistantInstrument,
                  "an instrument with a different frame per wavelength (for use with FitSKIRT)")
        ATTRIBUTE_ALLOWED_IF(MultiFrameInstrument, "OligoMonteCarloSimulation")

    PROPERTY_BOOL(writeTotal, "output the total flux")
        ATTRIBUTE_DEFAULT_VALUE(writeTotal, "false")

    PROPERTY_BOOL(writeStellarComps, "output the flux emitted from each stellar component seperately")
        ATTRIBUTE_DEFAULT_VALUE(writeStellarComps, "true")

    PROPERTY_ITEM_LIST(frames, InstrumentFrame, "the instrument frames (one for each wavelength)")
        ATTRIBUTE_DEFAULT_VALUE(frames, "InstrumentFrame")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that all attribute values have been appropriately set and performs
        setup for the instrument. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

    /** This function simulates the detection of a photon package by the instrument. It operates
        similarly to SimpleInstrument::detect(), except that the photon packages for different
        wavelengths are handed to different instrument frames. */
    void detect(PhotonPackage* pp) override;

    /** This function calibrates and outputs the instrument data. It operates similarly to
        SimpleInstrument::write(), except that a separate output file is written for each
        wavelength, using filenames that include the wavelength index \f$\ell\f$. */
    void write() override;
};

////////////////////////////////////////////////////////////////////

#endif
