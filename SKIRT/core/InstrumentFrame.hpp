/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INSTRUMENTFRAME_HPP
#define INSTRUMENTFRAME_HPP

#include "SimulationItem.hpp"
#include "ArrayTable.hpp"
class PhotonPackage;
class MultiFrameInstrument;

////////////////////////////////////////////////////////////////////

/** The InstrumentFrame class implements a single instrument frame with a number of pixels,
    field-of-view and center. It is used by the MultiFrameInstrument class to support a different
    frame for each wavelength. The position of the instrument frame is determined by the properties
    of its parent MultiFrameInstrument object (in fact, by the angle attributes of its
    DistantInstrument base class). It is assumed that the distance to the system is sufficiently
    large so that parallel projection can be used. Refer to the SingleInstrument class for more
    information on the extent and pixel size attributes offered by this class. */
class InstrumentFrame : public SimulationItem
{
    ITEM_CONCRETE(InstrumentFrame, SimulationItem, "a frame in the multi-frame instrument")

    PROPERTY_DOUBLE(fieldOfViewX, "the total field of view in the horizontal direction")
        ATTRIBUTE_QUANTITY(fieldOfViewX, "length")
        ATTRIBUTE_MIN_VALUE(fieldOfViewX, "]0")

    PROPERTY_INT(numPixelsX, "the number of pixels in the horizontal direction")
        ATTRIBUTE_MIN_VALUE(numPixelsX, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsX, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsX, "250")

    PROPERTY_DOUBLE(centerX, "the center of the frame in the horizontal direction")
        ATTRIBUTE_QUANTITY(centerX, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerX, "0")

    PROPERTY_DOUBLE(fieldOfViewY, "the total field of view in the vertical direction")
        ATTRIBUTE_QUANTITY(fieldOfViewY, "length")
        ATTRIBUTE_MIN_VALUE(fieldOfViewY, "]0")

    PROPERTY_INT(numPixelsY, "the number of pixels in the vertical direction")
        ATTRIBUTE_MIN_VALUE(numPixelsY, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsY, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsY, "250")

    PROPERTY_DOUBLE(centerY, "the center of the frame in the vertical direction")
        ATTRIBUTE_QUANTITY(centerY, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerY, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function performs setup for the instrument frame. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

private:
    /** This function returns the index of the spatial pixel on the detector that will be hit by a
        photon package, or -1 if the photon package does not hit the detector. It operates
        similarly to SingleFrameInstrument::pixelondetector(). */
    int pixelOnDetector(const PhotonPackage* pp) const;

public:
    /** This function simulates the detection of a photon package by the instrument frame. It
        operates similarly to SimpleInstrument::detect(), but for a single wavelength. If the
        parent multi-frame instrument has the writeTotal flag turned on, this function records the
        total flux. If the writeStellarComps flag is turned on, this function records the flux for
        each stellar component seperately. */
    void detect(PhotonPackage* pp);

    /** This function properly calibrates and outputs the instrument data. It operates similarly to
        SimpleInstrument::write(), but for the single wavelength specified through its wavelength
        index \f$\ell\f$. If the parent multi-frame instrument has the writeTotal flag turned on,
        this function outputs the total flux in a correspondingly named file. If the parent
        multi-frame instrument has the writeStellarComps flag turned on, this function writes the
        flux for each stellar component in a seperate output file, with a name that includes the
        stellar component index. In all cases, the name of each output file includes the wavelength
        index. */
    void calibrateAndWriteData(int ell);

private:
    /** This private function properly calibrates and outputs the instrument data. It is invoked
        from the public calibrateAndWriteData() function. */
    void calibrateAndWriteDataFrames(int ell, const vector<Array*>& farrays, const vector<string>& fnames);

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _fovxp{_fieldOfViewX};
    const int& _Nxp{_numPixelsX};
    const double& _xpc{_centerX};
    const double& _fovyp{_fieldOfViewY};
    const int& _Nyp{_numPixelsY};
    const double& _ypc{_centerY};

    // data members derived from the published attributes during setup
    size_t _Nframep{0}; // number of pixels in a frame; size_t so that size and index calculations happen in 64 bit
    double _xpmin{0.};
    double _xpmax{0.};
    double _xpsiz{0.};
    double _ypmin{0.};
    double _ypmax{0.};
    double _ypsiz{0.};

    // data members copied from the parent multi-frame instrument during setup
    MultiFrameInstrument* _instrument{nullptr};
    bool _writeTotal{false};
    bool _writeStellarComps{false};
    double _distance{0.};
    double _cosphi{0.}, _sinphi{0.};
    double _costheta{0.}, _sintheta{0.};
    double _cospa{0.}, _sinpa{0.};

    // total flux per pixel
    Array _ftotv;
    ArrayTable<2> _fcompvv;
};

////////////////////////////////////////////////////////////////////

#endif
