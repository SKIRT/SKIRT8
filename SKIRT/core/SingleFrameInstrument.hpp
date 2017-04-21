/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SINGLEFRAMEINSTRUMENT_HPP
#define SINGLEFRAMEINSTRUMENT_HPP

#include "DistantInstrument.hpp"

////////////////////////////////////////////////////////////////////

/** SingleFrameInstrument is an abstract class representing an integral field spectrograph, i.e. a
    three-dimensional detector with two spatial and one wavelength dimension. Depending on the
    subclass, this 3D detector can be subdivided into different subdetectors. In every dimension,
    the detector is specified by a fixed grid. The wavelength grid is determined by the
    simulation's wavelength grid; the field-of-view, number of pixels and image centre in the
    spatial dimensions can be set separately. The spatial grid consists of a linear grid of \f$N_x
    \times N_y\f$ rectangular pixels. Each pixel covers a surface \f$\Delta x \times \Delta y\f$
    with \f[ \Delta x = \frac{{\text{FOV}}_x}{N_x}, \Delta_y = \frac{{\text{FOV}}_y}{N_y} \f] where
    \f${\text{FOV}}_x\f$ and \f${\text{FOV}}_y\f$ are obviously the field-of-view in the X and Y
    direction, respectively. If we denote the centre of the image as \f$(x_c,y_c)\f$, we have \f[
    \begin{split} x_{\text{min}} &= x_c - \tfrac12\,{\text{FOV}}_x \\ x_{\text{max}} &= x_c +
    \tfrac12\,{\text{FOV}}_x \end{split} \f] and similar for the Y direction. The first pixel of
    the grid corresponds to \f[ \begin{split} x_{\text{min}} \leq &x < x_{\text{min}} + \Delta x \\
    y_{\text{min}} \leq &y < y_{\text{min}} + \Delta y \end{split} \f] The position of the
    observing instrument is determined by the properties of the DistantInstrument base class. It is
    assumed that the distance to the system is sufficiently large so that parallel projection can
    be used. */
class SingleFrameInstrument : public DistantInstrument
{
    ITEM_ABSTRACT(SingleFrameInstrument, DistantInstrument, "a single-frame instrument")

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
    /** This function performs setup for the instrument. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

protected:
    /** This convenience function returns the index of the spatial pixel on the detector that will
        be hit by a photon package, or -1 if the photon package does not hit the detector. Given
        the position \f${\boldsymbol{x}}=(x,y,z)\f$ of the last emission or scattering event of the
        photon package, the direction \f${\boldsymbol{k}}_{\text{obs}} = (\theta,\varphi)\f$
        towards the observer, and the position angle \f$\omega\f$ of the instrument, the impact
        coordinates \f$(x_{\text{p}},y_{\text{p}})\f$ are given by the following Euler-like
        transformation, where \f$z_{\text{p}}\f$ is ignored: \f[ \begin{bmatrix}x_{\text{p}} \\
        y_{\text{p}} \\ z_{\text{p}} \end{bmatrix} = \begin{bmatrix}\cos\omega & -\sin\omega & 0\\
        \sin\omega & \cos\omega & 0\\ 0 & 0 & 1 \end{bmatrix} \begin{bmatrix} 0 & 1 & 0 \\ -1 & 0 &
        0 \\ 0 & 0 & 1 \end{bmatrix} \begin{bmatrix} \cos\theta & 0 & -\sin\theta \\ 0 & 1 & 0 \\
        \sin\theta & 0 & \cos\theta \end{bmatrix} \begin{bmatrix}\cos\varphi & -\sin(-\varphi) & 0
        \\ \sin(-\varphi) & \cos\varphi & 0 \\ 0 & 0 & 1 \end{bmatrix}
        \begin{bmatrix}x\\y\\z\end{bmatrix} \f] In other words, the originating position is rotated
        about the Z-axis over the azimuth angle \f$\varphi\f$ (with a minus sign because the
        observer is looking towards the center rather than along the specified direction), then
        rotated about the new Y-axis over the inclination angle \f$\theta\f$, and finally rotated
        about the new Z-axis over the position angle \f$\omega\f$ reduced by 90 degrees (this
        constant transformation over -90 degrees is represented above as a separate matrix). The 90
        degree correction on the position angle is introduced so that it would be more natural to
        specify this angle; in most cases it can be left to its default value of 0. Given these
        impact coordinates, the pixel indices \f$i\f$ and \f$j\f$ are determined as \f[
        \begin{split} i &= \frac{{\text{floor}}(x_{\text{p}}-x_{\text{min}})}{\Delta x} \\ j &=
        \frac{{\text{floor}}(y_{\text{p}}-y_{\text{max}})}{\Delta y} \end{split} \f] where
        \f${\text{floor}}(z)\f$ is an operator that returns the largest integer that is not greater
        than \f$y\f$. The spatial pixel number \f$l\f$ is then determined as \f$l=i+j\,N_x\f$,
        asuming \f$i\f$ and \f$j\f$ are indeed within the detector range. */
    int pixelOnDetector(const PhotonPackage* pp) const;

    /** This convenience function calibrates one or more luminosity data cubes gathered by a
        DistantInstrument subclass and outputs each data cube as a FITS file. The incoming data is
        organized as a list of data arrays and a second list of corresponding human-readable names.
        Each array in the first list is a 3D data cube containing a limunosity value per pixel in
        the 2D frame and per wavelength in the simulation's wavelength grid. The strings in the
        second list are used as part of the output file name for the corresponding data cube using
        the pattern <tt>prefix_instrument_array.fits</tt>. The two lists must have the same size;
        if they are empty no output is generated. The calibration performed by this function takes
        care of the conversion from bolometric luminosity units to surface brightness units. The
        unit in which the surface brightness is written depends on the global units choice, but
        typically it is in \f$\text{W}\,\text{m}^{-2}\,\text{arcsec}^{-2}\f$. The calibration is
        performed in-place in the arrays, so the incoming data is overwritten. */
    void calibrateAndWriteDataCubes(const vector<Array*>& farrays, const vector<string>& fnames);

    //======================== Data Members ========================

protected:
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
};

////////////////////////////////////////////////////////////////////

#endif
