/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LUMFIT1_HPP
#define LUMFIT1_HPP

#include "Basics.hpp"
class Image;

////////////////////////////////////////////////////////////////////

/** The LumFit1 class contains all the information used to optimize a 1D luminosity problem.
    It uses the golden section method to constrain the single variable parameter.
    http://www.aip.de/groups/soe/local/numres/bookcpdf/c10-1.pdf */
class LumFit1
{
    //======================== Setters and Getters =======================

public:
    /** Sets the minimum luminosity multiplicator. */
    void setMinLum(double value);

    /** Sets the maximum luminosity multiplicator. */
    void setMaxLum(double value);

    //======================== Other Functions =======================

public:
    /** This function returns the luminosity multiplicator that best matches the input frame to the
        reference frame, and the corresponding \f$\chi^2\f$ value. The input frame is adjusted to
        contain the same mask as the reference frame. */
    void optimize(const Image& refframe, Image& frame, double& lum, double& chi2);

private:
    /** This function determines the \f$\chi^2\f$ value for a certain luminosity value x. */
    double function(Image& frame, double x);

    //======================== Data Members ========================

private:
    double _minLum{0.};
    double _maxLum{0.};
    const Image* _ref{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
