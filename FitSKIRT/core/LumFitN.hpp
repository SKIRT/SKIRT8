/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LUMFITN_HPP
#define LUMFITN_HPP

#include "Basics.hpp"
class Image;

////////////////////////////////////////////////////////////////////

/** The LumFitN class contains all the information needed to optimize an \f$N\f$-dimensional
    luminosity problem. It uses genetic algorithms provide by GALib to perform the optimization. */
class LumFitN
{
    //======================== Setters and Getters =======================

public:
    /** Sets the list of \f$N\f$ minimum luminosity multiplicators, one for each frame in the
        multi-dimensional luminosity problem. */
    void setMinLuminosities(const vector<double>& value);

    /** Sets the list of \f$N\f$ maximum luminosity multiplicators, one for each frame in the
        multi-dimensional luminosity problem. */
    void setMaxLuminosities(const vector<double>& value);

    //======================== Other Functions =======================

public:
    /** This function returns the \f$N\f$ luminosity multiplicators that best match the weighted
        sum of the input frames to the reference frame, and the corresponding \f$\chi^2\f$ value.
        The input frames are adjusted to contain the same mask as the reference frame. */
    void optimize(const Image& refframe, vector<Image>& frames, vector<double>& luminosities, double& chi2);

    //======================== Data Members ========================

private:
    vector<double> _minLum;
    vector<double> _maxLum;
};

////////////////////////////////////////////////////////////////////

#endif
