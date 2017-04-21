/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LumFit1.hpp"
#include "Image.hpp"

////////////////////////////////////////////////////////////////////

void LumFit1::setMinLum(double value)
{
    _minLum = value;
}

////////////////////////////////////////////////////////////////////

void LumFit1::setMaxLum(double value)
{
    _maxLum = value;
}

////////////////////////////////////////////////////////////////////

void LumFit1::optimize(const Image& refframe, Image& frame, double& lum, double& chi2)
{
    _ref = &refframe;
    Image sim = frame;
    double GOLD = 0.3819660113;

    double chi = 1;
    double lumvalue = 1;
    double prev_lumvalue=1;
    double a0 = _minLum;
    double b0 = _maxLum;

    // loop to constrain the best fitting luminosity by removing the worst boundary
    for (int j = 0; j < 300; j++)
    {
        double d = (b0-a0)*GOLD;
        double a1 = a0+d;
        double b1 = b0-d;
        double chia1 = function(sim, a1);
        double chib1 = function(sim, b1);

        if (chia1 < chib1)
        {
            b0=b1;
            lumvalue=a1;
            chi=chia1;
        }
        else
        {
            a0 = a1;
            lumvalue = b1;
            chi = chib1;
        }

        // condition to end fitting if the value does not change significantly
        if (abs(prev_lumvalue-lumvalue)/lumvalue <= 1e-8 && j>=20) j = 300;
        else prev_lumvalue = lumvalue;
    }
    chi2 = chi;
    lum = lumvalue;
    frame = sim;
}

////////////////////////////////////////////////////////////////////

double LumFit1::function(Image& frame, double x)
{
    double chi = 0;
    int arraysize = frame.size();

    // calculate the chi2 value and copy masked regions
    for (int m = 0; m < arraysize; m++)
    {
        double total_sim = x * frame[m];
        double sigma = sqrt( abs((*_ref)[m]) + total_sim);
        if ((*_ref)[m] == 0)
        {
            frame[m] = 0;
            total_sim = 0;
            sigma = 0;
        }
        else
        {
            chi += pow( ((*_ref)[m] - total_sim) / sigma, 2);
        }
    }
    return chi;
}

////////////////////////////////////////////////////////////////////
