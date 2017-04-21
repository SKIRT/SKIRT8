/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LumFit2.hpp"
#include "Image.hpp"

////////////////////////////////////////////////////////////////////

void LumFit2::setMinLumA(double value)
{
    _minLumA = value;
}

////////////////////////////////////////////////////////////////////

void LumFit2::setMaxLumA(double value)
{
    _maxLumA = value;
}

////////////////////////////////////////////////////////////////////

void LumFit2::setMinLumB(double value)
{
    _minLumB = value;
}

////////////////////////////////////////////////////////////////////

void LumFit2::setMaxLumB(double value)
{
    _maxLumB = value;
}

////////////////////////////////////////////////////////////////////

void LumFit2::optimize(const Image& refframe, Image& frameA, Image& frameB,
                       double& lumA, double& lumB, double& chi2)
{
    _ref = &refframe;

    double simplex[3][3];

    // Alpha, Beta, Gamma and Delta are simplex optimization parameters
    double Alpha = 1;
    double Beta = 0.5;
    double Gamma = 2;
    double Delta = 0.5;
    initialize(frameA, frameB, simplex);

    for (int i=0; i<200; i++)
    {
        // Store simplex to check recurrency
        double previous_simpl[3][3];
        for (int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++) previous_simpl[j][k] = simplex[j][k];
        }

        // Set center and reflected point
        double center[3],refl[3];
        setCenterReflected(frameA, frameB, simplex, center, refl, i, Alpha);

        // Determine if the simplex needs reflection, expansion or contraction
        if (simplex[2][0] <= refl[2] && refl[2] < simplex[2][1]) place(frameA, frameB, simplex, refl[0],refl[1]);
        else
        {
            if (refl[2] < simplex[2][0]) expand(frameA, frameB, simplex, center,refl,i,Gamma);
            else if (refl[2] >= simplex[2][1]) contract(frameA, frameB, simplex, center,refl,Beta,Delta);
        }

        // Recurrency corrections
        if (previous_simpl[0][0] == simplex[0][0] && previous_simpl[0][1] == simplex[1][0]
           && previous_simpl[1][0] == simplex[0][1] && previous_simpl[1][1] == simplex[1][1]
           && previous_simpl[2][0] == simplex[0][2] && previous_simpl[2][1] == simplex[1][2])
        {
            place( frameA, frameB, simplex, (simplex[0][2] + simplex[0][1] + simplex[0][0]) / 3,
                    (simplex[1][2] + simplex[1][1] + simplex[1][0]) / 3);
        }

        // End loop if there is hardly any improvement
        double x_diff = abs(simplex[0][0] - simplex[0][1]) + abs(simplex[0][0] - simplex[0][2]);
        double y_diff = abs(simplex[1][0] - simplex[1][1]) + abs(simplex[1][0] - simplex[1][2]);
        if (x_diff <= 1e-6 && y_diff <= 1e-6) i=200;
    }

    lumA = simplex[0][0];
    lumB = simplex[1][0];
    chi2 = simplex[2][0];
}

////////////////////////////////////////////////////////////////////

bool LumFit2::inSimplex(double simplex[3][3], double value, int x_y ) const
{
    bool present = false;
    for (int i=0; i<3; i++)
        if(simplex[x_y][i]==value) present = true;

    return present;
}

////////////////////////////////////////////////////////////////////

double LumFit2::function(Image& frameA, Image& frameB, double x, double y)
{
    double chi = 0;
    int arraysize = frameA.size();

    for (int m = 0; m < arraysize; m++)
    {
        double total_sim = x * frameA[m] + y * frameB[m];
        double sigma = sqrt( abs((*_ref)[m]) + total_sim);
        if ((*_ref)[m]==0)
        {
            frameA[m] = 0;
            frameB[m] = 0;
            total_sim = 0;
            sigma = 0;
        }
        else
        {
            chi += pow( ((*_ref)[m] - total_sim) / sigma,2);
        }
    }
    return chi;
}

////////////////////////////////////////////////////////////////////

void LumFit2::contract(Image& frameA, Image& frameB,
                       double simplex [3][3], double center[], double refl[], double beta, double delta)
{
    double point[3];

    if (simplex[2][1] <= refl[2] && refl[2] < simplex[2][2])
    {
        point[0] = center[0] + beta * (refl[0]-center[0]);
        point[1] = center[1] + beta * (refl[1]-center[1]);
        point[2] = function(frameA, frameB, point[0], point[1]);

        if (point[2] <= refl[2]) place(frameA, frameB, simplex, point[0], point[1]);
        else shrink(frameA, frameB, simplex, delta);
    }
    else
    {
        if (refl[2]>=simplex[2][2])
        {
            point[0] = center[0] + beta * (simplex[0][2] - center[0]);
            point[1] = center[1] + beta * (simplex[1][2] - center[1]);
            point[2] = function(frameA, frameB, point[0], point[1]);

            if (point[2]<simplex[2][2]) place(frameA, frameB, simplex, point[0], point[1]);
            else shrink(frameA, frameB, simplex, delta);
        }
        else shrink(frameA, frameB, simplex, delta);
    }
}

////////////////////////////////////////////////////////////////////

void LumFit2::expand(Image& frameA, Image& frameB,
                     double simplex[3][3], double center[], double refl[], int counter, double gamma)
{
    double point[3];
    point[0] = center[0] + gamma * (refl[0]-center[0]);
    point[1] = center[1] + gamma * (refl[1]-center[1]);
    nearEdgeCorrections(simplex, point, counter);
    point[2] = function(frameA, frameB, point[0], point[1]);

    if (point[2]<refl[2]) place(frameA, frameB, simplex, point[0], point[1]);
    else place(frameA, frameB, simplex, refl[0],refl[1]);
}

////////////////////////////////////////////////////////////////////

void LumFit2::initialize(Image& frameA, Image& frameB, double simplex[3][3])
{
    // Determine the initial simplex points
    double xle = _maxLumA - _minLumA;
    double yle = _maxLumB - _minLumB;
    double x_1 = _minLumA + 0.8*xle;
    double y_1 = _minLumB + 0.5*yle;
    double x_2 = _minLumA + 0.45*xle;
    double y_2 = _minLumB + 0.05*yle;
    double x_3 = _minLumA + 0.20*xle;
    double y_3 = _minLumB + 0.82*yle;
    double x_max = x_1;
    double y_max = y_1;

    // Determine the point with highest chi2
    if (function(frameA, frameB, x_1,y_1) > function(frameA, frameB, x_2,y_2) &&
        function(frameA, frameB, x_1,y_1) > function(frameA, frameB, x_3,y_3))
    {
        x_max = x_1;
        y_max = y_1;
    }
    if (function(frameA, frameB, x_2,y_2) > function(frameA, frameB, x_1,y_1) &&
        function(frameA, frameB, x_2,y_2) > function(frameA, frameB, x_3,y_3))
    {
        x_max = x_2;
        y_max = y_2;
    }
    if (function(frameA, frameB, x_3,y_3) > function(frameA, frameB, x_2,y_2) &&
        function(frameA, frameB, x_3,y_3) > function(frameA, frameB, x_1,y_1))
    {
        x_max = x_3;
        y_max = y_3;
    }

    // Fill the simplex with the worst value
    for (int l=0; l<3; l++)
    {
        simplex[0][l] = x_max;
        simplex[1][l] = y_max;
        simplex[2][l] = function( frameA, frameB, x_max,y_max);
    }


    // Add the other values, ranked from best to worst
    place(frameA, frameB, simplex, x_1,y_1);
    place(frameA, frameB, simplex, x_2,y_2);
    place(frameA, frameB, simplex, x_3,y_3);
}

////////////////////////////////////////////////////////////////////

void LumFit2::nearEdgeCorrections(double simplex[3][3], double Dpoint[], int counter) const
{
    double xle = _maxLumA - _minLumA;
    double yle = _maxLumB - _minLumB;
    int mod = counter/2;

    if (Dpoint[0]> _maxLumA)
    {
        if (inSimplex(simplex, _maxLumA,0)) Dpoint[0] = _maxLumA - (1+mod) * 0.01 * xle;
        else Dpoint[0] = _maxLumA;
    }
    if (Dpoint[0]< _minLumA)
    {
        if (inSimplex(simplex, _minLumA,0)) Dpoint[0] = _minLumA + (1+mod) * 0.01 * xle;
        else Dpoint[0] = _minLumA;
    }
    if (Dpoint[1]> _maxLumB)
    {
        if (inSimplex(simplex, _maxLumB,1)) Dpoint[1] = _maxLumB - (1+mod) * 0.01 * yle;
        else Dpoint[1] = _maxLumB;
    }
    if (Dpoint[1]<_minLumB)
    {
        if (inSimplex(simplex, _minLumB,1)) Dpoint[1] = _minLumB + (1+mod) * 0.01 * yle;
        else Dpoint[1] = _minLumB;
    }
}

////////////////////////////////////////////////////////////////////

void LumFit2::place(Image& frameA, Image& frameB, double simplex[3][3], double x, double y)
{
    for (int i=0; i<3; i++)
    {
        if (function(frameA, frameB, x,y) <= simplex[2][i])
        {
            for (int j=2; j>i; j--)
            {
                simplex[0][j] = simplex[0][j-1];
                simplex[1][j] = simplex[1][j-1];
                simplex[2][j] = simplex[2][j-1];
            }
            simplex[0][i] = x;
            simplex[1][i] = y;
            simplex[2][i] = function(frameA, frameB, x,y);
            i = 4;
        }
    }
}

////////////////////////////////////////////////////////////////////

void LumFit2::setCenterReflected(Image& frameA, Image& frameB, double simplex[3][3], double center[],
                                 double reflected[], int counter, double alpha)
{
    double averx = 0;
    double avery = 0;

    for (int i=0; i<2; i++)
    {
        averx += simplex[0][i];
        avery += simplex[1][i];
    }
    center[0] = averx/2;
    center[1] = avery/2;
    center[2] = function(frameA, frameB, center[0],center[1]);

    reflected[0] = center[0] + alpha * (center[0] - simplex[0][2]);
    reflected[1] = center[1] + alpha * (center[1] - simplex[1][2]);
    nearEdgeCorrections(simplex, reflected, counter);
    reflected[2] = function(frameA, frameB, reflected[0], reflected[1]);
}

////////////////////////////////////////////////////////////////////

void LumFit2::shrink(Image& frameA, Image& frameB, double simplex[3][3], double delta)
{
    double x_1 = simplex[0][0] + delta * (simplex[0][1] - simplex[0][0]);
    double x_2 = simplex[0][0] + delta * (simplex[0][2] - simplex[0][0]);
    double y_1 = simplex[1][0] + delta * (simplex[1][1] - simplex[1][0]);
    double y_2 = simplex[1][0] + delta * (simplex[1][2] - simplex[1][0]);

    place(frameA, frameB, simplex, x_1,y_1);
    place(frameA, frameB, simplex, x_2,y_2);
}

////////////////////////////////////////////////////////////////////
