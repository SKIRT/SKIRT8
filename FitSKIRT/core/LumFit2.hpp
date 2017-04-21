/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LUMFIT2_HPP
#define LUMFIT2_HPP

#include "Basics.hpp"
class Image;

////////////////////////////////////////////////////////////////////

/** The LumFit2 class contains all the information used to optimize a 2D luminosity problem.
    The optimization is performed using the simplex algorithm described in
    http://www.scholarpedia.org/article/Nelder-Mead_algorithm. */
class LumFit2
{
    //======================== Setters and Getters =======================

public:
    /** Sets the minimal value for the first luminosity multiplicator. */
    void setMinLumA(double value);

    /** Sets the maximal value for the first luminosity multiplicator. */
    void setMaxLumA(double value);

    /** Sets the minimal value for the second luminosity multiplicator. */
    void setMinLumB(double value);

    /** Sets the maximal value for the second luminosity multiplicator. */
    void setMaxLumB(double value);

    //======================== Other Functions =======================

public:
    /** This function returns the luminosity multiplicators that best match the weighted sum of the
        input frames to the reference frame, and the corresponding \f$\chi^2\f$ value. The input
        frames are adjusted to contain the same mask as the reference frame. */
    void optimize(const Image& refframe, Image& frameA, Image& frameB,
                  double& lumA, double& lumB, double& chi2);

private:
    /** This function determines if the x- or y-value is present in the simplex. */
    bool inSimplex(double simplex[3][3], double value, int x_y ) const;

    /** This function determines the \f$\chi^2\f$ value for certain luminosity values. */
    double function(Image& frameA, Image& frameB, double x, double y);

    /** This function determines if the simplex needs to be contracted or shrunk. */
    void contract(Image& frameA, Image& frameB,
                  double simplex[3][3], double center[], double refl[], double beta, double delta);

    /** This function determines if the simplex needs to be reflected or expanded. */
    void expand(Image& frameA, Image& frameB,
                double simplex[3][3], double center[], double refl[], int counter, double Gamma);

    /** This function determines the initial simplex and sorts it. */
    void initialize(Image& frameA, Image& frameB, double simplex[3][3]);

    /** This function checks if the simplex goes out of bound and corrects if necessary. The
        counter is used to make sure the corrections are not applied twice so the simplex collapses
        to a line. */
    void nearEdgeCorrections(double simplex[3][3], double Dpoint[], int counter) const;

    /** This function places the two values in the simplex in the correct place. The simplex is
        sorted from lowest \f$\chi^2\f$ value to highest. */
    void place(Image& frameA, Image& frameB, double simplex[3][3], double x, double y);

    /** This function determines and sets the center and reflected point. */
    void setCenterReflected(Image& frameA, Image& frameB,
                            double simplex[3][3], double center[], double reflected[], int counter, double alpha);

    /** This function determines if the simplex needs to be shrunk. */
    void shrink(Image& frameA, Image& frameB, double simplex[3][3], double delta);

    //======================== Data Members ========================

private:
    double _minLumA{0.};
    double _maxLumA{0.};
    double _minLumB{0.};
    double _maxLumB{0.};
    const Image* _ref{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
