/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DraineGraphiteGrainComposition.hpp"
#include "Constants.hpp"
#include "SpecialFunctions.hpp"

//////////////////////////////////////////////////////////////////////

DraineGraphiteGrainComposition::DraineGraphiteGrainComposition(SimulationItem *parent)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

namespace
{
    // function f_n(x) as defined in equation (10) of Draine & Li (2001)
    // expressed in terms of the Debye function (perform substitution t = y/x)
    double fn(int n, double x)
    {
        return SpecialFunctions::DebyeD(n, 1./x) * x / (n*n);
    }
}

//////////////////////////////////////////////////////////////////////

double DraineGraphiteGrainComposition::enthalpyFunction(double T)
{
    // equation (9) of Draine & Li (2001), integrated to obtain internal energy rather than heat capacity,
    // provides the enthalpy for a single atom (approximating N-2 in the equation by N)
    const double TD1 = 863.;
    const double TD2 = 2504.;
    double h = Constants::k() * ( TD1*fn(2,T/TD1) + 2.*TD2*fn(2,T/TD2) );

    // convert to enthalpy per kg
    const double Mcarbon = 12.015;     // mass of carbon in units of M_p
    return h / (Mcarbon * Constants::Mproton());
}

//////////////////////////////////////////////////////////////////////

void DraineGraphiteGrainComposition::setupSelfBefore()
{
    GrainComposition::setupSelfBefore();

    loadOpticalGrid(true, "GrainComposition/Draine/Gra_81.dat", true, false, false, false);
    calculateEnthalpyGrid(enthalpyFunction);
    setBulkDensity(2.24e3);
}

//////////////////////////////////////////////////////////////////////

string DraineGraphiteGrainComposition::name() const
{
    return "Draine_Graphite";
}

//////////////////////////////////////////////////////////////////////
