/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustMassDustCompNormalization.hpp"

////////////////////////////////////////////////////////////////////

double DustMassDustCompNormalization::normalizationFactor(const Geometry* /*geom*/, const DustMix* /*mix*/) const
{
    return _dustMass; // geometry always has total mass of 1
}

//////////////////////////////////////////////////////////////////////