/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "YDustCompNormalization.hpp"
#include "DustMix.hpp"
#include "FatalError.hpp"
#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

double YDustCompNormalization::normalizationFactor(const Geometry* geom, const DustMix* mix) const
{
    double sigma = geom->SigmaY();
    if (sigma<=0.) throw FATALERROR("Can't normalize dust mass for geometry with zero Y-axis surface density");
    return _opticalDepth / ( sigma * mix->kappaext(_wavelength) );
}

//////////////////////////////////////////////////////////////////////
