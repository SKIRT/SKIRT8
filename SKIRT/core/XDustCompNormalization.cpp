/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XDustCompNormalization.hpp"
#include "DustMix.hpp"
#include "FatalError.hpp"
#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

double XDustCompNormalization::normalizationFactor(const Geometry* geom, const DustMix* mix) const
{
    double sigma = geom->SigmaX();
    if (sigma<=0.) throw FATALERROR("Can't normalize dust mass for geometry with zero X-axis surface density");
    return _opticalDepth / ( sigma * mix->kappaext(_wavelength) );
}

//////////////////////////////////////////////////////////////////////
