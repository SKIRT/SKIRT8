/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ZDustCompNormalization.hpp"
#include "DustMix.hpp"
#include "FatalError.hpp"
#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

double ZDustCompNormalization::normalizationFactor(const Geometry* geom, const DustMix* mix) const
{
    double sigma = geom->SigmaZ();
    if (sigma<=0.) throw FATALERROR("Can't normalize dust mass for geometry with zero Z-axis surface density");
    return _opticalDepth / ( sigma * mix->kappaext(_wavelength) );
}

//////////////////////////////////////////////////////////////////////
