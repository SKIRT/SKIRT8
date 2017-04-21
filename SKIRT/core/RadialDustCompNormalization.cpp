/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RadialDustCompNormalization.hpp"
#include "DustMix.hpp"
#include "FatalError.hpp"
#include "SpheGeometry.hpp"

////////////////////////////////////////////////////////////////////

double RadialDustCompNormalization::normalizationFactor(const Geometry* geom, const DustMix* mix) const
{
    const SpheGeometry* sphegeom = dynamic_cast<const SpheGeometry*>(geom);
    if (!sphegeom) throw FATALERROR("Geometry is not spherically symmetric");
    double sigma = sphegeom->Sigmar();
    if (sigma<=0.) throw FATALERROR("Can't normalize dust mass for geometry with zero radial surface density");
    return _opticalDepth / ( sigma * mix->kappaext(_wavelength) );
}

//////////////////////////////////////////////////////////////////////
