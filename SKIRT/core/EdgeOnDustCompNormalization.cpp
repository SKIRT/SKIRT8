/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EdgeOnDustCompNormalization.hpp"
#include "AxGeometry.hpp"
#include "DustMix.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

double EdgeOnDustCompNormalization::normalizationFactor(const Geometry* geom, const DustMix* mix) const
{
    const AxGeometry* axgeom = dynamic_cast<const AxGeometry*>(geom);
    if (!axgeom) throw FATALERROR("Geometry is not axisymmetric");
    double sigma = axgeom->SigmaR();
    if (sigma<=0.) throw FATALERROR("Can't normalize dust mass for geometry with zero edge-on surface density");
    return _opticalDepth / ( sigma * mix->kappaext(_wavelength) );
}

//////////////////////////////////////////////////////////////////////
