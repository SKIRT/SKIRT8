/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FaceOnDustCompNormalization.hpp"
#include "AxGeometry.hpp"
#include "DustMix.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

double FaceOnDustCompNormalization::normalizationFactor(const Geometry* geom, const DustMix* mix) const
{
    const AxGeometry* axgeom = dynamic_cast<const AxGeometry*>(geom);
    if (!axgeom) throw FATALERROR("Geometry is not axisymmetric");
    double sigma = axgeom->SigmaZ();
    if (sigma<=0.) throw FATALERROR("Can't normalize dust mass for geometry with zero face-on surface density");
    return _opticalDepth / ( sigma * mix->kappaext(_wavelength) );
}

//////////////////////////////////////////////////////////////////////
