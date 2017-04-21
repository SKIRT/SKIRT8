/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SED.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "StringUtils.hpp"
#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

void SED::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    size_t Nlambda = find<WavelengthGrid>()->numWavelengths();
    if (_Lv.size() != Nlambda) throw FATALERROR("The luminosities in the SED have not been properly set");
}

//////////////////////////////////////////////////////////////////////

double SED::luminosity(int ell) const
{
    return _Lv[ell];
}

//////////////////////////////////////////////////////////////////////

const Array& SED::luminosities() const
{
    return _Lv;
}

//////////////////////////////////////////////////////////////////////

void SED::setLuminosities(const Array &Lv)
{
    // copy luminosities
    _Lv = Lv;

    // normalize luminosities to unity
    double sum = _Lv.sum();
    if (sum<=0) throw FATALERROR("The total luminosity in the SED is zero or negative ("
                                 + StringUtils::toString(sum) + ")");
    _Lv /= sum;
}

//////////////////////////////////////////////////////////////////////

void SED::setEmissivities(const Array& jv)
{
    WavelengthGrid* lambdagrid = find<WavelengthGrid>();
    setLuminosities(jv * lambdagrid->dlambdav());
}

//////////////////////////////////////////////////////////////////////

void SED::setEmissivities(const Array& lambdav, const Array& jv)
{
    setEmissivities(NR::resample<NR::interpolateLogLog>(find<WavelengthGrid>()->lambdav(), lambdav, jv));
}

//////////////////////////////////////////////////////////////////////
