/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SimpleOligoDustMix.hpp"
#include "FatalError.hpp"
#include "OligoWavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

void SimpleOligoDustMix::setupSelfBefore()
{
    DustMix::setupSelfBefore();

    // Verify that the wavelength grid (and thus the simulation) is of the Oligo type
    OligoWavelengthGrid* lambdagrid = find<OligoWavelengthGrid>();

    // Verify that the number of luminosities equals the number of wavelengths
    int Nlambda = lambdagrid->numWavelengths();
    if (_opacities.size() != static_cast<size_t>(Nlambda))
        throw FATALERROR("The number of extinction coefficients differs from the number of wavelengths");
    if (_albedos.size() != static_cast<size_t>(Nlambda))
        throw FATALERROR("The number of albedos differs from the number of wavelengths");
    if (_asymmetryParameters.size() != static_cast<size_t>(Nlambda))
        throw FATALERROR("The number of asymmetry parameters differs from the number of wavelengths");

    // Create temporary vectors with the absorption and scattering coefficients
    Array kappaabsv(Nlambda);
    Array kappascav(Nlambda);
    Array asymmparv(Nlambda);
    for (int ell=0; ell<Nlambda; ell++)
    {
        kappaabsv[ell] = _opacities[ell]*(1.-_albedos[ell]);
        kappascav[ell] = _opacities[ell]*_albedos[ell];
        asymmparv[ell] = _asymmetryParameters[ell];
    }

    // Add a dust population with these properties, providing 1 kg for the dust mass to fix the units to m2/kg
    addPopulation(1., kappaabsv, kappascav, asymmparv);
}

//////////////////////////////////////////////////////////////////////
