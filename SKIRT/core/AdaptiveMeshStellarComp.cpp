/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshStellarComp.hpp"
#include "AdaptiveMesh.hpp"
#include "Constants.hpp"
#include "BruzualCharlotSEDFamily.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "PhotonPackage.hpp"
#include "Random.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

AdaptiveMeshStellarComp::~AdaptiveMeshStellarComp()
{
    delete _mesh;
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshStellarComp::setupSelfBefore()
{
    BoxStellarComp::setupSelfBefore();

    // cache the random generator
    _random = find<Random>();

    // import the adaptive mesh
    _mesh = new AdaptiveMesh(_adaptiveMeshFile, vector<int>({_densityIndex, _metallicityIndex, _ageIndex}), extent());
    find<Log>()->info("Adaptive mesh data was successfully imported: " + std::to_string(_mesh->numCells()) + " cells.");

    // construct the library of SED models
    BruzualCharlotSEDFamily bc(this);

    find<Log>()->info("Filling the vectors with the SEDs of the cells... ");

    // local constants for units
    const double pc = Constants::pc();
    const double pc3 = pc*pc*pc;

    // the sizes of our vectors
    int Nlambda = find<WavelengthGrid>()->numWavelengths();
    int Ncells = _mesh->numCells();

    // construct a temporary matrix Lv with the luminosity of each cell at each wavelength
    // and also the permanent vector _Ltotv with the total luminosity for every wavelength bin
    ArrayTable<2> Lvv(Nlambda,Ncells);
    _Ltotv.resize(Nlambda);
    for (int m=0; m<Ncells; m++)
    {
        double rho = _mesh->value(_densityIndex, m);    // density in Msun / pc^3
        double V = _mesh->volume(m);                    // volume in m^3
        double M = rho * ( V/pc3 );                     // mass in Msun
        double Z = _mesh->value(_metallicityIndex, m);  // metallicity as dimensionless fraction
        double t = _mesh->value(_ageIndex, m);          // age in years

        const Array& Lv = bc.luminosities(M,Z,t);
        for (int ell=0; ell<Nlambda; ell++)
        {
            Lvv[ell][m] = Lv[ell];
            _Ltotv[ell] += Lv[ell];
        }
    }

    // construct the permanent vectors _Xvv with the normalized cumulative luminosities (per wavelength bin)
    _Xvv.resize(Nlambda,0);
    for (int ell=0; ell<Nlambda; ell++)
    {
        NR::cdf(_Xvv[ell], Lvv[ell]);
    }
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshStellarComp::luminosity(int ell) const
{
    return _Ltotv[ell];
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshStellarComp::launch(PhotonPackage* pp, int ell, double L) const
{
    int m = NR::locateClip(_Xvv[ell], _random->uniform());
    Position bfr = _mesh->randomPosition(_random, m);
    Direction bfk = _random->direction();
    pp->launch(L,ell,bfr,bfk);
}

//////////////////////////////////////////////////////////////////////
