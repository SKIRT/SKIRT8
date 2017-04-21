/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParallelDataCube.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "ProcessAssigner.hpp"
#include "StringUtils.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

ParallelDataCube::ParallelDataCube()
    : _partialCube(std::make_shared<Array>())
{
}

////////////////////////////////////////////////////////////////////

void ParallelDataCube::initialize(string name, size_t Nframep, SimulationItem* item)
{
    Log* log = item->find<Log>();
    _Nframep = Nframep;
    _comm = item->find<PeerToPeerCommunicator>();

    WavelengthGrid* wg = item->find<WavelengthGrid>();

    if(_comm->dataParallel())
    {
        _wavelengthAssigner = wg->assigner();
        _Nlambda = _wavelengthAssigner->assigned();
        log->info(name + " data cube is distributed; local size is "
                  + std::to_string(_Nframep) + "x" + std::to_string(_Nlambda)
                  + " (" + StringUtils::toMemSizeString(_Nframep*_Nlambda*sizeof(double)) + ")");
    }
    else
    {
        _wavelengthAssigner = nullptr;
        _Nlambda = wg->numWavelengths();
        log->info(name + " data cube is not distributed; size is "
                  + std::to_string(_Nframep) + "x" + std::to_string(_Nlambda)
                  + " (" + StringUtils::toMemSizeString(_Nframep*_Nlambda*sizeof(double)) + ")");
    }
    _partialCube->resize(_Nlambda*_Nframep);
}

////////////////////////////////////////////////////////////////////

std::shared_ptr<Array> ParallelDataCube::constructCompleteCube()
{
    // partial cube of equal size as total cube
    if (!_wavelengthAssigner || !_comm->isMultiProc())
    {
        // sum the data to root
        _comm->sum(*_partialCube);

        // give a handle to the summed cube at the root, and a dummy for the other processes
        return _comm->isRoot() ? _partialCube : std::make_shared<Array>();
    }
    // total cube is bigger than partial cube
    else
    {
        // allocate space for complete cube at root process
        auto completeCube = std::make_shared<Array>();
        if (_comm->isRoot()) completeCube->resize(_wavelengthAssigner->total()*_Nframep);

        // displacements parameters for gatherw
        std::vector<std::vector<int>> displacements;
        displacements.reserve(_comm->size());
        for (int i=0; i<_comm->size(); i++) displacements.push_back(_wavelengthAssigner->indicesForRank(i));

        // gather complete cube
        _comm->gatherWithPattern(&(*_partialCube)[0], _Nlambda*_Nframep, &(*completeCube)[0], 0, _Nframep, displacements);
        return completeCube;
    }
}

////////////////////////////////////////////////////////////////////

double& ParallelDataCube::operator()(int ell, int pixel)
{
    if (!_wavelengthAssigner)
        return (*_partialCube)[ell*_Nframep + pixel];
    else
    {
        if (!_wavelengthAssigner->validIndex(ell))
            throw FATALERROR("Wrong wavelength for this process!");

        return (*_partialCube)[_wavelengthAssigner->relativeIndex(ell)*_Nframep + pixel];
    }
}

////////////////////////////////////////////////////////////////////

double ParallelDataCube::operator()(int ell, int pixel) const
{
    if(!_wavelengthAssigner) return (*_partialCube)[ell*_Nframep + pixel];
    else
    {
        if (!_wavelengthAssigner->validIndex(ell))
            throw FATALERROR("Wrong wavelength for this process!");

        return (*_partialCube)[_wavelengthAssigner->relativeIndex(ell)*_Nframep + pixel];
    }
}

////////////////////////////////////////////////////////////////////
