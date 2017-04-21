/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SphericalAdaptiveMesh.hpp"
#include "AdaptiveMeshFile.hpp"
#include "AdaptiveMeshNode.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

SphericalAdaptiveMesh::SphericalAdaptiveMesh(AdaptiveMeshFile* meshfile, const vector<int>& fieldIndices,
                                             double rin, double rout)
    : _rin(rin), _rout(rout)
{
    // open the data file
    meshfile->open();

    // create a list of indices (g) without duplicates, ignoring negative values
    // create a hash table to map field indices (g) to storage indices (s)
    vector<int> uniqueIndices;
    for (int g : fieldIndices)
    {
        if (g >= 0 && !_storageIndices.count(g))
        {
            _storageIndices.emplace(g, uniqueIndices.size());
            uniqueIndices.push_back(g);
        }
    }

    // reserve room for the required number of fields
    _fieldvalues.resize(fieldIndices.size());

    // construct the root node, and recursively all other nodes
    // this also fills the _fieldvalues and _leafnodes vectors
    _root = new AdaptiveMeshNode(Box(_rin, 0, 0, _rout, M_PI, 2*M_PI),
                                 uniqueIndices, meshfile, _leafnodes, _fieldvalues);
    _Ncells = _leafnodes.size();

    // verify that all data was read and close the file
    if (meshfile->read()) throw FATALERROR("Superfluous data in mesh data after all nodes were read");
    meshfile->close();

    // determine small value relative to the domain extent
    _eps = 1e-12 * _rout;

    // clear the density distribution info
    _Ndistribs = 0;
    _integratedDensity = 0;
}

////////////////////////////////////////////////////////////////////

void SphericalAdaptiveMesh::addDensityDistribution(int densityField, int densityMultiplierField, double densityFraction)
{
    // verify indices
    if (!_storageIndices.count(densityField)) throw FATALERROR("Density field index out of range");
    if (densityMultiplierField >= 0 && (!_storageIndices.count(densityMultiplierField)
                                        || densityMultiplierField==densityField))
        throw FATALERROR("Density multiplier field index out of range");

    // map field indices to storage indices
    densityField = _storageIndices.at(densityField);
    densityMultiplierField = _storageIndices.count(densityMultiplierField)
                             ? _storageIndices.at(densityMultiplierField) : -1;

    // store the information for this density distribution
    _densityFields.push_back(densityField);
    _densityMultiplierFields.push_back(densityMultiplierField);
    _densityFractions.push_back(densityFraction);
    _Ndistribs++;

    // update the integrated density (ignore cells with negative density)
    double integratedDensity = 0.0;
    for (int m=0; m<_Ncells; m++)
    {
        double density = _fieldvalues[densityField][m] * densityFraction;
        if (densityMultiplierField >= 0) density *= _fieldvalues[densityMultiplierField][m];
        if (density > 0) integratedDensity += density*volume(m);
    }
    _integratedDensityv.push_back(integratedDensity);
    _integratedDensity += integratedDensity;
}

////////////////////////////////////////////////////////////////////

SphericalAdaptiveMesh::~SphericalAdaptiveMesh()
{
    delete _root;
}

////////////////////////////////////////////////////////////////////

int SphericalAdaptiveMesh::numCells() const
{
    return _Ncells;
}

////////////////////////////////////////////////////////////////////

int SphericalAdaptiveMesh::cellIndex(Position bfr) const
{
    // convert from cartesian to spherical coordinates
    double r, theta, phi;
    bfr.spherical(r, theta, phi);
    if (phi < 0) phi += 2*M_PI;
    const AdaptiveMeshNode* node = _root->whichNode(Position(r, theta, phi, Position::CoordinateSystem::CARTESIAN) );
    return node ? node->cellIndex() : -1;
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::volume() const
{
    return 4*M_PI/3 * (_rout*_rout*_rout - _rin*_rin*_rin);
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::volume(int m) const
{
    if (m < 0 || m >= _Ncells) throw FATALERROR("Cell index out of range: " + std::to_string(m));
    Box extent = _leafnodes[m]->extent();
    double r1 = extent.xmin();
    double r2 = extent.xmax();
    double t1 = extent.ymin();
    double t2 = extent.ymax();
    double f1 = extent.zmin();
    double f2 = extent.zmax();
    return 1./3. * (r2*r2*r2-r1*r1*r1) * (cos(t1)-cos(t2)) * (f2-f1);
}

////////////////////////////////////////////////////////////////////

Position SphericalAdaptiveMesh::randomPosition(Random* random, int m) const
{
    if (m < 0 || m >= _Ncells) throw FATALERROR("Cell index out of range: " + std::to_string(m));
    Box extent = _leafnodes[m]->extent();
    double r1 = extent.xmin();
    double r2 = extent.xmax();
    double t1 = extent.ymin();
    double t2 = extent.ymax();
    double f1 = extent.zmin();
    double f2 = extent.zmax();
    double r1s = r1*r1;
    double r2s = r2*r2;
    double r = sqrt( r1s + (r2s-r1s)*random->uniform() );
    double t = t1 + (t2-t1)*random->uniform();
    double f = f1 + (f2-f1)*random->uniform();
    return Position(r,t,f,Position::CoordinateSystem::SPHERICAL);
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::value(int g, int m) const
{
    if (!_storageIndices.count(g)) throw FATALERROR("Field index out of range: " + std::to_string(g));
    if (m < 0 || m >= _Ncells) throw FATALERROR("Cell index out of range: " + std::to_string(m));
    int s = _storageIndices.at(g);
    return _fieldvalues[s][m];
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::value(int g, Position bfr) const
{
    if (!_storageIndices.count(g)) throw FATALERROR("Field index out of range: " + std::to_string(g));
    int s = _storageIndices.at(g);
    int m = cellIndex(bfr);
    return m>=0 ? _fieldvalues[s][m] : 0;
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::density(int h, int m) const
{
    if (!_integratedDensity) throw FATALERROR("There is no density field");
    if (h < 0 || h >= _Ndistribs) throw FATALERROR("Density distribution index out of range: " + std::to_string(h));
    if (m < 0 || m >= _Ncells) throw FATALERROR("Cell index out of range: " + std::to_string(m));

    double density = _fieldvalues[_densityFields[h]][m] * _densityFractions[h];
    if (_densityMultiplierFields[h] >= 0) density *= _fieldvalues[_densityMultiplierFields[h]][m];
    return density > 0 ? density : 0;
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::density(int h, Position bfr) const
{
    int m = cellIndex(bfr);
    return m>=0 ? density(h, m) : 0;
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::density(int m) const
{
    if (!_integratedDensity) throw FATALERROR("There is no density field");
    if (m < 0 || m >= _Ncells) throw FATALERROR("Cell index out of range: " + std::to_string(m));

    double result = 0;
    for (int h=0; h<_Ndistribs; h++)
    {
        double density = _fieldvalues[_densityFields[h]][m] * _densityFractions[h];
        if (_densityMultiplierFields[h] >= 0) density *= _fieldvalues[_densityMultiplierFields[h]][m];
        if (density > 0) result += density;
    }
    return result;
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::density(Position bfr) const
{
    int m = cellIndex(bfr);
    return m>=0 ? density(m) : 0;
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::integratedDensity(int h) const
{
    return _integratedDensityv[h];
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::integratedDensity() const
{
    return _integratedDensity;
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::SigmaX() const
{
    const int NSAMPLES = 10000;
    double sum = 0;
    double xmin = -_rout;
    double xmax = _rout;
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(xmin + k*(xmax-xmin)/NSAMPLES, _eps, _eps));
    }
    return (sum/NSAMPLES)*(xmax-xmin);
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::SigmaY() const
{
    const int NSAMPLES = 10000;
    double sum = 0;
    double ymin = -_rout;
    double ymax = _rout;
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(_eps, ymin + k*(ymax-ymin)/NSAMPLES, _eps));
    }
    return (sum/NSAMPLES)*(ymax-ymin);
}

////////////////////////////////////////////////////////////////////

double SphericalAdaptiveMesh::SigmaZ() const
{
    const int NSAMPLES = 10000;
    double sum = 0;
    double zmin = -_rout;
    double zmax = _rout;
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(_eps, _eps, zmin + k*(zmax-zmin)/NSAMPLES));
    }
    return (sum/NSAMPLES)*(zmax-zmin);
}

////////////////////////////////////////////////////////////////////
