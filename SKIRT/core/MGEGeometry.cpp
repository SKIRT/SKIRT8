/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MGEGeometry.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "System.hpp"

//////////////////////////////////////////////////////////////////////

void MGEGeometry::setupSelfBefore()
{
    AxGeometry::setupSelfBefore();

    // read in the file with the raw MGE data

    string filepath = find<FilePaths>()->input(_filename);
    std::ifstream file = System::ifstream(filepath);
    if (! file.is_open()) throw FATALERROR("Could not open the MGE expansion data file " + filepath);
    find<Log>()->info("Reading MGE expansion data from file " + filepath + "...");
    double Mprev = -9e99;
    double sigmaprev = -9e99;
    double qprev = -9e99;
    double M, sigma, q;
    while (!file.eof())
    {
        file >> M >> sigma >> q;
        if ( M != Mprev || sigma != sigmaprev || q != qprev)
        {
            _Mv.push_back(M);
            Mprev = M;
            _sigmav.push_back(sigma);
            sigmaprev = sigma;
            _qv.push_back(q);
            qprev = q;
        }
    }
    file.close();
    find<Log>()->info("File " + filepath + " closed.");
    _Ncomp = _Mv.size();

    // convert from pixelscale to physical scale

    for (int i=0; i<_Ncomp; i++)
        _sigmav[i] *= _pixelscale;

    // convert the apparent flattening to real flattening
    // (see e.g. Bacon 1985, A&A, 143, 84)

    double cosi = cos(_inclination);
    double sini = sin(_inclination);
    for (int i=0; i<_Ncomp; i++)
    {
        if (_qv[i] < cosi) throw FATALERROR("MGE component with index " + std::to_string(i) + " can't be deprojected: "
                                "apparent flattening is smaller than cosine of inclination ("
                                + StringUtils::toString(_qv[i], 'f') + " < " + StringUtils::toString(cosi, 'f') + ")");
        _qv[i] = sqrt((_qv[i]-cosi)*(_qv[i]+cosi))/sini;
    }

    // convert the counts to normalized luminosity and set up a vector
    // with cumulative luminosities

    double Mtot = 0.0;
    _Mcumv.resize(_Ncomp);
    for (int i=0; i<_Ncomp; i++)
    {
        Mtot += _Mv[i];
        _Mcumv[i] = Mtot;
    }
    for (int i=0; i<_Ncomp; i++)
    {
        _Mv[i] /= Mtot;
        _Mcumv[i] /= Mtot;
    }
}

////////////////////////////////////////////////////////////////////

double MGEGeometry::density(double R, double z) const
{
    double rho = 0.0;
    for (int i=0; i<_Ncomp; i++)
    {
        double q = _qv[i];
        double M = _Mv[i];
        double sigma = _sigmav[i];
        double rho0 = M / pow(sqrt(2.0*M_PI)*sigma,3) / q;
        double m2 = R*R + z*z/(q*q);
        double sigma2 = sigma*sigma;
        rho += rho0 * exp(-0.5*m2/sigma2);
    }
    return rho;
}

////////////////////////////////////////////////////////////////////

Position MGEGeometry::generatePosition() const
{
    double X = random()->uniform();
    for (int i=0; i<_Ncomp; i++)
    {
        if (X<=_Mcumv[i])
        {
            double q = _qv[i];
            double sigma = _sigmav[i];
            double x = sigma * random()->gauss();
            double y = sigma * random()->gauss();
            double z = q * sigma * random()->gauss();
            return Position(x,y,z);
        }
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

double MGEGeometry::SigmaR() const
{
    double sum = 0.0;
    for (int i=0; i<_Ncomp; i++)
    {
        double sigma = _sigmav[i];
        double sigma2 = sigma*sigma;
        sum += _Mv[i]/(4.0*M_PI)/sigma2/_qv[i];
    }
    return sum;
}

//////////////////////////////////////////////////////////////////////

double MGEGeometry::SigmaZ() const
{
    double sum = 0.0;
    for (int i=0; i<_Ncomp; i++)
    {
        double sigma = _sigmav[i];
        double sigma2 = sigma*sigma;
        sum += _Mv[i]/(2.0*M_PI)/sigma2;
    }
    return sum;
}

////////////////////////////////////////////////////////////////////
