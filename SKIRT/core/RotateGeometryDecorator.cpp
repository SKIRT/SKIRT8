/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RotateGeometryDecorator.hpp"

//////////////////////////////////////////////////////////////////////

void RotateGeometryDecorator::setupSelfBefore()
{
    GenGeometry::setupSelfBefore();

    // cache frequently used values
    _sinalpha = sin(_eulerAlpha);
    _cosalpha = cos(_eulerAlpha);
    _sinbeta = sin(_eulerBeta);
    _cosbeta = cos(_eulerBeta);
    _singamma = sin(_eulerGamma);
    _cosgamma = cos(_eulerGamma);
    _R11 = _cosalpha*_cosgamma - _sinalpha*_cosbeta*_singamma;
    _R12 = _sinalpha*_cosgamma + _cosalpha*_cosbeta*_singamma;
    _R13 = _sinbeta*_singamma;
    _R21 = -_cosalpha*_singamma - _sinalpha*_cosbeta*_cosgamma;
    _R22 = -_sinalpha*_singamma + _cosalpha*_cosbeta*_cosgamma;
    _R23 = _sinbeta*_cosgamma;
    _R31 = _sinalpha*_sinbeta;
    _R32 = -_cosalpha*_sinbeta;
    _R33 = _cosbeta;
}

////////////////////////////////////////////////////////////////////

double RotateGeometryDecorator::density(Position bfr) const
{
    Position bfrorig = derotate(bfr);
    return _geometry->density(bfrorig);
}

////////////////////////////////////////////////////////////////////

Position RotateGeometryDecorator::generatePosition() const
{
    Position bfrorig = _geometry->generatePosition();
    Position bfr = rotate(bfrorig);
    return bfr;
}

////////////////////////////////////////////////////////////////////

double RotateGeometryDecorator::SigmaX() const
{
    return _geometry->SigmaX();
}

////////////////////////////////////////////////////////////////////

double RotateGeometryDecorator::SigmaY() const
{
    return _geometry->SigmaY();
}

////////////////////////////////////////////////////////////////////

double RotateGeometryDecorator::SigmaZ() const
{
    return _geometry->SigmaZ();
}

////////////////////////////////////////////////////////////////////

double RotateGeometryDecorator::probabilityForDirection(int ell, Position bfr, Direction bfk) const
{
    Position bfrorig = derotate(bfr);
    Direction bfkorig = derotate(bfk);
    return _geometry->probabilityForDirection(ell,bfrorig,bfkorig);
}

////////////////////////////////////////////////////////////////////

Direction RotateGeometryDecorator::generateDirection(int ell, Position bfr) const
{
    Position bfrorig = derotate(bfr);
    Direction bfkorig = _geometry->generateDirection(ell,bfrorig);
    Direction bfk = rotate(bfkorig);
    return bfk;
}

////////////////////////////////////////////////////////////////////

Position RotateGeometryDecorator::rotate(Position bfrorig) const
{
    double xorig = bfrorig.x();
    double yorig = bfrorig.y();
    double zorig = bfrorig.z();
    double x = _R11*xorig + _R12*yorig + _R13*zorig;
    double y = _R21*xorig + _R22*yorig + _R23*zorig;
    double z = _R31*xorig + _R32*yorig + _R33*zorig;
    return Position(x,y,z);
}

////////////////////////////////////////////////////////////////////

Position RotateGeometryDecorator::derotate(Position bfr) const
{
    double x = bfr.x();
    double y = bfr.y();
    double z = bfr.z();
    double xorig = _R11*x + _R21*y + _R31*z;
    double yorig = _R12*x + _R22*y + _R32*z;
    double zorig = _R13*x + _R23*y + _R33*z;
    return Position(xorig,yorig,zorig);
}

////////////////////////////////////////////////////////////////////

Direction RotateGeometryDecorator::rotate(Direction bfkorig) const
{
    double kxorig = bfkorig.x();
    double kyorig = bfkorig.y();
    double kzorig = bfkorig.z();
    double kx = _R11*kxorig + _R12*kyorig + _R13*kzorig;
    double ky = _R21*kxorig + _R22*kyorig + _R23*kzorig;
    double kz = _R31*kxorig + _R32*kyorig + _R33*kzorig;
    return Direction(kx,ky,kz);
}

////////////////////////////////////////////////////////////////////

Direction RotateGeometryDecorator::derotate(Direction bfk) const
{
    double kx = bfk.x();
    double ky = bfk.y();
    double kz = bfk.z();
    double kxorig = _R11*kx + _R21*ky + _R31*kz;
    double kyorig = _R12*kx + _R22*ky + _R32*kz;
    double kzorig = _R13*kx + _R23*ky + _R33*kz;
    return Direction(kxorig,kyorig,kzorig);
}

////////////////////////////////////////////////////////////////////
