/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SerializedData.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

// We keep the data vector empty as long as possible; this avoids allocations when SerializedData objects
// are created and then subseqeuently overwritten by a copy from another object.
//
// As soon as there is some data in the buffer, the first vector element specifies the number of elements actually
// in use, while the size of the data vector specifies the allocated number of elements.
//
// Pushing a data item causes the value or values of the item to be pushed to the data buffer, followed by
// the number of values pushed.

////////////////////////////////////////////////////////////////////

SerializedData::SerializedData()
{
}

////////////////////////////////////////////////////////////////////

void SerializedData::push(double d)
{
    canon();
    _data.push_back(d);
    _data.push_back(-1.);
    _data.front() += 2.;
}

////////////////////////////////////////////////////////////////////

void SerializedData::push(const vector<double>& v)
{
    canon();
    for (double d : v) _data.push_back(d);
    size_t n = v.size();
    _data.push_back(static_cast<double>(n));
    _data.front() += n+1.;
}

////////////////////////////////////////////////////////////////////

size_t SerializedData::used() const
{
    return _data.empty() ? 0 : static_cast<size_t>(_data.front());
}

////////////////////////////////////////////////////////////////////

const double* SerializedData::data() const
{
    if (_data.empty()) throw FATALERROR("Can't provide pointer to empty serialized data buffer");
    return &_data[0];
}

////////////////////////////////////////////////////////////////////

SerializedData::SerializedData(size_t allocate)
    : _data( max(allocate, static_cast<size_t>(1)) )
{
    _data.front() = 1.;
}

////////////////////////////////////////////////////////////////////

size_t SerializedData::allocated() const
{
    return _data.size();
}

////////////////////////////////////////////////////////////////////

double* SerializedData::data()
{
    if (_data.empty()) _data.push_back(1.);
    return &_data[0];
}

////////////////////////////////////////////////////////////////////

void SerializedData::shrink()
{
    if (used() < allocated()) _data.resize(used());
}

////////////////////////////////////////////////////////////////////

double SerializedData::pop()
{
    canon();
    if (_data.front() <= 1.) throw FATALERROR("Attempt to pop data item from empty serialized data buffer");
    if (_data.back() != -1.) throw FATALERROR("Attempt to pop a scalar but serialized data buffer has a vector");
    _data.pop_back();
    double d = _data.back();
    _data.pop_back();
    _data.front() -= 2.;
    return d;
}

////////////////////////////////////////////////////////////////////

void SerializedData::pop(double& d)
{
    canon();
    if (_data.front() <= 1.) throw FATALERROR("Attempt to pop data item from empty serialized data buffer");
    if (_data.back() != -1.) throw FATALERROR("Attempt to pop a scalar but serialized data buffer has a vector");
    _data.pop_back();
    d = _data.back();
    _data.pop_back();
    _data.front() -= 2.;
}

////////////////////////////////////////////////////////////////////

void SerializedData::pop(vector<double>& v)
{
    canon();
    if (_data.front() <= 1.) throw FATALERROR("Attempt to pop data item from empty serialized data buffer");
    if (_data.back() < 0.) throw FATALERROR("Attempt to pop a vector but serialized data buffer has a scalar");
    size_t n = static_cast<size_t>(_data.back());
    _data.pop_back();

    v.assign(_data.end()-n, _data.end());
    _data.erase(_data.end()-n, _data.end());
    _data.front() -= n+1.;
}

////////////////////////////////////////////////////////////////////

void SerializedData::canon()
{
    if (_data.empty()) _data.push_back(1.);
    else shrink();
}

////////////////////////////////////////////////////////////////////
