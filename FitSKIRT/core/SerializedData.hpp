/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SERIALIZEDDATA_HPP
#define SERIALIZEDDATA_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** This class allows to serialize a sequence of data items of type double or vector<double> into a
    linear buffer of double values that can be transmitted between processes, and to subsequently
    deserialize the transmitted buffer back to a copy of the original sequence of data items. Data
    items must be pushed onto a SerializedData object and then popped from the same object or a
    (transmitted) copy in reverse order for retrieval. The sender and receiver must agree on the
    type, number and precise order of the transmitted data items.

    Because the total number of double values in the buffer is stored as one of double values in
    the buffer, the size of the buffer cannot exceed \f$2^{53} \approx 9\times10^{15}\f$ double
    values. In practice, this is not much of a limitation. */
class SerializedData final
{
    // ============= Sending ================

public:
    /** The default constructor creates a SerializedData object with an empty buffer, which will
        grow as needed while data items are being pushed. */
    SerializedData();

    /** This function pushes the specified double value onto the serialized data. */
    void push(double d);

    /** This function pushes the specified vector of double values onto the serialized data. */
    void push(const vector<double>& v);

    /** This function returns the number of double values currently used in the serialized data. */
    size_t used() const;

    /** This function returns a read-only pointer to the raw serialized data. The pointer is
        invalidated by invoking any of the non-const functions on this object except data(). */
    const double* data() const;

    // ============= Receiving ================

public:
    /** This constructor creates a SerializedData object with a buffer including the specified
        number of double values, initialized such that no values are actually in use. This
        constructor is typically used to prepare a SerializedData object for receiving transmitted
        data. The buffer allocation is invalidated by invoking any of the non-const functions on
        this object except data(). */
    SerializedData(size_t allocate);

    /** This function returns the number of double values currently allocated for the serialized
        data. */
    size_t allocated() const;

    /** This function returns a writable pointer to the raw serialized data of the object. The
        pointer is invalidated by invoking any of the non-const functions on this object except
        data(). */
    double* data();

    /** This function shrinks the buffer allocation to match the number of values actually in use.
        The function is usually called after the buffer was externally overwritten with data
        received from a remote SerializedData object. */
    void shrink();

    /** This function pops the most recent data item from the serialized data, and returns its
        value as a double. */
    double pop();

    /** This function pops the most recent data item from the serialized data and stores it in the
        specified double argument, overwriting any previous contents. */
    void pop(double& d);

    /** This function pops the most recent data item from the serialized data and stores it in the
        specified vector of doubles argument, overwriting any previous contents. */
    void pop(vector<double>& v);

    // ============= Private functions ================

private:
    /** This function transforms the data buffer to canonical form: the size specified in the first
        item matches the size of the vector. This implies that the vector is nonempty and that the
        functions used() and allocated() return the same value. */
    void canon();

    // ============= Data Members ================

private:
    // the serialized data buffer
    vector<double> _data;
};

////////////////////////////////////////////////////////////////////

#endif
