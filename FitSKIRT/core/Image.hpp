/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMAGE_HPP
#define IMAGE_HPP

#include "Array.hpp"
class SimulationItem;

////////////////////////////////////////////////////////////////////

/** An instance of this class represents a single 2D image frame where each pixel is stored as a
    (double precision) floating point value. The class offers facilities to load and save an image
    from and to a FITS file, and to access and perform calculations with the pixel data. */
class Image
{
    //=============== Constructing and Loading  ==================

public:
    /** The default constructor creates an empty image. */
    Image();

    /** This constructor creates an image from the specified FITS file. The arguments of this
        constructor are the same as those for the load() function; see the description of that
        function for more information. */
    Image(const SimulationItem* item, string filename, string directory = string(), bool log = true);

    /** This function loads an image from a FITS file in the context of the simulation item
        hierarchy specified through the first argument, replacing any previous information. This
        allows the function to issue log messages and to determine, for example, the simulation's
        input path. The second argument specifies the filename of the FITS file to be loaded. If
        needed, the function adds the ".fits" filename extension to the filename. If the third
        argument is a non-empty string, it is taken as the path to the directory in which the file
        resides. Otherwise, the filename is interpreted relative to the simulation's input path. If
        the last argument is provided and set to false, informational log messages are suppressed.
        */
    void load(const SimulationItem* item, string filename, string directory = string(), bool log = true);

    /** This function resizes the image to a specified width and height, erasing any previous pixel
        data. */
    void resize(int nx, int ny);

    /** This function copies the data from the specified array into the image's own pixel data
        array. The size of the specified array must correspond to the current dimensions of the
        receiving image. */
    void setData(const Array& data);

    /** This function moves the data from the specified array into the image's own pixel data array
        (invalidating the contents of the incoming array in the process). The size of the specified
        array must correspond to the current dimensions of the receiving image. */
    void moveData(Array&& data);

    //========================= Saving ==========================

public:
    /** This function saves the image to a FITS file in the context of the simulation item
        hierarchy specified through the first argument. This allows the function to issue log
        messages using the human-readable description of the data given as the second argument, to
        skip writing in non-root processes of a multi-process similation, and to determine a full
        path from the output filename specified as the third argument, relative to the simulation's
        output path. The function also adds the simulation output prefix and, if needed, the
        ".fits" filename extension to the filename. The remaining arguments specify the physical
        properties of the image, which will be included in the FITS header information. */
    void save(const SimulationItem* item, string description, string filename,
              double incx, double incy, string dataUnits, string xyUnits);

    //===================== Accessing properties =======================

public:
    /** This function returns the number of pixels in the x direction. */
    int sizeX() const
    {
        return _nx;
    }

    /** This function returns the number of pixels in the y direction. */
    int sizeY() const
    {
        return _ny;
    }

    /** This function returns the total number of pixels in this image. */
    int size() const
    {
        return _nx*_ny;
    }

    //===================== Accessing pixel data =======================

public:
    /** This function returns a constant reference to the pixel data array. */
    const Array& data() const
    {
        return _data;
    }

    /** This operator returns a writable reference to the pixel value at a given index. */
    double& operator[](int i)
    {
        return _data[i];
    }

    /** This operator returns the pixel value at a given index. */
    double operator[](int i) const
    {
        return _data[i];
    }

    /** This operator returns a writable reference to the pixel value defined by the given x and y
        coordinates in the image. */
    double& operator()(int x, int y)
    {
        return _data[x + _nx*y];
    }

    /** This operator returns the pixel value defined by the given x and y coordinates in the
        image. */
    double operator()(int x, int y) const
    {
        return _data[x + _nx*y];
    }

    //===================== Calculating with pixel data =======================

public:
    /** This operator multiplies each pixel in the image by the given value. */
    Image& operator*= (double value)
    {
        _data *= value;
        return *this;
    }

    /** This operator divides each pixel in the image by the given value. */
    Image& operator/= (double value)
    {
        _data /= value;
        return *this;
    }

    /** This operator adds another image to this image, pixel by pixel. The two images must have
        the same size. */
    Image& operator+= (const Image& other)
    {
        _data += other.data();
        return *this;
    }

    /** This operator subtracts another image from this image, pixel by pixel. The two images must
        have the same size. */
    Image& operator-= (const Image& other)
    {
        _data -= other.data();
        return *this;
    }

    /** This operator divides this image by another image, pixel by pixel. The two images must have
        the same size. */
    Image& operator/= (const Image& other)
    {
        _data /= other.data();
        return *this;
    }

    /** This function returns the sum of all the pixel values in the image. */
    double sum() const
    {
        return _data.sum();
    }

    /** This function applies the absolute value operator to each pixel in this image. */
    Image& abs()
    {
        for (size_t i = 0; i<_data.size(); ++i) _data[i] = std::abs(_data[i]);
        return *this;
    }

    //======================== Data Members ========================

private:
    // the pixel data array
    Array _data;

    // the pixel dimensions of the image
    int _nx{0};
    int _ny{0};
};

//===================== Calculating with pixel data =======================

/** This operator returns a new image created by multiplying each pixel in the given image by the
    given value. */
inline Image operator* (const Image& lhs, double rhs)
{
    Image result = lhs;
    return result *= rhs;
}

/** This operator returns a new image created by adding the given images pixel by pixel. The two
    images must have the same size. */
inline Image operator+ (const Image& lhs, const Image& rhs)
{
    Image result = lhs;
    return result += rhs;
}

/** This operator returns a new image created by subtracting the given images pixel by pixel. The
    two images must have the same size. */
inline Image operator- (const Image& lhs, const Image& rhs)
{
    Image result = lhs;
    return result -= rhs;
}

/** This operator returns a new image created by dividing the first image by the second pixel by
    pixel. The two images must have the same size. */
inline Image operator/ (const Image& lhs, const Image& rhs)
{
    Image result = lhs;
    return result /= rhs;
}

////////////////////////////////////////////////////////////////////

#endif
