/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEGRAINCOMPOSITION_HPP
#define FILEGRAINCOMPOSITION_HPP

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The FileGrainComposition class represents the optical and calorimetric properties of dust
    grains loaded from two data files provided by the user. The file names are specified as
    attributes of this class.

    The first file contains the optical grain properties. The file must have a simple text format
    as described here. Any initial lines that start with a # character are considered to be part of
    a header and are thus ignored. The first number on the first non-header line specifies the
    number of grain size grid points \f$N_a\f$; the first number on the second non-header line
    specifies the number of wavelength grid points \f$N_\lambda\f$. Subsequently there must be
    \f$N_a\f$ blocks containing \f$N_\lambda+1\f$ lines each. The first number on the first line in
    each block provides the grain size value for this block (in micron); blocks must be in order of
    increasing grain size. Subsequent lines in the block provide information for the various
    wavelengths at this grain size. These lines are in order of increasing wavelength. Each data
    line in a block must have four columns in the following order: \f$\lambda\f$,
    \f$Q^\text{abs}\f$, \f$Q^\text{sca}\f$, \f$g\f$. The wavelength must be given in micron; the
    other three values are dimensionless. For all lines discussed above, any additional information
    at the end of the line is ignored.

    The second file contains the calorimetric grain properties. The file should have a simple text
    format as described here. Any initial lines that start with a # character are considered to be
    part of a header and are thus ignored. The first number on the first non-header line specifies
    the number of temperature grid points \f$N_T\f$. The remaining lines tabulate the enthalpy for
    increasing temperature. Each line must have two columns specifying the temperature \f$T_t\f$
    (in K) and the corresponding enthalpy \f$h_t\f$ (in erg/g). Any additional information at the
    end of a line is ignored. */
class FileGrainComposition : public GrainComposition
{
    ITEM_CONCRETE(FileGrainComposition, GrainComposition, "a dust grain composition based on user-supplied data files")

    PROPERTY_STRING(opticalFilename, "the name of the file containing the optical grain properties")

    PROPERTY_STRING(calorimetricFilename, "the name of the file containing the calorimetric grain properties")

    PROPERTY_DOUBLE(bulkMassDensity, "the bulk mass density for this grain material")
        ATTRIBUTE_QUANTITY(bulkMassDensity, "bulkmassdensity")
        ATTRIBUTE_MIN_VALUE(bulkMassDensity, "[100 kg/m3")
        ATTRIBUTE_MAX_VALUE(bulkMassDensity, "10000 kg/m3]")
        ATTRIBUTE_DEFAULT_VALUE(bulkMassDensity, "2240 kg/m3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the optical and calorimetric properties from the specified input files,
        and sets the bulk mass density to the value of the corresponding attribute. */
    void setupSelfBefore() override;

    //====================== Identifying =====================

public:
    /** This function returns a brief human-readable identifier for the type of grain composition
        represented by the instance. For this class, the function returns a concatenation of
        the optical and calorimetric file names. */
    string name() const override;
};

////////////////////////////////////////////////////////////////////

#endif
