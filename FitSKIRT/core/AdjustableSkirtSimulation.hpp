/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADJUSTABLESKIRTSIMULATION_HPP
#define ADJUSTABLESKIRTSIMULATION_HPP

#include "SimulationItem.hpp"
#include <unordered_map>

////////////////////////////////////////////////////////////////////

/** The AdjustableSkirtSimulation class allows performing a SKIRT simulation loaded from a ski
    file. The contents of the ski file can be adjusted before the simulation hierarchy is actually
    created, as described for the performWith() function. */
class AdjustableSkirtSimulation : public SimulationItem
{
    ITEM_CONCRETE(AdjustableSkirtSimulation, SimulationItem, "an adjustable SKIRT simulation")

    PROPERTY_STRING(skiFilename, "the name of the ski file specifying the SKIRT simulation")

    ITEM_END()

    //======== Construction - Setup - Run - Destruction  ===========

protected:
    /** This function reads the specified ski file into memory, constructs a single simulation
        using the default values provided in the ski file, and extracts and remembers relevant
        information from the default simulation hierarchy. It also verifies that the default
        simulation hierarchy includes a multi-frame instrument and that the fit scheme in which
        this object resides has the same type of unit system, so that the implementation of this
        class can use the fit scheme unit system to replace values in the default simulation
        hierarchy. */
    void setupSelfBefore() override;

    //====================== Other functions =======================

public:
    /** A shorthand type definition for a replacement dictionary as passed to performWith(). Refer
        to the description of the performWith() function for more information. */
    typedef std::unordered_map<string, std::pair<double,string>> ReplacementDict;

    /** This function returns the number of luminosity components in the simulation. */
    size_t numComponents() const;

    /** This function returns the number of wavelengths in the simulation, which is the same as the
        number of instrument frames in the simulation. */
    size_t numWavelengths() const;

    /** This function returns the name of the instrument used in the simulation. */
    string instrumentName() const;

    /** This function returns the number of horizontal pixels in the frame at index \em ell. */
    int pixelsX(int ell) const;

    /** This function returns the number of vertical pixels in the frame at index \em ell. */
    int pixelsY(int ell) const;

    /** This function returns the x increment of the frame at index \em ell. */
    double pixelSizeX(int ell) const;

    /** This function returns the y increment of the frame at index \em ell. */
    double pixelSizeY(int ell) const;

    /** This function returns the luminosity specified in the simulation for the stellar component
        at index \em k at wavelength index \em ell. */
    double luminosity(int k, int ell) const;

    /** This function runs the SKIRT simulation specified by the previously loaded ski file after
        adjusting its contents as specified through the \em replacements dictionary argument (see
        detailed description below). The prefix string, if present, is appended to the filename
        prefix for all output files of this simulation run.

        The \em replacements dictionary contains a set of key/value pairs controlling replacement
        of labeled attribute values in the ski file. The key for each dictionary item is a string
        matching an attribute label in the ski file, as explained below. The value for each
        dictionary item is in fact a pair of values: the numeric replacement value in SI units, and
        a physical quantity specifier such as "length" (or the empty string if the value is
        dimensionless).

        To mark a numeric attribute value for replacement in the ski file, enclose the value in
        square brackets and provide a label, as in the following example:
        \verbatim
        <SersicGeometry index="3" radius="[stellar_scale:1500 pc]"/>
        \endverbatim
        The brackets must be just within the quotes delimiting the attribute value. The label must
        start with a letter and contain only letters, digits and underscores. It must be
        immediately followed by a colon and then the regular attribute value, possibly including a
        unit specifier. Spaces are not allowed except between the value and the unit specifier,
        where a single space is required.

        If the label matches one of the keys in the replacement dictionary handed to this function,
        the corresponding value is substituted in the ski file, removing the brackets and the
        label. If the label does not match one of the keys in the replacement dictionary, a fatal
        error is thrown. */
    void performWith(const ReplacementDict& replacements, string prefix=string());

private:
    /** This private function performs the specified adjustments on the previously loaded ski
        content as described for the performWith() function, and returns the result. */
    string adjustedSkiContent(const ReplacementDict& replacements);

    //======================== Data Members ========================

private:
    // initialized during setup
    string _skiContent;                 // the content of the ski file, without modifications

    // information extracted from the default simulation hierarchy during setup
    size_t _ncomponents{0};             // the number of stellar components
    size_t _nwavelengths{0};            // the number of wavelentghs/instrument frames
    string _instrname;                  // the name of the multi-frame instrument
    vector<size_t> _pixelsX;            // the number of x pixels for each frame
    vector<size_t> _pixelsY;            // the number of y pixels for each frame
    vector<double> _pixelSizeX;         // the x increment for each frame
    vector<double> _pixelSizeY;         // the y increment for each frame
    vector<vector<double>> _luminosities; // the luminosities per wavelength (inner) for each stellar component (outer)
};

////////////////////////////////////////////////////////////////////

#endif
