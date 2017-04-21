/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTSYSTEM_HPP
#define DUSTSYSTEM_HPP

#include "DustDistribution.hpp"
#include "DustGrid.hpp"
#include "Array.hpp"
#include "Position.hpp"
#include "Table.hpp"
#include <mutex>
class DustGridDensityInterface;
class DustMix;
class PhotonPackage;
class ProcessAssigner;

//////////////////////////////////////////////////////////////////////

/** DustSystem is an abstract class for representing complete dust systems, including both a dust
    distribution (i.e. a complete description of the spatial distribution and optical properties of
    the dust) and a grid on which this distribution is discretized. There are specialized
    subclasses for use with oligochromatic and panchromatic simulations respectively. A DustSystem
    object contains a vector of dust cells, each of which contain all useful information on the
    dust within that particular piece of the configuration space. Furthermore, a DustSystem object
    contains pointers a DustDistribution object and to a DustGrid object). A subclass may of course
    maintain additional information depending on its needs. */
class DustSystem : public SimulationItem
{
    ITEM_ABSTRACT(DustSystem, SimulationItem, "a dust system")

    PROPERTY_ITEM(dustDistribution, DustDistribution, "the dust distribution")
        ATTRIBUTE_DEFAULT_VALUE(dustDistribution, "CompDustDistribution")

    PROPERTY_ITEM(dustGrid, DustGrid, "the dust grid")
        ATTRIBUTE_DEFAULT_VALUE(dustGrid, "OctTreeDustGrid")

    PROPERTY_INT(numSamples, "the number of random density samples for determining cell mass")
        ATTRIBUTE_MIN_VALUE(numSamples, "10")
        ATTRIBUTE_MAX_VALUE(numSamples, "1000")
        ATTRIBUTE_DEFAULT_VALUE(numSamples, "100")

    PROPERTY_BOOL(writeConvergence, "output a data file with convergence checks on the dust system")
        ATTRIBUTE_DEFAULT_VALUE(writeConvergence, "true")

    PROPERTY_BOOL(writeDensity, "output FITS files displaying the dust density distribution")
        ATTRIBUTE_DEFAULT_VALUE(writeDensity, "true")

    PROPERTY_BOOL(writeDepthMap, "output FITS file with a V-band optical depth map seen from the center")
        ATTRIBUTE_DEFAULT_VALUE(writeDepthMap, "false")

    PROPERTY_BOOL(writeQuality, "calculate and output quality metrics for the dust grid")
        ATTRIBUTE_DEFAULT_VALUE(writeQuality, "false")
        ATTRIBUTE_SILENT(writeQuality)

    PROPERTY_BOOL(writeCellProperties, "output a data file with relevant properties for all dust cells")
        ATTRIBUTE_DEFAULT_VALUE(writeCellProperties, "false")
        ATTRIBUTE_SILENT(writeCellProperties)

    PROPERTY_BOOL(writeCellsCrossed, "output statistics on the number of cells crossed per path")
        ATTRIBUTE_DEFAULT_VALUE(writeCellsCrossed, "false")
        ATTRIBUTE_SILENT(writeCellsCrossed)

    PROPERTY_BOOL(writeStellarDensity, "output normalized stellar density on dust grid")
        ATTRIBUTE_DEFAULT_VALUE(writeStellarDensity, "false")
        ATTRIBUTE_SILENT(writeStellarDensity)

    ITEM_END()

    /** \fn writeStellarDensity
        The writeStellarDensity flag is intended for rare specialty
        situations where one wants to get information on the discretization of the stellar
        distribution over the dust grid. This makes sense and is supported only when the stellar
        system consists of a single component of type GeometricStellarComp. In other cases, one
        would have to take into account the lumonosity normalization of the respective components
        and/or particles or cells of the input distribution. If the flag is enabled with an
        unsupported stellar system, a fatal error occurs during setup. */

    //============= Construction - Setup - Destruction =============

protected:
    /** This function performs setup for the dust system, which includes several tasks. First, the
        function verifies that either all dust mixes in the dust system support polarization, or
        none of them do. The next task consists of calculating and storing the volume and the dust
        density (for every dust component) of all the cells. To calculate the volume of a given
        dust cell, we just call the corresponding function of the dust grid. To calculate and set
        the density corresponding to given dust component in a given cell, a number of random
        positions are generated within the cell (see sampleCount()). The density in the cell is
        calculated as the mean of the density values (found using a call to the corresponding
        function of the dust distribution) in these points. The calculation of both volume and
        density is parallellized. Finally, the function optionally invokes various writeXXX()
        functions depending on the state of the corresponding write flags. */
    void setupSelfAfter() override;

private:
    /** This function serves as the parallelization body for calculating the volume of each cell.
        */
    void setVolumeBody(size_t m);

    /** This function serves as the parallelization body for setting the density value of each cell
        through the DustGridDensityInterface interface, if available. */
    void setGridDensityBody(size_t m);

    /** This function serves as the parallelization body for setting the density value of each cell
        by taking random density sample. */
    void setSampleDensityBody(size_t m);

    /** This function is used to assemble the container that stores the densities of all dust cells
        for each dust component. If multiprocessing is enabled, the calculation of these densities
        can be performed in parallel by the different processes, depending on the type of
        ProcessAssigner that is used for this dust system. When a ProcessAssigner subclass is used
        which distributes the calculation of the dust cell densities amongst the different parallel
        processes, each process contains only the densities for a particular set of dust cells (but
        for each dust component). Therefore, this function uses the PeerToPeerCommunicator object
        to broadcast the densities of the dust cells assigned to a particular process to all other
        processes and storing them in the appropriate place in the container. This is implemented
        as an element-wise summation of this container across all processes, where the density for
        a particular dust cell and dust component will be added to a series of zeros, coming from
        the processes that were not assigned to that dust cell. The end result will be an assembly
        of the densities, where each process stores the density over the entire dust grid. */
    void assemble();

    /** This function writes out a simple text file, named <tt>prefix_ds_convergence.dat</tt>,
        providing a convergence check on the dust system. The function calculates the total dust
        mass, the face-on surface density and the edge-on surface density by directly integrating
        over the dust grid, and compares these values with the expected "theoretical" values found
        by just calling the corresponding functions of the dust distribution. The results are
        written to the file and can be studied to see whether the chosen dust grid is adequate for
        the chosen dust distribution. */
    void doWriteConvergence() const;

    /** This function writes out FITS files with the theoretical dust density and the
        grid-discretized dust density in the coordinate planes. Each of these maps contains 1024 x
        1024 pixels, and covers as a field of view the total extension of the grid. The number of
        data files written depends on the dimension of the dust system's geometry: for spherical
        symmetry only the intersection with the xy plane is written, for axial symmetry the
        intersections with the xy and xz planes are written, and for general geometries all three
        intersections are written. The difference between the theoretical dust density maps (named
        <tt>prefix_ds_trhoXX.fits</tt>) and the grid-discretized dust density maps (named
        <tt>prefix_ds_grhoXX.fits</tt>) is the following: the theoretical dust density is the total
        dust density of the dust distribution, i.e.\ the actual dust density that would correspond
        to an infinitely fine dust grid. The grid-discretized dust density maps on the other hand
        give the value of the dust as read from the finite-resolution dust grid. A comparison of
        both sets of maps can reveal whether the chosen dust grid is suitable (in the ideal case,
        there would be no difference between both sets of maps). */
    void doWriteDensity() const;

    /** This function writes out a FITS file named <tt>prefix_ds_tau.fits</tt> with an all-sky
        V-band optical depth map as seen from the coordinate origin. The map has 1600 x 800 pixels
        and uses the Mollweide projection to project the complete sky onto a proportional 2:1
        ellipse. The values outside of the ellipse are set to zero. The direction
        \f$\bf{k}=(\theta,\phi)\f$ corresponding to the pixel in the map with horizontal and
        vertical indices \f$(i,j)\f$ can be found through the inverse Mollweide projection, which
        in this case can be written as follows: \f[ x=(i+\frac{1}{2})/N_\mathrm{pixels,x} \f] \f[
        y=(j+\frac{1}{2})/N_\mathrm{pixels,y} \f] \f[ \alpha=\arcsin(2y-1) \f] \f[
        \theta=\arccos(\frac{2\alpha+\sin 2\alpha}{\pi}) \f] \f[ \phi=\frac{\pi(2x-1)}{\cos\alpha}
        \f] */
    void doWriteDepthMap() const;

    /** This function writes out a simple text file, named <tt>prefix_ds_quality.dat</tt>,
        providing some basic quality metrics for the dust grid. The first metric consists of the
        mean value and the standard deviation for the absolute difference \f$|\rho_g-\rho_t|\f$
        between the theoretical and grid density in a large number of randomly chosen points,
        uniformly distributed over the dust grid volume. The second metric consists of the mean
        value and the standard deviation for the difference \f$|\tau_g-\tau_t|\f$ between the
        theoretical and grid optical depth, calculated for a large number of line segments with
        random end points uniformly distributed over the dust grid volume. */
    void doWriteQuality() const;

    /** This function writes out a text data file, named <tt>prefix_ds_cellprops.dat</tt>, which
        contains a line for each cell in the dust grid. Each line contains columns representing the
        following dust cell properties: cell index, x,y,z coordinates of the cell center, volume,
        density, mass, mass fraction relative to total dust mass in the distribution, and V-band
        optical depth of the cell diagonal. The values are listed in the corresponding output
        units, except for the last two which are dimensionless quantities.

        In addition, this function calculates some statistics on the optical depth of the cells,
        and outputs them both to the text data file and to the message log. */
    void doWriteCellProperties() const;

    /** This function outputs a text data file, named <tt>prefix_ds_stellar.dat</tt>, which
        contains a line for each cell in the dust grid. Each line contains the dust cell index and
        the normalized density of the stellar distribution sampled at the cell center. This makes
        sense and is supported only when the stellar system consists of a single component of type
        GeometricStellarComp. In other cases, one would have to take into account the lumonosity
        normalization of the respective components and/or particles or cells of the input
        distribution. If invoked for a simulation configuration with an unsupported stellar system,
        a fatal error occurs during setup. */
    void doWriteStellarDensity() const;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the dust system, which depends on the (lack of)
        symmetry in the geometry of its distribution. A value of 1 means spherical symmetry, 2
        means axial symmetry and 3 means none of these symmetries. */
    int dimension() const;

    /** This function returns the number of dust cells. */
    int numCells() const;

    /** This function returns the number of dust components. */
    int numComponents() const;

    /** This function returns true if the dust mixes in this dust system support polarization;
        false otherwise. During setup it is verified that either all dust mixes support
        polarization, or none of them do. */
    bool polarization() const;

    /** This function returns a pointer to the dust mixture corresponding to the \f$h\f$'th dust
        component. */
    DustMix* mix(int h) const;

    /** This function returns a pointer to a dust mixture that is selected randomly among the dust
        mixes of the dust components in the dust system. If we have just a single dust component,
        this is simple. If there are multiple dust components, the relative probability of
        selecting the dust mix for a particular dust component \f$h\f$ is given by
        \f$\kappa_\ell^{\text{sca}}(h)\,\rho({\bf{r}},h)\f$ at wavelength index \f$\ell\f$. */
    DustMix* randomMixForPosition(Position bfr, int ell) const;

    /** This function returns the number of the dust cell that contains the position
        \f${\boldsymbol{r}}\f$. The function just passes the call to corresponding function of the
        dust grid. */
    int whichCell(Position bfr) const;

    /** This function returns a random location in the dust cell with cell number \f$m\f$. The
        function just passes the call to corresponding function of the dust grid. */
    Position randomPositionInCell(int m) const;

    /** This function returns the volume of the dust cell with cell number \f$m\f$. */
    double volume(int m) const;

    /** This function returns the dust density corresponding to dust component \f$h\f$ of the dust
        cell with cell number \f$m\f$. If \f$m=-1\f$, i.e. if the cell number corresponds to a
        non-existing cell outside the grid, the value zero is returned. */
    double density(int m, int h) const;

    /** This function returns the total dust density of the dust cell with cell number \f$m\f$. If
        \f$m=-1\f$, i.e. if the cell number corresponds to a non-existing cell outside the grid,
        the value zero is returned. */
    double density(int m) const;

    /** This function returns a vector with the mean radiation field \f$J_{\ell,m}\f$ at all
        wavelength indices in the dust cell with cell number \f$m\f$. It is calculated as \f[
        J_{\ell,m} = \frac{ L_{\ell,m}^{\text{abs}} }{ 4\pi\, V_m\, (\Delta\lambda)_\ell \sum_h
        \kappa_{\ell,h}^{\text{abs}}\, \rho_{m,h} } \f] with \f$L_{\ell,m}^{\text{abs}}\f$ the
        absorbed luminosity, \f$\kappa_{\ell,h}^{\text{abs}}\f$ the absorption coefficient
        corresponding to the \f$h\f$'th dust component, \f$\rho_{m,h}\f$ the dust density
        corresponding to the \f$h\f$'th dust component, and \f$V_m\f$ the volume of the cell. */
    Array meanIntensity(int m) const;

    /** This function calculates the optical depth
        \f$\tau_{\ell,{\text{path}}}({\boldsymbol{r}},{\boldsymbol{k}})\f$ at wavelength index
        \f$\ell\f$ along a path through the dust system starting at the position
        \f${\boldsymbol{r}}\f$ into the direction \f${\boldsymbol{k}}\f$, where \f$\ell\f$,
        \f${\boldsymbol{r}}\f$ and \f${\boldsymbol{k}}\f$ are obtained from the specified
        PhotonPackage object, and it stores the resulting details back into the photon package
        object. The hard work is done by calling the DustGrid::path() function which stores the
        geometrical information on the path through the dust grid into the photon package: the cell
        numbers \f$m\f$ of the cells that are crossed by the path, the pathlength \f$(\Delta
        s)_m\f$ covered in that particular cell and a total path length counter \f$s_m\f$ that
        gives the total path length covered between the starting point \f${\boldsymbol{r}}\f$ and
        the boundary of the cell. With this information given, the calculation of the optical depth
        is rather straightforward: it is calculated as \f[
        \tau_{\ell,{\text{path}}}({\boldsymbol{r}},{\boldsymbol{k}}) = \sum_m (\Delta s)_m \sum_h
        \kappa_{\ell,h}^{\text{ext}}\, \rho_m, \f] where \f$\kappa_{\ell,h}^{\text{abs}}\f$ is the
        extinction coefficient corresponding to the \f$h\f$'th dust component at wavelength index
        \f$\ell\f$ and \f$\rho_{m,h}\f$ the dust density in the cell with cell number \f$m\f$
        corresponding to the \f$h\f$'th dust component. The function also stores the details on the
        calculation of the optical depth in the photon package, specifically it stores the optical
        depth covered within the \f$m\f$'th dust cell, \f[ (\Delta\tau_\ell)_m = (\Delta s)_m
        \sum_h \kappa_{\ell,h}^{\text{ext}}\, \rho_m, \f] and the total optical depth
        \f$\tau_{\ell,m}\f$ covered between the starting point \f${\boldsymbol{r}}\f$ and the
        boundary of the cell. */
    void fillOpticalDepth(PhotonPackage* pp);

    /** This function returns the optical depth
        \f$\tau_{\ell,{\text{d}}}({\boldsymbol{r}},{\boldsymbol{k}})\f$ at wavelength index
        \f$\ell\f$ along a path through the dust system starting at the position
        \f${\boldsymbol{r}}\f$ into the direction \f${\boldsymbol{k}}\f$ for a distance \f$d\f$,
        where \f$\ell\f$, \f${\boldsymbol{r}}\f$ and \f${\boldsymbol{k}}\f$ are obtained from the
        specified PhotonPackage object. The function first determines the photon package's path
        through the dust grid, storing the geometric information about the path segments through
        each cell into the photon package, and then calculates the optical depth at the specified
        distance. The calculation proceeds as described for the fillOpticalDepth() function; the
        differences being that the path length is limited to the specified distance, and that this
        function does not store the optical depth information back into the PhotonPackage object.
        */
    double opticalDepth(PhotonPackage* pp, double distance);

    /** If the writeCellsCrossed attribute is true, this function writes out a data file (named
        <tt>prefix_ds_crossed.dat</tt>) with statistics on the number of dust grid cells crossed
        per path calculated through the grid. The first column on each line specifies a particular
        number of cells crossed; the second column indicates the number of paths that crossed this
        precise number of cells. In effect this provides a histogram for the distribution of the
        path length (measured in the number of cells crossed). This virtual function can be
        overridden in a subclass to write out additional results of the simulation stored in the
        dust system. In that case, the overriding function must also call the implementation in
        this base class. */
    virtual void write() const;

    /** This pure virtual function must be implemented in each subclass to indicate whether dust
        emission is turned on for this dust system. The function returns true if dust emission is
        turned on, and false otherwise. It is provided in this base class because it is invoked
        from the general MonteCarloSimulation class. */
    virtual bool hasDustEmission() const = 0;

    /** This pure virtual function must be implemented in each subclass to indicate whether the
        absorption rates in each cell need to be stored for this dust system. This is needed if
        dust emission is turned on, but it can also be chosen if the user wants to study the mean
        intensity of the radiation field. It is provided in this base class because it is invoked
        from the general MonteCarloSimulation class. */
    virtual bool hasDustAbsorption() const = 0;

    /** This pure virtual function must be implemented in each subclass to simulate absorption of
        of a monochromatic luminosity package in the specified dust cell. The function should be
        called only if dustemission() returns true. It is provided in this base class because it is
        referenced from the general MonteCarloSimulation class (although it is actually invoked
        only for panchromatic simulations). */
    virtual void absorb(int m, int ell, double DeltaL, bool ynstellar) = 0;

    /** This pure virtual function returns the absorbed luminosity \f$L_{\ell,m}\f$ at wavelength index
        \f$\ell\f$ in the dust cell with cell number \f$m\f$. */
    virtual double absorbedLuminosity(int m, int ell) const = 0;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    DustDistribution*& _dd{_dustDistribution};
    DustGrid*& _grid{_dustGrid};

    // data members initialized during setup
    DustGridDensityInterface* _gdi{nullptr};
    ProcessAssigner* _setupAssigner{nullptr};
    int _Ncomp{0};      // cached number of components (index h) in dust distribution
    int _Ncells{0};     // cached number of cells (index m) in dust grid
    Array _volumev;     // volume for each cell (indexed on m)
    Table<2> _rhovv;    // density for each cell and each dust component (indexed on m,h)
    vector<int64_t> _crossed;
    std::mutex _crossedMutex;
};

//////////////////////////////////////////////////////////////////////

#endif
