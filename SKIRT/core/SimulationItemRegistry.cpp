/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SimulationItemRegistry.hpp"
#include "ItemRegistry.hpp"
#include "SkirtUnitDef.hpp"

// ---> add new items below in alphabetical order

#include "AdaptiveMeshAmrvacFile.hpp"
#include "AdaptiveMeshAsciiFile.hpp"
#include "AdaptiveMeshDustDistribution.hpp"
#include "AdaptiveMeshDustGrid.hpp"
#include "AdaptiveMeshGeometry.hpp"
#include "AdaptiveMeshStellarComp.hpp"
#include "AllCellsDustLib.hpp"
#include "AmHydrocarbonGrainComposition.hpp"
#include "Benchmark1DDustMix.hpp"
#include "Benchmark2DDustMix.hpp"
#include "BinTreeDustGrid.hpp"
#include "BlackBodySED.hpp"
#include "BolLuminosityStellarCompNormalization.hpp"
#include "BoxClipGeometryDecorator.hpp"
#include "BrokenExpDiskGeometry.hpp"
#include "BruzualCharlotSED.hpp"
#include "BruzualCharlotSEDFamily.hpp"
#include "CartesianDustGrid.hpp"
#include "ClumpyGeometryDecorator.hpp"
#include "CombineGeometryDecorator.hpp"
#include "CompDustDistribution.hpp"
#include "ConfigurableDustMix.hpp"
#include "ConicalShellGeometry.hpp"
#include "CubBackgroundGeometry.hpp"
#include "CubicSplineSmoothingKernel.hpp"
#include "Cylinder2DDustGrid.hpp"
#include "CylindricalClipGeometryDecorator.hpp"
#include "Dim1DustLib.hpp"
#include "Dim2DustLib.hpp"
#include "DraineGraphiteGrainComposition.hpp"
#include "DraineIonizedPAHGrainComposition.hpp"
#include "DraineLiDustMix.hpp"
#include "DraineNeutralPAHGrainComposition.hpp"
#include "DraineSilicateGrainComposition.hpp"
#include "DustComp.hpp"
#include "DustEmGrainComposition.hpp"
#include "DustMassDustCompNormalization.hpp"
#include "DustMixPopulation.hpp"
#include "EdgeOnDustCompNormalization.hpp"
#include "EinastoGeometry.hpp"
#include "ElectronDustMix.hpp"
#include "EnstatiteGrainComposition.hpp"
#include "ExpDiskGeometry.hpp"
#include "ExtragalacticUnits.hpp"
#include "FaceOnDustCompNormalization.hpp"
#include "FileGrainComposition.hpp"
#include "FileMesh.hpp"
#include "FileSED.hpp"
#include "FileWavelengthGrid.hpp"
#include "FoamGeometryDecorator.hpp"
#include "ForsteriteGrainComposition.hpp"
#include "FrameInstrument.hpp"
#include "FSPSVariableIMFSEDFamily.hpp"
#include "FullInstrument.hpp"
#include "GammaGeometry.hpp"
#include "GaussianGeometry.hpp"
#include "GreyBodyDustEmissivity.hpp"
#include "HyperboloidGeometry.hpp"
#include "HyperboloidShellGeometry.hpp"
#include "InstrumentFrame.hpp"
#include "InstrumentSystem.hpp"
#include "InterstellarDustMix.hpp"
#include "KuruczSED.hpp"
#include "LaserGeometry.hpp"
#include "LinMesh.hpp"
#include "LogMesh.hpp"
#include "LogNormalGrainSizeDistribution.hpp"
#include "LogWavelengthGrid.hpp"
#include "MappingsSEDFamily.hpp"
#include "MarastonSED.hpp"
#include "MappingsSED.hpp"
#include "MeanZubkoDustMix.hpp"
#include "MeshDustComponent.hpp"
#include "MGEGeometry.hpp"
#include "MRNDustMix.hpp"
#include "MieSilicateGrainComposition.hpp"
#include "MinSilicateGrainComposition.hpp"
#include "ModifiedLogNormalGrainSizeDistribution.hpp"
#include "ModifiedPowerLawGrainSizeDistribution.hpp"
#include "MultiFrameInstrument.hpp"
#include "NestedLogWavelengthGrid.hpp"
#include "NetzerAccretionDiskGeometry.hpp"
#include "OctTreeDustGrid.hpp"
#include "OffsetGeometryDecorator.hpp"
#include "OligoDustSystem.hpp"
#include "OligoMonteCarloSimulation.hpp"
#include "OligoStellarComp.hpp"
#include "OligoWavelengthGrid.hpp"
#include "PanDustSystem.hpp"
#include "PanMonteCarloSimulation.hpp"
#include "PanStellarComp.hpp"
#include "PanWavelengthGrid.hpp"
#include "ParaboloidGeometry.hpp"
#include "ParaboloidShellGeometry.hpp"
#include "ParticleTreeDustGrid.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "PegaseSED.hpp"
#include "PerspectiveInstrument.hpp"
#include "PinteBenchmarkDustMix.hpp"
#include "PinteBenchmarkGeometry.hpp"
#include "PlummerGeometry.hpp"
#include "PointGeometry.hpp"
#include "PolarizedGraphiteGrainComposition.hpp"
#include "PolarizedSilicateGrainComposition.hpp"
#include "PowMesh.hpp"
#include "PowerLawGrainSizeDistribution.hpp"
#include "PseudoSersicGeometry.hpp"
#include "QuasarSED.hpp"
#include "RadialDustCompNormalization.hpp"
#include "Random.hpp"
#include "ReadFitsGeometry.hpp"
#include "RingGeometry.hpp"
#include "RotateGeometryDecorator.hpp"
#include "SEDInstrument.hpp"
#include "SIUnits.hpp"
#include "SPHDustDistribution.hpp"
#include "SPHGeometry.hpp"
#include "SPHStellarComp.hpp"
#include "SersicGeometry.hpp"
#include "ShellGeometry.hpp"
#include "SimpleInstrument.hpp"
#include "SimpleOligoDustMix.hpp"
#include "SingleGrainSizeDistribution.hpp"
#include "SmoothingKernel.hpp"
#include "SpectralLuminosityStellarCompNormalization.hpp"
#include "SpheBackgroundGeometry.hpp"
#include "SpheGeometry.hpp"
#include "Sphere1DDustGrid.hpp"
#include "Sphere2DDustGrid.hpp"
#include "SphericalAdaptiveMeshDustDistribution.hpp"
#include "SphericalClipGeometryDecorator.hpp"
#include "SpheroidalGeometryDecorator.hpp"
#include "SpiralStructureGeometryDecorator.hpp"
#include "StarburstSED.hpp"
#include "Starburst99SED.hpp"
#include "Starburst99SEDFamily.hpp"
#include "StellarSurfaceGeometry.hpp"
#include "StellarSystem.hpp"
#include "StellarUnits.hpp"
#include "SunSED.hpp"
#include "SymPowMesh.hpp"
#include "TTauriDiskGeometry.hpp"
#include "ThemisDustMix.hpp"
#include "TorusGeometry.hpp"
#include "TransientDustEmissivity.hpp"
#include "TriaxialGeometryDecorator.hpp"
#include "Trust1Geometry.hpp"
#include "Trust2Geometry.hpp"
#include "Trust6Geometry.hpp"
#include "Trust7aGeometry.hpp"
#include "Trust7bGeometry.hpp"
#include "TrustDustMix.hpp"
#include "TrustGraphiteGrainComposition.hpp"
#include "TrustMeanDustMix.hpp"
#include "TrustNeutralPAHGrainComposition.hpp"
#include "TrustPolarizedMeanDustMix.hpp"
#include "TrustSilicateGrainComposition.hpp"
#include "TwoPhaseDustGrid.hpp"
#include "UniformCuboidGeometry.hpp"
#include "UniformSmoothingKernel.hpp"
#include "VoronoiDustDistribution.hpp"
#include "VoronoiDustGrid.hpp"
#include "VoronoiGeometry.hpp"
#include "VoronoiMeshAsciiFile.hpp"
#include "VoronoiStellarComp.hpp"
#include "WeingartnerDraineDustMix.hpp"
#include "XDustCompNormalization.hpp"
#include "YDustCompNormalization.hpp"
#include "ZDustCompNormalization.hpp"
#include "ZubkoDustMix.hpp"
#include "ZubkoGraphiteGrainSizeDistribution.hpp"
#include "ZubkoPAHGrainSizeDistribution.hpp"
#include "ZubkoSilicateGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

SimulationItemRegistry::SimulationItemRegistry(string version, string format)
{
    // start a new schema
    ItemRegistry::beginSchema("SKIRT", "a SKIRT parameter file", version, "ski",
                              "skirt-simulation-hierarchy", "MonteCarloSimulation", format,
                              "http://www.skirt.ugent.be/skirt");

    // add the SKIRT unit definitions
    ItemRegistry::addUnitDef<SkirtUnitDef>();

    // add the SKIRT simulation items
    ItemRegistry::add<SimulationItem>();

    // ---> add new items in the order you want them to appear in choice lists for the user

    // basic building blocks
    ItemRegistry::add<Simulation>();
    ItemRegistry::add<Random>();
    ItemRegistry::add<Units>();
    ItemRegistry::add<SIUnits>();
    ItemRegistry::add<StellarUnits>();
    ItemRegistry::add<ExtragalacticUnits>();

    // Monte Carlo simulations
    ItemRegistry::add<MonteCarloSimulation>();
    ItemRegistry::add<OligoMonteCarloSimulation>();
    ItemRegistry::add<PanMonteCarloSimulation>();

    // instruments
    ItemRegistry::add<InstrumentSystem>();
    ItemRegistry::add<Instrument>();
    ItemRegistry::add<DistantInstrument>();
    ItemRegistry::add<SingleFrameInstrument>();
    ItemRegistry::add<SEDInstrument>();
    ItemRegistry::add<FrameInstrument>();
    ItemRegistry::add<SimpleInstrument>();
    ItemRegistry::add<FullInstrument>();
    ItemRegistry::add<PerspectiveInstrument>();
    ItemRegistry::add<MultiFrameInstrument>();
    ItemRegistry::add<InstrumentFrame>();

    // wavelength grids
    ItemRegistry::add<WavelengthGrid>();
    ItemRegistry::add<OligoWavelengthGrid>();
    ItemRegistry::add<PanWavelengthGrid>();
    ItemRegistry::add<LogWavelengthGrid>();
    ItemRegistry::add<NestedLogWavelengthGrid>();
    ItemRegistry::add<FileWavelengthGrid>();

    // stellar systems
    ItemRegistry::add<StellarSystem>();
    ItemRegistry::add<StellarComp>();
    ItemRegistry::add<BoxStellarComp>();
    ItemRegistry::add<GeometricStellarComp>();
    ItemRegistry::add<OligoStellarComp>();
    ItemRegistry::add<PanStellarComp>();
    ItemRegistry::add<StellarCompNormalization>();
    ItemRegistry::add<BolLuminosityStellarCompNormalization>();
    ItemRegistry::add<SpectralLuminosityStellarCompNormalization>();
    ItemRegistry::add<SPHStellarComp>();
    ItemRegistry::add<AdaptiveMeshStellarComp>();
    ItemRegistry::add<VoronoiStellarComp>();

    // isotropic geometries
    ItemRegistry::add<Geometry>();
    ItemRegistry::add<PointGeometry>();
    ItemRegistry::add<SpheGeometry>();
    ItemRegistry::add<PlummerGeometry>();
    ItemRegistry::add<GammaGeometry>();
    ItemRegistry::add<SersicGeometry>();
    ItemRegistry::add<PseudoSersicGeometry>();
    ItemRegistry::add<EinastoGeometry>();
    ItemRegistry::add<GaussianGeometry>();
    ItemRegistry::add<ShellGeometry>();
    ItemRegistry::add<AxGeometry>();
    ItemRegistry::add<SepAxGeometry>();
    ItemRegistry::add<ExpDiskGeometry>();
    ItemRegistry::add<BrokenExpDiskGeometry>();
    ItemRegistry::add<MGEGeometry>();
    ItemRegistry::add<TTauriDiskGeometry>();
    ItemRegistry::add<RingGeometry>();
    ItemRegistry::add<TorusGeometry>();
    ItemRegistry::add<ConicalShellGeometry>();
    ItemRegistry::add<ParaboloidGeometry>();
    ItemRegistry::add<ParaboloidShellGeometry>();
    ItemRegistry::add<HyperboloidGeometry>();
    ItemRegistry::add<HyperboloidShellGeometry>();
    ItemRegistry::add<GenGeometry>();
    ItemRegistry::add<BoxGeometry>();
    ItemRegistry::add<UniformCuboidGeometry>();
    ItemRegistry::add<SPHGeometry>();
    ItemRegistry::add<AdaptiveMeshGeometry>();
    ItemRegistry::add<VoronoiGeometry>();
    ItemRegistry::add<ReadFitsGeometry>();

    // anistropic geometries
    ItemRegistry::add<LaserGeometry>();
    ItemRegistry::add<NetzerAccretionDiskGeometry>();
    ItemRegistry::add<StellarSurfaceGeometry>();
    ItemRegistry::add<SpheBackgroundGeometry>();
    ItemRegistry::add<CubBackgroundGeometry>();

    // benchmark geometries
    ItemRegistry::add<PinteBenchmarkGeometry>();
    ItemRegistry::add<Trust1Geometry>();
    ItemRegistry::add<Trust2Geometry>();
    ItemRegistry::add<Trust6Geometry>();
    ItemRegistry::add<Trust7aGeometry>();
    ItemRegistry::add<Trust7bGeometry>();

    // geometry decorators
    ItemRegistry::add<OffsetGeometryDecorator>();
    ItemRegistry::add<RotateGeometryDecorator>();
    ItemRegistry::add<SpheroidalGeometryDecorator>();
    ItemRegistry::add<TriaxialGeometryDecorator>();
    ItemRegistry::add<ClipGeometryDecorator>();
    ItemRegistry::add<SphericalClipGeometryDecorator>();
    ItemRegistry::add<CylindricalClipGeometryDecorator>();
    ItemRegistry::add<BoxClipGeometryDecorator>();
    ItemRegistry::add<SpiralStructureGeometryDecorator>();
    ItemRegistry::add<ClumpyGeometryDecorator>();
    ItemRegistry::add<CombineGeometryDecorator>();
    ItemRegistry::add<FoamGeometryDecorator>();

    // smoothing kernels
    ItemRegistry::add<SmoothingKernel>();
    ItemRegistry::add<UniformSmoothingKernel>();
    ItemRegistry::add<CubicSplineSmoothingKernel>();

    // dust systems
    ItemRegistry::add<DustSystem>();
    ItemRegistry::add<OligoDustSystem>();
    ItemRegistry::add<PanDustSystem>();

    // dust components and corresponding normalizations
    ItemRegistry::add<DustComp>();
    ItemRegistry::add<MeshDustComponent>();

    // dust component normalizations
    ItemRegistry::add<DustCompNormalization>();
    ItemRegistry::add<DustMassDustCompNormalization>();
    ItemRegistry::add<RadialDustCompNormalization>();
    ItemRegistry::add<FaceOnDustCompNormalization>();
    ItemRegistry::add<EdgeOnDustCompNormalization>();
    ItemRegistry::add<XDustCompNormalization>();
    ItemRegistry::add<YDustCompNormalization>();
    ItemRegistry::add<ZDustCompNormalization>();

    // dust distributions
    ItemRegistry::add<DustDistribution>();
    ItemRegistry::add<BoxDustDistribution>();
    ItemRegistry::add<CompDustDistribution>();
    ItemRegistry::add<SPHDustDistribution>();
    ItemRegistry::add<AdaptiveMeshDustDistribution>();
    ItemRegistry::add<SphericalAdaptiveMeshDustDistribution>();
    ItemRegistry::add<VoronoiDustDistribution>();

    // mesh file representations
    ItemRegistry::add<AdaptiveMeshFile>();
    ItemRegistry::add<AdaptiveMeshAsciiFile>();
    ItemRegistry::add<AdaptiveMeshAmrvacFile>();
    ItemRegistry::add<VoronoiMeshFile>();
    ItemRegistry::add<VoronoiMeshAsciiFile>();

    // meshes for the dust grids
    ItemRegistry::add<Mesh>();
    ItemRegistry::add<MoveableMesh>();
    ItemRegistry::add<AnchoredMesh>();
    ItemRegistry::add<LinMesh>();
    ItemRegistry::add<PowMesh>();
    ItemRegistry::add<SymPowMesh>();
    ItemRegistry::add<LogMesh>();
    ItemRegistry::add<FileMesh>();

    // dust grids
    ItemRegistry::add<DustGrid>();
    ItemRegistry::add<SphereDustGrid>();
    ItemRegistry::add<Sphere1DDustGrid>();
    ItemRegistry::add<Sphere2DDustGrid>();
    ItemRegistry::add<CylinderDustGrid>();
    ItemRegistry::add<Cylinder2DDustGrid>();
    ItemRegistry::add<BoxDustGrid>();
    ItemRegistry::add<CartesianDustGrid>();
    ItemRegistry::add<TwoPhaseDustGrid>();
    ItemRegistry::add<TreeDustGrid>();
    ItemRegistry::add<OctTreeDustGrid>();
    ItemRegistry::add<BinTreeDustGrid>();
    ItemRegistry::add<ParticleTreeDustGrid>();
    ItemRegistry::add<VoronoiDustGrid>();
    ItemRegistry::add<AdaptiveMeshDustGrid>();

    // dust mixtures
    ItemRegistry::add<DustMix>();
    ItemRegistry::add<InterstellarDustMix>();
    ItemRegistry::add<DraineLiDustMix>();
    ItemRegistry::add<MeanZubkoDustMix>();
    ItemRegistry::add<SimpleOligoDustMix>();
    ItemRegistry::add<Benchmark1DDustMix>();
    ItemRegistry::add<Benchmark2DDustMix>();
    ItemRegistry::add<PinteBenchmarkDustMix>();
    ItemRegistry::add<TrustMeanDustMix>();
    ItemRegistry::add<TrustPolarizedMeanDustMix>();
    ItemRegistry::add<MultiGrainDustMix>();
    ItemRegistry::add<MRNDustMix>();
    ItemRegistry::add<WeingartnerDraineDustMix>();
    ItemRegistry::add<ZubkoDustMix>();
    ItemRegistry::add<TrustDustMix>();
    ItemRegistry::add<ThemisDustMix>();
    ItemRegistry::add<ConfigurableDustMix>();
    ItemRegistry::add<ElectronDustMix>();

    // grain compositions
    ItemRegistry::add<GrainComposition>();
    ItemRegistry::add<DraineGraphiteGrainComposition>();
    ItemRegistry::add<DraineSilicateGrainComposition>();
    ItemRegistry::add<DraineNeutralPAHGrainComposition>();
    ItemRegistry::add<DraineIonizedPAHGrainComposition>();
    ItemRegistry::add<MieSilicateGrainComposition>();
    ItemRegistry::add<MinSilicateGrainComposition>();
    ItemRegistry::add<AmHydrocarbonGrainComposition>();
    ItemRegistry::add<EnstatiteGrainComposition>();
    ItemRegistry::add<ForsteriteGrainComposition>();
    ItemRegistry::add<PolarizedGraphiteGrainComposition>();
    ItemRegistry::add<PolarizedSilicateGrainComposition>();
    ItemRegistry::add<TrustGraphiteGrainComposition>();
    ItemRegistry::add<TrustSilicateGrainComposition>();
    ItemRegistry::add<TrustNeutralPAHGrainComposition>();
    ItemRegistry::add<DustEmGrainComposition>();
    ItemRegistry::add<FileGrainComposition>();

    // grain size distributions
    ItemRegistry::add<GrainSizeDistribution>();
    ItemRegistry::add<RangeGrainSizeDistribution>();
    ItemRegistry::add<PowerLawGrainSizeDistribution>();
    ItemRegistry::add<ModifiedPowerLawGrainSizeDistribution>();
    ItemRegistry::add<LogNormalGrainSizeDistribution>();
    ItemRegistry::add<ModifiedLogNormalGrainSizeDistribution>();
    ItemRegistry::add<SingleGrainSizeDistribution>();
    ItemRegistry::add<ZubkoGraphiteGrainSizeDistribution>();
    ItemRegistry::add<ZubkoSilicateGrainSizeDistribution>();
    ItemRegistry::add<ZubkoPAHGrainSizeDistribution>();

    // dust mix population
    ItemRegistry::add<DustMixPopulation>();

    // dust emissivity calculators
    ItemRegistry::add<DustEmissivity>();
    ItemRegistry::add<GreyBodyDustEmissivity>();
    ItemRegistry::add<TransientDustEmissivity>();

    // dust libraries
    ItemRegistry::add<DustLib>();
    ItemRegistry::add<AllCellsDustLib>();
    ItemRegistry::add<Dim1DustLib>();
    ItemRegistry::add<Dim2DustLib>();

    // stellar SEDs
    ItemRegistry::add<SED>();
    ItemRegistry::add<StellarSED>();
    ItemRegistry::add<BlackBodySED>();
    ItemRegistry::add<SunSED>();
    ItemRegistry::add<PegaseSED>();
    ItemRegistry::add<BruzualCharlotSED>();
    ItemRegistry::add<Starburst99SED>();
    ItemRegistry::add<MarastonSED>();
    ItemRegistry::add<MappingsSED>();
    ItemRegistry::add<StarburstSED>();
    ItemRegistry::add<KuruczSED>();
    ItemRegistry::add<QuasarSED>();
    ItemRegistry::add<FileSED>();

    // SED families
    ItemRegistry::add<SEDFamily>();
    ItemRegistry::add<BruzualCharlotSEDFamily>();
    ItemRegistry::add<Starburst99SEDFamily>();
    ItemRegistry::add<FSPSVariableIMFSEDFamily>();
    ItemRegistry::add<MappingsSEDFamily>();
}

////////////////////////////////////////////////////////////////////

const SchemaDef* SimulationItemRegistry::getSchemaDef()
{
    return ItemRegistry::getSchemaDef("SKIRT");
}

////////////////////////////////////////////////////////////////////

SimulationItemRegistry::~SimulationItemRegistry()
{
    ItemRegistry::finalize();
}

////////////////////////////////////////////////////////////////////

