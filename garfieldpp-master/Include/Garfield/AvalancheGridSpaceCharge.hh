#ifndef GARFIELD_AVALANCHEGRIDSPACECHARGE_HH
#define GARFIELD_AVALANCHEGRIDSPACECHARGE_HH

#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "AvalancheMicroscopic.hh"
#include "ComponentParallelPlate.hh"
#include "GarfieldConstants.hh"
#include "Sensor.hh"

namespace Garfield {
/// Propagates avalanches with the 2d (axi-symmetric) space-charge routine from
/// Lippmann, Riegler (2004) in uniform background fields. Different options to
/// calculate space-charge-fields can be chosen.
class AvalancheGridSpaceCharge {
 public:
  /// Constructor
  AvalancheGridSpaceCharge();

  /// Destructor
  ~AvalancheGridSpaceCharge() = default;

  /// Reset grid i.e. ImportElectrons can be called again.
  void Reset();

  /// Enable/disable debugging ( = log messages) (default off)
  void EnableDebugging(const bool option = true) { m_bDebug = option; }

  /// Enable sticky anode (default on)
  void EnableStickyAnode(const bool option = true) { m_bStick = option; }

  /// Enable to use TOF swarm parameters (default on)
  void EnableTOF(const bool option = true) { m_bUseTOF = option; }

  /// Enable diffusion (default off)
  void EnableDiffusion(const bool option = true) { m_bDiffusion = option; }

  /// Enable space charge calculations (default on)
  void EnableSpaceChargeEffect(const bool option = true) {
    m_bSpaceCharge = option;
  }

  /// Enable adaptive time stepping (default on)
  void EnableAdaptiveTimeStepping(const bool option = true) {
    m_bAdaptiveTime = option;
  }

  /// Enable Monte Carlo (gain/diffusion) up to a certain total number of
  /// electrons e.g. (1e5) (default on)
  ///   if disabled: Mean values are considered instead of the statistical
  ///   processes
  void EnableMC(const bool option = true) { m_bMC = option; }

  /// Set electron saturation value applied for each gas gap if space charge
  /// effect is turned off (default 1e8)
  void SetNCrit(const long NCrit = 1e8) { m_lNCrit = NCrit; }

  /// Sets the different options to calculate the space charge field
  ///   coulomb: free field approximation
  ///   (relaxation: fdm-relaxation methode (with a condition to stop))
  ///   mirror: symmetric 3 layer single gap rpc with metal - resistive layer -
  ///   gas gap - r. l. - m.
  void SetFieldCalculation(const std::string &option = "coulomb",
                           const int nof_approx = 1.) {
    m_sFieldOption = std::move(option);
    m_iFieldApprox = std::move(nof_approx);
  }

  /// Set the streamer-inception criterion constant K in the interval (0,
  /// &infin;) s.t. 1 = 100%
  void SetK(float option = 0.95) { m_fStreamerK = option; }

  /// Stop the avalanche if K % field is reached
  void SetStopAtK(bool option = true) { m_bStopAtK = option; }

  /// Set the sensor (+ determines if base Cmp is CmpParallelPlate (MRPCS)).
  void SetSensor(Sensor *sensor) {
    m_sensor = sensor;
    // Determine if any component is CmpParallelPlate (if not it will stay
    // nullptr).
    m_ParallelPlate = nullptr;
    const size_t nofCmp = m_sensor->GetNumberOfComponents();
    for (size_t i = 0; i < nofCmp; i++) {
      if (!m_ParallelPlate) {
        m_ParallelPlate =
            dynamic_cast<ComponentParallelPlate *>(m_sensor->GetComponent(i));
      }
    }
  }

  /**
   *
   * @param zmin coordinate along direction of background field
   * @param zmax should more or less be at the end of the gap
   * @param zsteps number of steps
   * @param rmax coordinate perpendicular to background field
   * @param rsteps number of steps
   */
  void Set2dGrid(double zmin, double zmax, int zsteps, double rmax, int rsteps);

  /// Import electron (no ions) data from AvalancheMicroscopic class to
  /// axi-symmetric grid
  void ImportElectronsFromAvalancheMicroscopic(AvalancheMicroscopic *avmc);

  /// Set n electrons onto the grid
  void AvalancheElectron(double x, double y, double z, double t = 0, int n = 1);

  /// (After calling AvalancheElectron) add more electrons on the same
  /// transversal line (y-freedom)
  void AddExtraAvalancheElectron(double y, int n = 1);

  /// Starts the simulation with the imported electrons for a time step dt
  /// (dt = -1 until there are no electrons left in the gap).
  void StartGridAvalanche(double dtime = -1);

  /// Returns the total positive charge in the gap's
  long GetAvalancheSize() { return m_nTotPosIons; }

  /// Return current mean distance of the electrons on the grid
  double GetMeanDistance();

  /// Returns if 100 * K % background charge has been reached
  [[nodiscard]] long ReachedKPercent() const {
    if (m_bFieldK)
      return m_lElectronsK;
    else
      return -1;
  }

  /// Returns the total electron number evolution
  [[nodiscard]] const std::vector<std::pair<double, long>>
  GetElectronEvolution() const {
    return m_vNElectronEvolution;
  }

  /// Export the current grid to a txt file (electron, ions numbers and field
  /// magnitude)
  ///  Take care: only where electrons are located is the space-charge field
  ///  evaluated!
  void ExportGrid(const std::string &filename);

 private:
  struct GridNode {
    long nElectron = 0;  ///< electrons on node
    double nPosIon = 0;  ///< pos ion on node (smeared values allowed)
    double nNegIon = 0;  ///< neg ion on node (smeared values allowed)
    // holder memories for stepping in time:
    long nElectronHolder = 0;  ///< at t+dt
    double nPosIonHolder = 0;  ///< at t+dt
    double nNegIonHolder = 0;  ///< at t+dt

    double townsend = 0;    ///< townsend at this node 1/cm
    double attachment = 0;  ///< attachment at this node 1/cm
    /// Magnitude of velocity of the node (not negative) cm/ns
    double velocity = 0.;
    /// Diffusion along E.
    double dSigmaL = 0;  
    /// Diffusion transverse to E (radial, phi dir is net 0).
    double dSigmaT = 0;          
   
    double Wv = 0;  ///< flux drift cm/ns
    double Wr = 0;  ///< bulk drift cm/ns
    /// Ionization rate from TOF experiment 1/ns -> 1/cm
    double townsendPT = 0;  
    /// Attachment rate from TOF experiment 1/ns -> 1/cm
    double attachmentPT = 0;  
    /// Space-charge electric field in R direction (can be negative)
    double eFieldR = 0;  
    /// Space-charge electric field in Z direction (can be negative)
    double eFieldZ = 0;  

    double time = 0.;  ///< Node clock.

    bool anode = false;  ///< init the anode
    /// LayerIndex in ParallelPlate convention != gas gap index
    int layerIndex = 0;  
    /// Gas gap index: -1 if not gas gap; starts with 0, 1, ...
    int gasGapIndex = 0;  
    bool isGasGap = true;
  };

  struct Point {
    double x, y, z;  ///< coordinates
    double t;        ///< time

    int gasLayerIndex;
  };

  // Prepare grid and place stored electrons from AvalancheMicroscopic import
  void PrepareElectronsFromMicroscopicAvalanche();

  // Assign electron to the closest grid point
  bool SnapTo2dGrid(double x, double y, double z, long n = 1, int gasLayer = 0);

  // Prepare the mesh with the ComponentParallelPlate
  void Prepare2dMesh();

  // Transports the electrons/nodes a timestep
  bool TransportTimeStep();

  // Diffuses the electrons/nodes a timestep
  void DiffuseTimeStep(double dx, long nElectron, double nPosIon,
                       double nNegIon, int iz, int ir, int gasGap);

  // Redistributes the charges
  void DistributeCharges(long nElectron, double nPosIon, double nNegIon, int iz,
                         int ir, double stepZ, double stepR, int gasGap);

  // Calculate the field from all the contributions to the bin of interest. May
  // need much more functionalities/tables.
  void GetLocalField(int iz, int ir, double &eFieldZ, double &eFieldR,
                     const std::string &fieldOption, int gasGap);

  // Calculate the field of charged ring in vacuum using coulomb potentials and
  // indices
  void GetFreeChargedRing(int iz, int ir, int fz, int fr, double &eFieldZ,
                          double &eFieldR);

  // Calculate the field of charged ring in vacuum using coulomb potentials and
  // coordinates
  void GetFreeChargedRing(double zi, double ri, double zf, double rf,
                          double &eFieldZ, double &eFieldR);

  // Get field at (zi, ri) from N charges at (zf, rf) either as a ring or a
  // coulomb ball (rf = 0) if i and f are too close it is considered as self
  // interaction and not included
  bool AddFieldFromChargeAt(int iz, int ir, int fz, int fr, double N,
                            double &eFieldZ, double &eFieldR);

  // Get field at (zi, ri) from N charges at (zf, rf) either as a ring or a
  // coulomb ball (rf = 0) if i and f are too close it is considered as self
  // interaction and not included
  bool AddFieldFromChargeAt(int iz, int ir, double zf, double rf, double N,
                            double &eFieldZ, double &eFieldR);

  // Get swarm parameters at electric field magnitude
  void GetSwarmParameters(double MagEField, double &alpha, double &eta,
                          double &drift, double &dSigmaL, double &dSigmaT,
                          double &wv, double &wr, double &alphaPT,
                          double &etaPT, int gasGap);

  // Change from 2dGrid to Global coordinates
  void GetGlobalCoordinates(double r, double z, double phi, double &xg,
                            double &yg, double &zg, int gasGap);

  // Import elliptic integral values
  void ImportEllipticIntegralValues(const std::string &filename);

  // Gets elliptic integrals via list
  void GetEllipticIntegrals(double x, double &K, double &E);

  // Get from index the gas gap number, else -1
  int GetGasGapNumber(int layerIndex) {
    auto it =
        std::find(m_vIndexGasGaps.begin(), m_vIndexGasGaps.end(), layerIndex);
    return (it != m_vIndexGasGaps.end())
               ? std::distance(m_vIndexGasGaps.begin(), it)
               : -1;
  }

 private:
  std::string m_className = "AvalancheGridSpaceCharge";

  bool m_bDebug = false;
  bool m_bDiffusion = false;
  bool m_bStick = true;
  // boolean for AvalancheElectron
  bool m_bDriftAvalanche = false;
  // boolean for ImportElectronsFromAvalancheMicroscopic
  bool m_bImportAvalanche = false;
  bool m_bPreparedImportAvalanche = false;
  long m_lNCrit = 1e8;
  bool m_bSpaceCharge = true;

  float m_fStreamerK = 0.95;
  bool m_bStopAtK = false;
  bool m_bFieldK = false;
  long m_lElectronsK;

  bool m_bAdaptiveTime = true;
  bool m_bImportElliptic = false;
  /// Flag if TOF parameters should be used, else Magboltz
  /// drift and SST spatial coefficients
  bool m_bUseTOF = true;  
  /// Flag if bulk drift velocity is available to the simulation
  bool m_bWrAvailable = true;  
  /// Flag if temporal rates are available to the simulation
  bool m_bRatesAvailable = true;  

  bool m_bMC = true;

  int m_iFieldApprox = 1;    //< order of approximation in Set(1,2,3,...)
  double m_dMinGroups = 50;  //< same values as lippmann

  ComponentParallelPlate *m_ParallelPlate = nullptr;
  Sensor *m_sensor = nullptr;

  std::vector<double> m_zGrid;  ///< Grid points of z-coordinate.
  int m_zSteps = 0;             ///< Number of grid points.
  double m_zStepSize = 0.;      /// Distance between the grid points.

  std::vector<double> m_rGrid;  ///< Grid points of r-coordinate.
  int m_rSteps = 0.;            ///< Number of grid points
  double m_rStepSize = 0.;      /// Distance between the grid points.

  bool m_isgridset = false;  ///< Keeps track if the grid has been defined.
  long m_nTotElectron = 0;   ///< Total amount of electrons at time step.
  long m_nTotPosIons = 0;    ///< total amount of charge created

  double m_time = 0.;   ///< Clock.
  double m_time0 = 0.;  ///< Initial time.
  double m_dt = 0.;     ///< Time step.
  /// Tracking if the charges are still in the drift gap.
  bool m_run = true;

  std::vector<std::vector<int>> m_zGasGapBoundaries;  ///< [k] -> {izLeft, ..., izRight}

  std::vector<std::vector<GridNode>> m_grid;  ///< grid with nodes on it
  /// Electrons to transfer onto grid
  std::vector<std::vector<Point>> m_vElectrons;  

  std::vector<std::pair<double, long>> m_vNElectronEvolution;
  std::vector<long> m_vGroupSizes = {
      1500, 800, 400, 200, 100,
      50,   20,  10,  5,   2};  ///< same values as Lippmann & Riegler
  /// Which layer indices are gas layers
  std::vector<int> m_vIndexGasGaps = {0};  

  /// Coordinates of center of electron number.
  /// Required: y in [zmin, zmax]
  std::vector<std::vector<double>> m_vCoNGasLayer{};  
  /// Example point (y-coord) in each gas gap
  std::vector<double> m_vYPointInGasGap{};  

  /// Uniform background field in z direction, can be negative.
  std::vector<double> m_ezBkg = {0};
  /// Which gas gaps are saturated if saturation is on
  std::vector<int> m_vSaturatedGaps{};  

  std::string m_sFieldOption = "coulomb";
  std::vector<double> m_vXElliptic;
  std::vector<double> m_vKElliptic;
  std::vector<double> m_vEElliptic;
};

}  // namespace Garfield

#endif  // GARFIELD_AVALANCHEGRIDSPACECHARGE_HH
