#ifndef G_AVALANCHE_MICROSCOPIC_H
#define G_AVALANCHE_MICROSCOPIC_H

#include <string>
#include <vector>

#include <TH1.h>

#include "GarfieldConstants.hh"
#include "Sensor.hh"
#include "ViewDrift.hh"
#include "GPUInterface.hh"

#include "MultiProcessInterface.hh"

namespace Garfield {

  class AvalancheMicroscopicGPU;

/// Calculate electron drift lines and avalanches using microscopic tracking.

class AvalancheMicroscopic {
 public:
  /// Default constructor
  AvalancheMicroscopic() : AvalancheMicroscopic(nullptr) {}
  /// Constructor
  AvalancheMicroscopic(Sensor* sensor);
  /// Destructor
  ~AvalancheMicroscopic() {}

  /// Set the sensor.
  void SetSensor(Sensor* sensor);

  /// Switch on drift line plotting.
  void EnablePlotting(ViewDrift* view, const size_t nColl = 100);
  /// Switch off drift line plotting.
  void DisablePlotting() { m_viewer = nullptr; }
  /// Draw a marker at every excitation or not.
  void EnableExcitationMarkers(const bool on = true) { m_plotExcitations = on; }
  /// Draw a marker at every ionising collision or not.
  void EnableIonisationMarkers(const bool on = true) { m_plotIonisations = on; }
  /// Draw a marker at every attachment or not.
  void EnableAttachmentMarkers(const bool on = true) { m_plotAttachments = on; }

  /// Switch calculation of induced currents on or off (default: enabled).
  void EnableSignalCalculation(const bool on = true) { m_doSignal = on; }
  /// Use the weighting potential (as opposed to the weighting field)
  /// for calculating the induced current.
  void UseWeightingPotential(const bool on = true) {
    m_useWeightingPotential = on;
  }
  /// Integrate the weighting field over a drift line step when
  /// calculating the induced current (default: off).
  void EnableWeightingFieldIntegration(const bool on = true) {
    m_integrateWeightingField = on;
  }

  /// Switch on calculation of the total induced charge (default: off).
  void UseInducedCharge(const bool on = true) { m_doInducedCharge = on; }

  /// Compute and store the path length of each trajectory (default: off).
  void EnablePathLengthComputation(const bool on = true) {
    m_computePathLength = on;
  } 
  /// Fill a histogram with the electron energy distribution.
  void EnableElectronEnergyHistogramming(TH1* histo);
  /// Stop histogramming the electron energy distribution.
  void DisableElectronEnergyHistogramming() { m_histElectronEnergy = nullptr; }
  /// Fill a histogram with the hole energy distribution.
  void EnableHoleEnergyHistogramming(TH1* histo);
  /// Stop histogramming the hole energy distribution.
  void DisableHoleEnergyHistogramming() { m_histHoleEnergy = nullptr; }

  /** Fill histograms of the distance between successive collisions.
   * \param histo
   pointer to the histogram to be filled
   * \param opt
   direction ('x', 'y', 'z', 'r')
   */
  void SetDistanceHistogram(TH1* histo, const char opt = 'r');
  /// Fill distance distribution histograms for a given collision type.
  void EnableDistanceHistogramming(const int type);
  /// Stop filling distance distribution histograms for a given collision type.
  void DisableDistanceHistogramming(const int type);
  /// Stop filling distance distribution histograms.
  void DisableDistanceHistogramming();
  /// Fill histograms of the energy of electrons emitted in ionising collisions.
  void EnableSecondaryEnergyHistogramming(TH1* histo);
  /// Stop histogramming the secondary electron energy distribution.
  void DisableSecondaryEnergyHistogramming() { m_histSecondary = nullptr; }

  /// Switch on storage of drift lines (default: off).
  void EnableDriftLines(const bool on = true) { m_storeDriftLines = on; }

  /** Switch on photon transport.
   * \remark This feature has not been tested thoroughly. */
  void EnablePhotonTransport(const bool on = true) { m_usePhotons = on; }

  /// Switch on stepping according to band structure E(k), for semiconductors.
  void EnableBandStructure(const bool on = true) {
    m_useBandStructure = on;
  }

  /// Switch on update of coordinates for null-collision steps (default: off).
  void EnableNullCollisionSteps(const bool on = true, const int nSteps = 1) {
    m_useNullCollisionSteps = on;
    m_nullCollScale = nSteps <= 1 ? 1. : 1. / nSteps;
  }

  /// Switch on Runge-Kutta-Nystrom stepping (default: off).
  void EnableRKNSteps(const bool on = true) {
    m_rknSteps = on;
  }
  /// Set error tolerance and minimum step size on Runge-Kutta-Nystrom method 
  /// (default: 1.e-10 and 1.e-5).
  void SetRKNTolerance(const double sTol = 1.e-10,
                       const double sMinStep = 1.e-5) {
    m_rknsteperrortol = sTol; 
    m_rknMinh = sMinStep;
  }

  /** Set a (lower) energy threshold for electron transport.
   * This can be useful for simulating delta electrons. */
  void SetElectronTransportCut(const double cut) { m_deltaCut = cut; }
  /// Retrieve the value of the energy threshold.
  double GetElectronTransportCut() const { return m_deltaCut; }

  /// Set an energy threshold for photon transport.
  void SetPhotonTransportCut(const double cut) { m_gammaCut = cut; }
  /// Retrieve the energy threshold for transporting photons.
  double GetPhotonTransportCut() const { return m_gammaCut; }

  /** Set a max. avalanche size (i. e. ignore ionising collisions
   once this size has been reached). */
  void EnableAvalancheSizeLimit(const unsigned int size) { m_sizeCut = size; }
  /// Do not apply a limit on the avalanche size.
  void DisableAvalancheSizeLimit() { m_sizeCut = 0; }
  /// Retrieve the currently set size limit.
  int GetAvalancheSizeLimit() const { return m_sizeCut; }

  /// Switch on/off using the magnetic field in the stepping algorithm.
  void EnableMagneticField(const bool on = true) { 
    m_useBfieldAuto = false;
    m_useBfield = on; 
  }

  /// Set number of collisions to be skipped for storing drift lines.
  void SetCollisionSteps(const unsigned int n) { m_nCollSkip = n; }

  /// Define a time interval (only carriers inside the interval are simulated).
  void SetTimeWindow(const double t0, const double t1);
  /// Do not restrict the time interval within which carriers are simulated.
  void UnsetTimeWindow() { m_hasTimeWindow = false; }

  /// Return the number of electrons and ions in the avalanche.
  void GetAvalancheSize(int& ne, int& ni) const {
    ne = m_nElectrons;
    ni = m_nIons;
  }
  void GetAvalancheSize(int& ne, int& nh, int& ni) const {
    ne = m_nElectrons;
    nh = m_nHoles;
    ni = m_nIons;
  }
  /// Return the number of electrons and ions in the avalanche.
  void GetAvalancheSizeGPU(int& ne, int& ni) const {
    ne = m_nElectronsGPU;
    ni = m_nIonsGPU;
  }

  struct Point {
    double x, y, z;    ///< Coordinates.
    double t;          ///< Time.
    double energy;     ///< Kinetic energy.
    double kx, ky, kz; ///< Direction/wave vector.
    int band;          ///< Band.
  };

  struct Electron {
    int status = 0;                ///< Status.
    std::vector<Point> path;       ///< Drift line.
    double pathLength = 0.;        ///< Path length.
  };

  const std::vector<Electron>& GetElectrons() const { return m_electrons; }
  const std::vector<Electron>& GetHoles() const { return m_holes; }
  /** Return the number of electron trajectories in the last
   * simulated avalanche (including captured electrons). */
  size_t GetNumberOfElectronEndpoints() const {
    return m_electrons.size();
  }
  unsigned int GetNumberOfElectronEndpointsGPU() const {
    return m_electrons_gpu.size();
  }
  /** Return the coordinates and time of start and end point of a given
   * electron drift line.
   * \param i index of the drift line
   * \param x0,y0,z0,t0 coordinates and time of the starting point
   * \param x1,y1,z1,t1 coordinates and time of the end point
   * \param e0,e1 initial and final energy
   * \param status status code (see GarfieldConstants.hh)
   */
  void GetElectronEndpoint(const size_t i, double& x0, double& y0,
                           double& z0, double& t0, double& e0, double& x1,
                           double& y1, double& z1, double& t1, double& e1,
                           int& status) const;
  void GetElectronEndpointGPU(const size_t i, double& x0, double& y0,
                           double& z0, double& t0, double& e0, double& x1,
                           double& y1, double& z1, double& t1, double& e1,
                           int& status) const;
  void GetElectronEndpoint(const size_t i, double& x0, double& y0,
                           double& z0, double& t0, double& e0, double& x1,
                           double& y1, double& z1, double& t1, double& e1,
                           double& dx1, double& dy1, double& dz1,
                           int& status) const;
  size_t GetNumberOfElectronDriftLinePoints(const size_t i = 0) const;
  size_t GetNumberOfHoleDriftLinePoints(const size_t i = 0) const;
  void GetElectronDriftLinePoint(double& x, double& y, double& z, double& t,
                                 const size_t ip, const size_t ie = 0) const;
  void GetHoleDriftLinePoint(double& x, double& y, double& z, double& t,
                             const size_t ip, const size_t ih = 0) const;

  size_t GetNumberOfHoleEndpoints() const { return m_holes.size(); }
  void GetHoleEndpoint(const size_t i, double& x0, double& y0, double& z0,
                       double& t0, double& e0, double& x1, double& y1,
                       double& z1, double& t1, double& e1, int& status) const;

  size_t GetNumberOfPhotons() const { return m_photons.size(); }
  // Status codes:
  //   -2: photon absorbed by gas molecule
  void GetPhoton(const size_t i, double& e, double& x0, double& y0,
                 double& z0, double& t0, double& x1, double& y1, double& z1,
                 double& t1, int& status) const;

  /** Calculate an electron drift line.
   * \param x,y,z,t starting point of the electron
   * \param e initial energy of the electron
   * \param dx,dy,dz initial direction vector of the electron
   * If the initial direction is not specified, it is sampled randomly.
   * Secondary electrons are not transported. */
  bool DriftElectron(const double x, const double y, const double z,
                     const double t, const double e, const double dx = 0.,
                     const double dy = 0., const double dz = 0.);

  /// Calculate an avalanche initiated by a given electron.
  bool AvalancheElectron(const double x, const double y, const double z,
                         const double t, const double e,
                         const double dx = 0., const double dy = 0.,
                         const double dz = 0.);
  /// Add an electron to the list of particles to be transported.
  void AddElectron(const double x, const double y, const double z,
                   const double t, const double e, 
                   const double dx = 0., const double dy = 0., 
                   const double dz = 0.);
  /// Continue the avalanche simulation from the current set of electrons.
  bool ResumeAvalanche();

  /// Set a callback function to be called at every step.
  void SetUserHandleStep(void (*f)(double x, double y, double z, double t,
                                   double e, double dx, double dy, double dz,
                                   bool hole));
  /// Deactivate the user handle called at every step.
  void UnsetUserHandleStep() { m_userHandleStep = nullptr; }
  /// Set a callback function to be called at every (real) collision.
  void SetUserHandleCollision(void (*f)(double x, double y, double z, double t,
                                        int type, int level, Medium* m,
                                        double e0, double e1, double dx0,
                                        double dy0, double dz0, double dx1,
                                        double dy1, double dz1));
  /// Deactivate the user handle called at every collision.
  void UnsetUserHandleCollision() { m_userHandleCollision = nullptr; }
  /// Set a user handling procedure, to be called at every attachment.
  void SetUserHandleAttachment(void (*f)(double x, double y, double z, double t,
                                         int type, int level, Medium* m));
  /// Deactivate the user handle called at every attachment.
  void UnsetUserHandleAttachment() { m_userHandleAttachment = nullptr; }
  /// Set a user handling procedure, to be called at every inelastic collision.
  void SetUserHandleInelastic(void (*f)(double x, double y, double z, double t,
                                        int type, int level, Medium* m));
  /// Deactivate the user handle called at every inelastic collision.
  void UnsetUserHandleInelastic() { m_userHandleInelastic = nullptr; }
  /// Set a user handling procedure, to be called at every ionising collision
  /// or excitation followed by Penning transfer.
  void SetUserHandleIonisation(void (*f)(double x, double y, double z, double t,
                                         int type, int level, Medium* m));
  /// Deactivate the user handle called at every ionisation.
  void UnsetUserHandleIonisation() { m_userHandleIonisation = nullptr; }

  /// Switch on debugging messages.
  void EnableDebugging() { m_debug = true; }
  void DisableDebugging() { m_debug = false; }
  /// Empirical scaling applied only to the first actual multiplication
  /// of the primary carrier.
  void SetFirstPrimaryMultiplicationScale(const double scale) {
    m_firstPrimaryMultiplicationScale =
        scale < 0. ? 0. : (scale > 1. ? 1. : scale);
  }
  double GetFirstPrimaryMultiplicationScale() const {
    return m_firstPrimaryMultiplicationScale;
  }

  // class to store CPU/GPU benchmark comparisons
  struct Statistics {
    std::vector<double> gpu_stack_process_time;
    std::vector<double> cpu_stack_process_time;
    std::vector<double> gpu_stack_transport_time;
    std::vector<double> cpu_stack_transport_time;

    std::vector<size_t> stack_old_size;
    std::vector<size_t> stack_new_size;
  };

  Statistics GetStatistics() { return m_stats; }
  
  void SetRunModeOptions(MPRunMode mode, int device = -1);
  void SetMaxNumShowerLoops(int max_loops) {m_maxNumShowerLoops = max_loops;}
  void SetShowProgress(bool show_progress) {m_showProgress = show_progress;}
  /// Switch on/off per-shower timing/statistics collection.
  void EnableTransportStatistics(const bool on = true) {
    m_collectStatistics = on;
  }
  /// Switch on/off retaining terminal carrier trajectories.
  /// Disable this only when GetAvalancheSize is the required result.
  void EnableEndpointStorage(const bool on = true) { m_storeEndpoints = on; }
  void SetDebugShowerIterationAndElectronID(int iter_num, int elec_id) {m_debugShowerLoopNum = iter_num; m_debugElectronID = elec_id;}

 private:
  std::string m_className = "AvalancheMicroscopic";

  MPRunMode m_runMode{MPRunMode::Normal};
  int m_cudaDevice{-1};
  int m_maxNumShowerLoops{-1};
  int m_debugElectronID{-1};
  int m_debugShowerLoopNum{-1};
  bool m_showProgress{false};
  bool m_collectStatistics{true};
  bool m_storeEndpoints{true};
  AvalancheMicroscopicGPU *m_gpuInterface{nullptr};

  Statistics m_stats;
  
  Sensor* m_sensor = nullptr;

  std::vector<Electron> m_electrons;
  std::vector<Electron> m_electrons_gpu;
  std::vector<Electron> m_holes;

  std::vector<std::pair<Point, Particle> > m_stackStoreCPU;
  std::vector<Electron> m_stackStoreGPU;

  struct photon {
    int status;             ///< Status
    double energy;          ///< Energy
    double x0, y0, z0, t0;  ///< Starting point and time.
    double x1, y1, z1, t1;  ///< End point and time.
  };
  std::vector<photon> m_photons;

  /// Number of electrons produced
  int m_nElectrons = 0;
  /// Number of holes produced
  int m_nHoles = 0;
  /// Number of ions produced
  int m_nIons = 0;

  /// Number of electrons produced
  int m_nElectronsGPU = 0;
  /// Number of holes produced
  int m_nHolesGPU = 0;
  /// Number of ions produced
  int m_nIonsGPU = 0;

  ViewDrift* m_viewer = nullptr;
  bool m_plotExcitations = true;
  bool m_plotIonisations = true;
  bool m_plotAttachments = true;

  TH1* m_histElectronEnergy = nullptr;
  TH1* m_histHoleEnergy = nullptr;
  TH1* m_histDistance = nullptr;
  char m_distanceOption = 'r';
  std::vector<int> m_distanceHistogramType;

  TH1* m_histSecondary = nullptr;

  bool m_doSignal = true;
  bool m_useWeightingPotential = true;
  bool m_integrateWeightingField = false;
  bool m_doInducedCharge = false;

  bool m_computePathLength = false;
  bool m_storeDriftLines = false;
  bool m_usePhotons = false;
  bool m_useBandStructure = true;
  bool m_useNullCollisionSteps = false;
  bool m_useBfieldAuto = true;
  bool m_useBfield = false;
  
  bool m_rknSteps = false;
  double m_rknsteperrortol = 1.e-10;
  double m_rknMinh = 1.e-5;
  double m_nullCollScale = 1.;

  // Transport cuts
  double m_deltaCut = 0.;
  double m_gammaCut = 0.;

  // Max. avalanche size
  unsigned int m_sizeCut = 0;

  size_t m_nCollSkip = 100;
  size_t m_nCollPlot = 100;

  bool m_hasTimeWindow = false;
  double m_tMin = 0.;
  double m_tMax = 0.;

  // User procedures
  void (*m_userHandleStep)(double x, double y, double z, double t, double e,
                           double dx, double dy, double dz,
                           bool hole) = nullptr;
  void (*m_userHandleCollision)(double x, double y, double z, double t,
                                int type, int level, Medium* m, double e0,
                                double e1, double dx0, double dy0, double dz0,
                                double dx1, double dy1, double dz1) = nullptr;
  void (*m_userHandleAttachment)(double x, double y, double z, double t,
                                 int type, int level, Medium* m) = nullptr;
  void (*m_userHandleInelastic)(double x, double y, double z, double t,
                                int type, int level, Medium* m) = nullptr;
  void (*m_userHandleIonisation)(double x, double y, double z, double t,
                                 int type, int level, Medium* m) = nullptr;

  // Switch on/off debugging messages
  bool m_debug = false;
  double m_firstPrimaryMultiplicationScale = 1.;

  bool TransportElectrons(std::vector<std::pair<Point, Particle> >& stack,
                          const bool aval);
  int TransportElectron(const Point& p0, const bool hole, 
                        const bool aval, const bool signal,
                        const bool primaryGeneration,
                        std::vector<double>& ts, 
                        std::vector<std::array<double, 3> >& xs,
                        std::vector<Point>& path,
                        std::vector<std::pair<Point, Particle> >& newParticles);
  int TransportElectronBfield(const Point& p0, const bool hole, 
                              const bool aval, 
                              const bool signal,
                              const bool primaryGeneration, 
                              std::vector<double>& ts, 
                              std::vector<std::array<double, 3> >& xs,
                              std::vector<Point>& path,
                              std::vector<std::pair<Point, Particle> >& newParticles);
  int TransportElectronSc(const Point& p0, const bool hole, 
                          const bool aval, 
                          const bool signal,
                          const bool primaryGeneration, 
                          std::vector<double>& ts, 
                          std::vector<std::array<double, 3> >& xs,
                          std::vector<Point>& path,
                          std::vector<std::pair<Point, Particle> >& newParticles);
  void TransportPhoton(const double x, const double y, const double z,
                       const double t, const double e,
                       std::vector<std::pair<Point, Particle> >& newParticles);

  bool transportParticleStack(const bool aval,
                              std::vector<std::pair<Point, Particle> > &particles, 
                              std::vector<std::pair<Point, Particle> > &newParticles,
                              const bool primaryGeneration,
                              const bool signal,
                              const bool useBfield,
                              const bool sc);

  static bool IsInactive(const Electron& item) {
    return item.status == StatusLeftDriftMedium ||
           item.status == StatusBelowTransportCut ||
           item.status == StatusOutsideTimeWindow ||
           item.status == StatusLeftDriftArea || 
           item.status == StatusAttached ||
           item.status == StatusHitPlane;
  }
  void Update(std::vector<Electron>::iterator it, const double x,
              const double y, const double z, const double t,
              const double energy, const double kx, const double ky,
              const double kz, const int band);
  void AddToEndPoints(const Electron& item, const bool hole) {
    Electron electron;
    electron.status = item.status;
    electron.path = item.path;
    if (hole) {
      m_holes.push_back(std::move(electron));
    } else {
      m_electrons.push_back(std::move(electron));
    }
  }

  void Terminate(double x0, double y0, double z0, double t0, double& x1,
                 double& y1, double& z1, double& t1) const;

  void PlotCollision(const int cstype, const size_t did,
                     const double x, const double y, const double z,
                     size_t& nCollPlot) const;
  void FillDistanceHistogram(const int cstype,
                             const double x, const double y, const double z,
                             double& xLast, double& yLast, double& zLast) const;

public:
  std::vector<std::pair<Point, Particle> > GetStackOld() { return m_stackStoreCPU; }
  std::vector<Electron> GetStackOldGPU() {return m_stackStoreGPU; }

friend class AvalancheMicroscopicGPU;
};
}

#endif
