#ifndef G_AVALANCHE_GRID_H
#define G_AVALANCHE_GRID_H

#include <iostream>
#include <string>
#include <vector>

#include "AvalancheMicroscopic.hh"
#include "Garfield/ComponentParallelPlate.hh"
#include "GarfieldConstants.hh"
#include "Sensor.hh"

namespace Garfield {

/// Calculate avalanches in a uniform electric field using avalanche
/// statistics.

class AvalancheGrid {
 public:
  /// Constructor
  AvalancheGrid() {}
  /// Destructor
  ~AvalancheGrid() {}
  /// Set the sensor.
  void SetSensor(Sensor *sensor) { m_sensor = sensor; }

  /** Start grid based avalanche simulation.
   *
   * \param zmin,zmax z-coordinate range of grid [cm].
   * \param zsteps amount of z-coordinate points in grid.
   * \param xmin,xmax x-coordinate range of grid [cm].
   * \param xsteps amount of x-coordinate points in grid.
   */
  void StartGridAvalanche();
  /// Set the electron drift velocity (in cm / ns).
  void SetElectronVelocity(const double vx, const double vy, const double vz) {
    double vel = sqrt(vx * vx + vy * vy + vz * vz);
    if (vel != std::abs(vx) && vel != std::abs(vy) && vel != std::abs(vz))
      return;
    int nx = (int)vx / vel;
    int ny = (int)vy / vel;
    int nz = (int)vz / vel;
    m_velNormal = {nx, ny, nz};
    m_Velocity = -std::abs(vel);
  }
  /// Set the electron Townsend coefficient (in 1 / cm).
  void SetElectronTownsend(const double town) { m_Townsend = town; }
  /// Set the electron attachment coefficient (in 1 / cm).
  void SetElectronAttachment(const double att) { m_Attachment = att; }
  /// Set the maximum avalanche size (1e7 by default).
  void SetMaxAvalancheSize(const double size) { m_MaxSize = size; }
  /// Enable transverse diffusion of electrons with transverse diffusion
  /// coefficients (in √cm).
  void EnableDiffusion(const double diffSigma) {
    m_diffusion = true;
    m_DiffSigma = diffSigma;
  }
  /** Add an electron to the initial configuration.
   *
   * \param x x-coordinate of initial electron.
   * \param y y-coordinate of initial electron.
   * \param z t-coordinate of initial electron.
   * \param t starting time of avalanche.
   * \param n number of electrons at this point.
   */
  void AvalancheElectron(const double x, const double y, const double z,
                         const double t = 0, const int n = 1);
  /// Import electron data from AvalancheMicroscopic class
  void ImportElectronsFromAvalancheMicroscopic(AvalancheMicroscopic *avmc);

  /// Import electron data from AvalancheMicroscopic class
  void SetGrid(const double xmin, const double xmax, const int xsteps,
               const double ymin, const double ymax, const int ysteps,
               const double zmin, const double zmax, const int zsteps);

  /// Returns the initial number of electrons in the avalanche.
  int GetAmountOfStartingElectrons() { return m_nestart; }
  /// Returns the final number of electrons in the avalanche.
  int GetAvalancheSize() { return m_nTotal; }

  /// Asigning layer index to all Avalanche nodes.
  void AsignLayerIndex(ComponentParallelPlate *RPC);

  void EnableDebugging() { m_debug = true; }

  void Reset();

 private:
  bool m_debug = false;

  double m_Townsend = -1;  // [1/cm]

  double m_Attachment = -1;  // [1/cm]

  double m_Velocity = 0.;  // [cm/ns]

  std::vector<int> m_velNormal = {0, 0, 0};

  double m_MaxSize = 1.6e7;  // Saturations size
  // Check if avalanche has reached maximum size
  bool m_Saturated = false;  
  // Time when the avalanche has reached maximum size
  double m_SaturationTime = -1.;  


  bool m_diffusion = false;  // Check if transverse diffusion is enabled.

  double m_DiffSigma = 0.;  // Transverse diffusion coefficients (in √cm).

  int m_nestart = 0.;

  bool m_driftAvalanche = false;
  bool m_importAvalanche = false;

  bool m_layerIndix = false;
  std::vector<double> m_nLayer;

  std::string m_className = "AvalancheGrid";

  Sensor *m_sensor = nullptr;

  bool m_printPar = false;

  std::vector<double> m_zgrid;  ///< Grid points of z-coordinate.
  int m_zsteps = 0;             ///< Number of grid points.
  double m_zStepSize = 0.;      ///< Distance between the grid points.

  std::vector<double> m_ygrid;  ///< Grid points of y-coordinate.
  int m_ysteps = 0;             ///< Number of grid points.
  double m_yStepSize = 0.;      ///< Distance between the grid points.

  std::vector<double> m_xgrid;  ///< Grid points of x-coordinate.
  int m_xsteps = 0;             ///< Number of grid points.
  double m_xStepSize = 0.;      ///< Distance between the grid points.

  bool m_gridset = false; ///< Keeps track if the grid has been defined.
  int m_nTotal = 0;       ///< Total amount of charge.
  double m_time = 0;      ///< Clock.
  bool m_run = true;      ///< Tracking if the charges are still in the drift gap.
 
 struct Path {
    std::vector<double> ts ={};
    std::vector< std::array<double, 3> > xs ={};
    std::vector<double> qs = {};
 };

  struct AvalancheNode {
    double ix = 0;
    double iy = 0;
    double iz = 0;

    int n = 1;

    int layer = 0;

    double townsend = 0;
    double attachment = 0;
    double velocity = 0;

    double stepSize = 0;
    std::vector<int> velNormal = {0, 0, 0};

    double time = 0.;  ///< Clock.
    double dt = -1.;   ///< time step.

    bool active = true;
    double dSigmaL = 0;
    double dSigmaT = 0;
      
    Path path;
  };

  std::vector<AvalancheNode> m_activeNodes = {};

  // Assign electron to the closest grid point.
  bool SnapToGrid(const double x, const double y, const double z,
                  const double v, const int n = 1);
  // Go to next time step.
  void NextAvalancheGridPoint();
  // Obtain the Townsend coef., Attachment coef. and velocity vector from
  // sensor class.
  bool GetParameters(AvalancheNode &node);

  void DeactivateNode(AvalancheNode &node);
};
}  // namespace Garfield

#endif
