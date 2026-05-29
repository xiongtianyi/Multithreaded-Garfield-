#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include <iomanip>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"

#define __GPUCOMPILE__
#include "AvalancheMicroscopicGPU.h"
#undef __GPUCOMPILE__

using highres_clock_t = std::chrono::high_resolution_clock;
using second_t = std::chrono::duration<double, std::ratio<1> >;

namespace {

double Mag(const double x, const double y, const double z) {
  return sqrt(x * x + y * y + z * z);
}

void Normalise(double& x, double& y, double& z) {
  const double d = Mag(x, y, z);
  if (d > 0.) {
    const double scale = 1. / d;
    x *= scale;
    y *= scale;
    z *= scale;
  }
}

void ToLocal(const std::array<std::array<double, 3>, 3>& rot, const double xg,
             const double yg, const double zg, double& xl, double& yl,
             double& zl) {
  xl = rot[0][0] * xg + rot[0][1] * yg + rot[0][2] * zg;
  yl = rot[1][0] * xg + rot[1][1] * yg + rot[1][2] * zg;
  zl = rot[2][0] * xg + rot[2][1] * yg + rot[2][2] * zg;
}

void ToGlobal(const std::array<std::array<double, 3>, 3>& rot, const double xl,
              const double yl, const double zl, double& xg, double& yg,
              double& zg) {
  xg = rot[0][0] * xl + rot[1][0] * yl + rot[2][0] * zl;
  yg = rot[0][1] * xl + rot[1][1] * yl + rot[2][1] * zl;
  zg = rot[0][2] * xl + rot[1][2] * yl + rot[2][2] * zl;
}

void RotationMatrix(double bx, double by, double bz, const double bmag,
                    const double ex, const double ey, const double ez,
                    std::array<std::array<double, 3>, 3>& rot) {
  // Adopting the Magboltz convention, the stepping is performed
  // in a coordinate system with the B field along the x axis
  // and the electric field at an angle btheta in the x-z plane.

  // Calculate the first rotation matrix (to align B with x axis).
  std::array<std::array<double, 3>, 3> rB = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  if (bmag > Garfield::Small) {
    bx /= bmag;
    by /= bmag;
    bz /= bmag;
    const double bt = by * by + bz * bz;
    if (bt > Garfield::Small) {
      const double btInv = 1. / bt;
      rB[0][0] = bx;
      rB[0][1] = by;
      rB[0][2] = bz;
      rB[1][0] = -by;
      rB[2][0] = -bz;
      rB[1][1] = (bx * by * by + bz * bz) * btInv;
      rB[2][2] = (bx * bz * bz + by * by) * btInv;
      rB[1][2] = rB[2][1] = (bx - 1.) * by * bz * btInv;
    } else if (bx < 0.) {
      // B field is anti-parallel to x.
      rB[0][0] = -1.;
      rB[1][1] = -1.;
    }
  }
  // Calculate the second rotation matrix (rotation around x axis).
  const double fy = rB[1][0] * ex + rB[1][1] * ey + rB[1][2] * ez;
  const double fz = rB[2][0] * ex + rB[2][1] * ey + rB[2][2] * ez;
  const double ft = sqrt(fy * fy + fz * fz);
  std::array<std::array<double, 3>, 3> rX = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  if (ft > Garfield::Small) {
    // E and B field are not parallel.
    rX[1][1] = rX[2][2] = fz / ft;
    rX[1][2] = -fy / ft;
    rX[2][1] = -rX[1][2];
  }
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      rot[i][j] = 0.;
      for (unsigned int k = 0; k < 3; ++k) {
        rot[i][j] += rX[i][k] * rB[k][j];
      }
    }
  }
}

Garfield::AvalancheMicroscopic::Point MakePoint(
    const double x, const double y, const double z,  const double t, 
    const double energy, const double dx, const double dy, const double dz,
    const int band) {

  Garfield::AvalancheMicroscopic::Point p;
  p.x = x;
  p.y = y;
  p.z = z;
  p.t = t;
  p.energy = energy;
  p.kx = dx;
  p.ky = dy;
  p.kz = dz;
  p.band = band;
  return p;
}

Garfield::AvalancheMicroscopic::Point MakePoint(
    const double x, const double y, const double z,  const double t, 
    const double energy) {
  // Randomise the direction.
  double dx = 0., dy = 0., dz = 1.;
  Garfield::RndmDirection(dx, dy, dz);
  return MakePoint(x, y, z, t, energy, dx, dy, dz, 0);
}

void PrintStatus(const std::string& hdr, const std::string& status,
                 const double x, const double y, const double z,
                 const bool hole) {
  const std::string eh = hole ? "Hole " : "Electron ";
  std::cout << hdr << eh << status << " at " << x << ", " << y << ", " << z
            << "\n";
}

}  // namespace

namespace Garfield {

AvalancheMicroscopic::AvalancheMicroscopic(Sensor* sensor) :
    m_sensor(sensor) {
  m_electrons.reserve(10000);
  m_holes.reserve(10000);
  m_photons.reserve(1000);
}

void AvalancheMicroscopic::SetSensor(Sensor* s) {
  if (!s) {
    std::cerr << m_className << "::SetSensor: Null pointer.\n";
    return;
  }
  m_sensor = s;
}

void AvalancheMicroscopic::EnablePlotting(ViewDrift* view,
                                          const size_t nColl) {
  if (!view) {
    std::cerr << m_className << "::EnablePlotting: Null pointer.\n";
    return;
  }
  m_viewer = view;
  m_nCollPlot = std::max(nColl, 1ul);
}

void AvalancheMicroscopic::EnableElectronEnergyHistogramming(TH1* histo) {
  if (!histo) {
    std::cerr << m_className << "::EnableElectronEnergyHistogramming:\n"
              << "    Null pointer.\n";
    return;
  }

  m_histElectronEnergy = histo;
}

void AvalancheMicroscopic::EnableHoleEnergyHistogramming(TH1* histo) {
  if (!histo) {
    std::cerr << m_className << "::EnableHoleEnergyHistogramming:\n"
              << "    Null pointer.\n";
    return;
  }

  m_histHoleEnergy = histo;
}

void AvalancheMicroscopic::SetDistanceHistogram(TH1* histo, const char opt) {
  if (!histo) {
    std::cerr << m_className << "::SetDistanceHistogram: Null pointer.\n";
    return;
  }

  m_histDistance = histo;

  if (opt == 'x' || opt == 'y' || opt == 'z' || opt == 'r') {
    m_distanceOption = opt;
  } else {
    std::cerr << m_className << "::SetDistanceHistogram:";
    std::cerr << "    Unknown option " << opt << ".\n";
    std::cerr << "    Valid options are x, y, z, r.\n";
    std::cerr << "    Using default value (r).\n";
    m_distanceOption = 'r';
  }

  if (m_distanceHistogramType.empty()) {
    std::cout << m_className << "::SetDistanceHistogram:\n";
    std::cout << "    Don't forget to call EnableDistanceHistogramming.\n";
  }
}

void AvalancheMicroscopic::EnableDistanceHistogramming(const int type) {
  // Check if this type of collision is already registered
  // for histogramming.
  const unsigned int nDistanceHistogramTypes = m_distanceHistogramType.size();
  if (nDistanceHistogramTypes > 0) {
    for (unsigned int i = 0; i < nDistanceHistogramTypes; ++i) {
      if (m_distanceHistogramType[i] != type) continue;
      std::cout << m_className << "::EnableDistanceHistogramming:\n";
      std::cout << "    Collision type " << type
                << " is already being histogrammed.\n";
      return;
    }
  }

  m_distanceHistogramType.push_back(type);
  std::cout << m_className << "::EnableDistanceHistogramming:\n";
  std::cout << "    Histogramming of collision type " << type << " enabled.\n";
  if (!m_histDistance) {
    std::cout << "    Don't forget to set the histogram.\n";
  }
}

void AvalancheMicroscopic::DisableDistanceHistogramming(const int type) {
  if (std::find(m_distanceHistogramType.begin(), m_distanceHistogramType.end(),
                type) == m_distanceHistogramType.end()) {
    std::cerr << m_className << "::DisableDistanceHistogramming:\n"
              << "    Collision type " << type << " is not histogrammed.\n";
    return;
  }

  m_distanceHistogramType.erase(
      std::remove(m_distanceHistogramType.begin(),
                  m_distanceHistogramType.end(), type),
      m_distanceHistogramType.end());
}

void AvalancheMicroscopic::DisableDistanceHistogramming() {
  m_histDistance = nullptr;
  m_distanceHistogramType.clear();
}

void AvalancheMicroscopic::EnableSecondaryEnergyHistogramming(TH1* histo) {
  if (!histo) {
    std::cerr << m_className << "::EnableSecondaryEnergyHistogramming:\n"
              << "    Null pointer.\n";
    return;
  }

  m_histSecondary = histo;
}

void AvalancheMicroscopic::SetTimeWindow(const double t0, const double t1) {
  if (fabs(t1 - t0) < Small) {
    std::cerr << m_className << "::SetTimeWindow:\n";
    std::cerr << "    Time interval must be greater than zero.\n";
    return;
  }

  m_tMin = std::min(t0, t1);
  m_tMax = std::max(t0, t1);
  m_hasTimeWindow = true;
}

void AvalancheMicroscopic::GetElectronEndpoint(const size_t i, 
    double& x0, double& y0, double& z0, double& t0, double& e0,
    double& x1, double& y1, double& z1, double& t1, double& e1, 
    int& status) const {
  if (i >= m_electrons.size()) {
    std::cerr << m_className << "::GetElectronEndpoint: Index out of range.\n";
    status = -3;
    return;
  }
  if (m_electrons[i].path.empty()) {
    std::cerr << m_className << "::GetElectronEndpoint: Empty drift line.\n";
    status = -3;
    return;
  }
  x0 = m_electrons[i].path[0].x;
  y0 = m_electrons[i].path[0].y;
  z0 = m_electrons[i].path[0].z;
  t0 = m_electrons[i].path[0].t;
  e0 = m_electrons[i].path[0].energy;
  x1 = m_electrons[i].path.back().x;
  y1 = m_electrons[i].path.back().y;
  z1 = m_electrons[i].path.back().z;
  t1 = m_electrons[i].path.back().t;
  e1 = m_electrons[i].path.back().energy;
  status = m_electrons[i].status;
}

void AvalancheMicroscopic::GetElectronEndpointGPU(const size_t i, double& x0,
                                               double& y0, double& z0,
                                               double& t0, double& e0,
                                               double& x1, double& y1,
                                               double& z1, double& t1,
                                               double& e1, int& status) const {
  if (i >= m_electrons_gpu.size()) {
    std::cerr << m_className << "::GetElectronEndpoint: Index out of range.\n";
    status = -3;
    return;
  }

  if (m_electrons_gpu[i].path.empty()) {
    std::cerr << m_className << "::GetElectronEndpoint: Empty drift line.\n";
    status = -3;
  }

  x0 = m_electrons_gpu[i].path[0].x;
  y0 = m_electrons_gpu[i].path[0].y;
  z0 = m_electrons_gpu[i].path[0].z;
  t0 = m_electrons_gpu[i].path[0].t;
  e0 = m_electrons_gpu[i].path[0].energy;
  x1 = m_electrons_gpu[i].path.back().x;
  y1 = m_electrons_gpu[i].path.back().y;
  z1 = m_electrons_gpu[i].path.back().z;
  t1 = m_electrons_gpu[i].path.back().t;
  e1 = m_electrons_gpu[i].path.back().energy;
  status = m_electrons_gpu[i].status;
}

void AvalancheMicroscopic::GetHoleEndpoint(const size_t i, 
    double& x0, double& y0, double& z0, double& t0, double& e0, 
    double& x1, double& y1, double& z1, double& t1, double& e1,
    int& status) const {
  if (i >= m_holes.size()) {
    std::cerr << m_className << "::GetHoleEndpoint: Index out of range.\n";
    status = -3;
    return;
  }
  if (m_electrons[i].path.empty()) {
    std::cerr << m_className << "::GetHoleEndpoint: Empty drift line.\n";
    status = -3;
    return;
  }
  x0 = m_holes[i].path[0].x;
  y0 = m_holes[i].path[0].y;
  z0 = m_holes[i].path[0].z;
  t0 = m_holes[i].path[0].t;
  e0 = m_holes[i].path[0].energy;
  x1 = m_holes[i].path.back().x;
  y1 = m_holes[i].path.back().y;
  z1 = m_holes[i].path.back().z;
  t1 = m_holes[i].path.back().t;
  e1 = m_holes[i].path.back().energy;
  status = m_holes[i].status;
}

size_t AvalancheMicroscopic::GetNumberOfElectronDriftLinePoints(
    const size_t i) const {
  if (i >= m_electrons.size()) {
    std::cerr << m_className << "::GetNumberOfElectronDriftLinePoints: "
              << "Index out of range.\n";
    return 0;
  }
  return m_electrons[i].path.size();
}

size_t AvalancheMicroscopic::GetNumberOfHoleDriftLinePoints(
    const size_t i) const {
  if (i >= m_holes.size()) {
    std::cerr << m_className << "::GetNumberOfHoleDriftLinePoints: "
              << "Index out of range.\n";
    return 0;
  }
  return m_holes[i].path.size();
}

void AvalancheMicroscopic::GetElectronDriftLinePoint(
    double& x, double& y, double& z, double& t, const size_t ip,
    const size_t ie) const {
  if (ie >= m_electrons.size()) {
    std::cerr << m_className << "::GetElectronDriftLinePoint:\n"
              << "    Endpoint index (" << ie << ") out of range.\n";
    return;
  }
  if (ip >= m_electrons[ie].path.size()) {
    std::cerr << m_className << "::GetElectronDriftLinePoint:\n"
              << "    Drift line point index (" << ip << ") out of range.\n";
    return;
  }
  x = m_electrons[ie].path[ip].x;
  y = m_electrons[ie].path[ip].y;
  z = m_electrons[ie].path[ip].z;
  t = m_electrons[ie].path[ip].t;
}

void AvalancheMicroscopic::GetHoleDriftLinePoint(double& x, double& y,
                                                 double& z, double& t,
                                                 const size_t ip,
                                                 const size_t ih) const {
  if (ih >= m_holes.size()) {
    std::cerr << m_className << "::GetHoleDriftLinePoint:\n"
              << "    Endpoint index (" << ih << ") out of range.\n";
    return;
  }
  if (ip >= m_holes[ih].path.size()) {
    std::cerr << m_className << "::GetHoleDriftLinePoint:\n"
              << "    Drift line point index (" << ip << ") out of range.\n";
    return;
  }
  x = m_holes[ih].path[ip].x;
  y = m_holes[ih].path[ip].y;
  z = m_holes[ih].path[ip].z;
  t = m_holes[ih].path[ip].t;
}

void AvalancheMicroscopic::GetPhoton(const size_t i, double& e,
    double& x0, double& y0, double& z0, double& t0, 
    double& x1, double& y1, double& z1, double& t1, int& status) const {
  if (i >= m_photons.size()) {
    std::cerr << m_className << "::GetPhoton: Index out of range.\n";
    return;
  }

  x0 = m_photons[i].x0;
  x1 = m_photons[i].x1;
  y0 = m_photons[i].y0;
  y1 = m_photons[i].y1;
  z0 = m_photons[i].z0;
  z1 = m_photons[i].z1;
  t0 = m_photons[i].t0;
  t1 = m_photons[i].t1;
  status = m_photons[i].status;
  e = m_photons[i].energy;
}

void AvalancheMicroscopic::SetUserHandleStep(
    void (*f)(double x, double y, double z, double t, double e, double dx,
              double dy, double dz, bool hole)) {
  if (!f) {
    std::cerr << m_className << "::SetUserHandleStep: Null pointer.\n";
    return;
  }
  m_userHandleStep = f;
}

void AvalancheMicroscopic::SetUserHandleCollision(
    void (*f)(double x, double y, double z, double t, int type, int level,
              Medium* m, double e0, double e1, double dx0, double dy0,
              double dz0, double dx1, double dy1, double dz1)) {
  m_userHandleCollision = f;
}

void AvalancheMicroscopic::SetUserHandleAttachment(void (*f)(
    double x, double y, double z, double t, int type, int level, Medium* m)) {
  m_userHandleAttachment = f;
}

void AvalancheMicroscopic::SetUserHandleInelastic(void (*f)(
    double x, double y, double z, double t, int type, int level, Medium* m)) {
  m_userHandleInelastic = f;
}

void AvalancheMicroscopic::SetUserHandleIonisation(void (*f)(
    double x, double y, double z, double t, int type, int level, Medium* m)) {
  m_userHandleIonisation = f;
}

bool AvalancheMicroscopic::DriftElectron(
    const double x, const double y, const double z, const double t,
    const double e, const double dx, const double dy, const double dz) {
  std::vector<std::pair<Point, Particle> > particles;
  Point p = MakePoint(x, y, z, t, e, dx, dy, dz, 0);
  particles.emplace_back(std::make_pair(std::move(p), Particle::Electron));
  return TransportElectrons(particles, false);
}

bool AvalancheMicroscopic::AvalancheElectron(
    const double x, const double y, const double z, const double t,
    const double e, const double dx, const double dy, const double dz) {

  std::vector<std::pair<Point, Particle> > particles;
  Point p = MakePoint(x, y, z, t, e, dx, dy, dz, 0);
  particles.emplace_back(std::make_pair(std::move(p), Particle::Electron));
  return TransportElectrons(particles, true);
}

void AvalancheMicroscopic::AddElectron(
    const double x, const double y, const double z, const double t,
    const double e, const double dx, const double dy, const double dz) {

  Electron electron;
  electron.status = StatusAlive;
  electron.path.emplace_back(MakePoint(x, y, z, t, e, dx, dy, dz, 0));
  m_electrons.push_back(std::move(electron));
}

bool AvalancheMicroscopic::ResumeAvalanche() {
  std::vector<std::pair<Point, Particle> > particles;
  for (const auto& p : m_electrons) {
    if (p.status == StatusAlive || p.status == StatusOutsideTimeWindow) { 
      particles.emplace_back(std::make_pair(p.path.back(), Particle::Electron));
    }
  }
  for (const auto& p : m_holes) {
    if (p.status == StatusAlive || p.status == StatusOutsideTimeWindow) { 
      particles.emplace_back(std::make_pair(p.path.back(), Particle::Hole));
    }
  }
  return TransportElectrons(particles, true);
}

bool AvalancheMicroscopic::TransportElectrons(
    std::vector<std::pair<Point, Particle> >& particles, const bool aval) {

  // Clear the list of electrons, holes and photons.
  m_electrons.clear();
  m_holes.clear();
  m_photons.clear();

  // Reset the particle counters.
  m_nElectrons = m_nHoles = m_nIons = 0;

  // Make sure that the sensor is defined.
  if (!m_sensor) {
    std::cerr << m_className 
              << "::TransportElectrons: Sensor is not defined.\n";
    return false;
  }

  // Do we need to consider the magnetic field?
  const bool useBfield = m_useBfieldAuto ? m_sensor->HasMagneticField() : 
                         m_useBfield;

  // Do we need to compute the induced signal?
  const bool signal = m_doSignal && (m_sensor->GetNumberOfElectrodes() > 0);
  bool sc = false;
  // Loop over the initial set of electrons/holes.
  for (auto& p : particles) {
    // Make sure that the starting point is inside the active area.
    const double x0 = p.first.x;
    const double y0 = p.first.y;
    const double z0 = p.first.z;
    if (!m_sensor->IsInArea(x0, y0, z0)) {
      std::cerr << m_className << "::TransportElectrons: "
                << "Starting point is outside the active area.\n";
      return false;
    }
    // Make sure that the starting point is inside a "driftable" 
    // microscopic medium.
    Medium* medium = m_sensor->GetMedium(x0, y0, z0);
    if (!medium || !medium->IsDriftable() || !medium->IsMicroscopic()) {
      std::cerr << m_className << "::TransportElectrons: "
                << "Starting point is not in a valid medium.\n";
      return false;
    } 
    // Make sure the initial energy is positive.
    const double e0 = std::max(p.first.energy, Small);

    if (medium->IsSemiconductor() && m_useBandStructure) {
      sc = true;
      if (p.first.band < 0) {
        // Sample the initial momentum and band.
        medium->GetElectronMomentum(e0, p.first.kx, p.first.ky, p.first.kz,
                                    p.first.band);
      }
    } else {
      p.first.band = 0;
      const double kmag = Mag(p.first.kx, p.first.ky, p.first.kz);
      if (fabs(kmag) < Small) {
        // Direction has zero norm, draw a random direction.
        RndmDirection(p.first.kx, p.first.ky, p.first.kz);
      } else {
        // Normalise the direction to 1.
        const double scale = 1. / kmag;
        p.first.kx *= scale;
        p.first.ky *= scale;
        p.first.kz *= scale;
      }
    }
  }
  std::vector<std::pair<Point, Particle> > newParticles;
  m_stats.gpu_stack_process_time.clear();
  m_stats.cpu_stack_process_time.clear();
  m_stats.gpu_stack_transport_time.clear();
  m_stats.cpu_stack_transport_time.clear();
  m_stats.stack_old_size.clear();
  m_stats.stack_new_size.clear();

  // GPUREMOVE: Note that the change below means that the transferred variables (medium, id, etc.) are not updated within the transport
  // loops so will have to be recalculated more than normal. This will affect both CPU and GPU versions.
  int loop_count = 0;

  if (m_runMode == MPRunMode::GPUExclusive) {
#ifdef USEGPU
    if (!m_gpuInterface) {
      m_gpuInterface = new AvalancheMicroscopicGPU;
      m_gpuInterface->SetCUDADevice(m_cudaDevice);
      m_gpuInterface->TransferClassInternalInfo(this);
    }
    m_gpuInterface->TransferStackFromCPUToGPU(particles);
#else
    std::cout << "ERROR: GPU use requested but Garfield has not been built with GPU support" << std::endl;
    return false;
#endif
  }

  std::chrono::time_point< highres_clock_t > start;
  double process_time_cpu{0};
  double stack_time_cpu{0};
  double process_time_gpu{0};
  double stack_time_gpu{0};
  unsigned int num_new_particles{0};
  unsigned int num_curr_particles{0};
  unsigned int num_new_particles_gpu{0};
  unsigned int num_curr_particles_gpu{0};
  bool primaryGeneration = true;

  while ((m_maxNumShowerLoops == -1) || (loop_count < m_maxNumShowerLoops)) {
    if (m_showProgress) {
      std::cout << "--------------------------------------------------------------" << std::endl;
      std::cout << "Starting Shower iteration:         " << loop_count << std::endl;
    }

    // --------------------------------------------
    // Process and transport the particle stack depending on GPU config
    if (m_runMode == MPRunMode::Normal) {
      if (m_collectStatistics) start = highres_clock_t::now();
      num_new_particles = newParticles.size();

      // this is all the processing of the particle stack that's needed
      if (loop_count) {
        if (!aval) break;
        particles.swap(newParticles);
      }

      if (particles.size() == 0)
        break;

      if (m_collectStatistics) {
        process_time_cpu = std::chrono::duration_cast<second_t>(
            highres_clock_t::now() - start).count();
        start = highres_clock_t::now();
      }

      num_curr_particles = particles.size();

      if (!transportParticleStack(aval, particles, newParticles,
                                  primaryGeneration, signal, useBfield, sc))
        return false;
      primaryGeneration = false;

      if (m_collectStatistics) {
        stack_time_cpu = std::chrono::duration_cast<second_t>(
            highres_clock_t::now() - start).count();
      }
    } 
    else if (m_runMode == MPRunMode::GPUExclusive) {
#ifdef USEGPU
      if (m_collectStatistics) start = highres_clock_t::now();

      if (m_gpuInterface->processParticleStack(num_curr_particles_gpu, num_new_particles_gpu) == 0) {
        break;
      }

      if (m_collectStatistics) {
        process_time_gpu = std::chrono::duration_cast<second_t>(
            highres_clock_t::now() - start).count();
        start = highres_clock_t::now();
      }

      // TODO: TN GPU: Fix arguments (medium, ID, useBandStructure, Flim, Finv all set
      // to constant values)
      if (!m_gpuInterface->transportParticleStack(aval, this, 0, false, 0, 0, 0, 0, useBfield, sc))
        return false;

      if (m_collectStatistics) {
        stack_time_gpu = std::chrono::duration_cast<second_t>(
            highres_clock_t::now() - start).count();
      }
#endif
    }

    loop_count++;    

    if (m_showProgress) {
      if (m_runMode == MPRunMode::Normal)
        std::cout << "    - Current particle stack size (CPU): " << num_curr_particles << std::endl;
      if (m_runMode == MPRunMode::GPUExclusive)
        std::cout << "    - Current particle stack size (GPU): " << num_curr_particles_gpu << std::endl;

      if (m_stats.cpu_stack_transport_time.size() > 0) {
        std::cout << "CPU Stats: Stk (" << *(m_stats.cpu_stack_process_time.end()-1) << 
                  "),  Tpt (" << *(m_stats.cpu_stack_transport_time.end()-1) << ")" << std::endl;
      }
      if (m_stats.gpu_stack_transport_time.size() > 0) {
        std::cout << "GPU Stats: Stk (" << *(m_stats.gpu_stack_process_time.end()-1) << 
                  "),  Tpt (" << *(m_stats.gpu_stack_transport_time.end()-1) << ")" << std::endl;
      }
      std::cout << "--------------------------------------------------------------" << std::endl;
    }

    // Statistics are diagnostic output and are not needed in rate benchmarks.
    if (m_collectStatistics) {
      if (m_runMode == MPRunMode::GPUExclusive) {
        m_stats.stack_old_size.push_back(num_curr_particles_gpu);
        m_stats.stack_new_size.push_back(num_new_particles_gpu);
      } else {
        m_stats.stack_old_size.push_back(num_curr_particles);
        m_stats.stack_new_size.push_back(num_new_particles_gpu);
      }
      
      if (stack_time_cpu > 0) {
        m_stats.cpu_stack_transport_time.push_back(stack_time_cpu);
        m_stats.cpu_stack_process_time.push_back(process_time_cpu);
      }

      if (stack_time_gpu > 0) {
        m_stats.gpu_stack_transport_time.push_back(stack_time_gpu);
        m_stats.gpu_stack_process_time.push_back(process_time_gpu);
      }
    }
  }

  // Multiprocessor clean up
  if (m_runMode == MPRunMode::GPUExclusive) {
    #ifdef USEGPU
    // copy over stack if there's any to compare
    if (m_maxNumShowerLoops > -1) {
      m_gpuInterface->TransferStackFromGPUToCPU(m_stackStoreGPU, false);  
      m_stackStoreCPU = particles;
    }

    // Copy endpoints over
    m_gpuInterface->TransferStackFromGPUToCPU(m_electrons_gpu, true);
    m_sensor->TransferGPUElectrodeSignals(m_gpuInterface->m_sensor);
    #endif
  }

  // Calculate the induced charge.
  if (m_doInducedCharge) {
    for (const auto& p : m_electrons) {
      m_sensor->AddInducedCharge(-1, p.path[0].x, p.path[0].y, p.path[0].z, 
                                 p.path.back().x, p.path.back().y, p.path.back().z);
    }
    for (const auto& p : m_holes) {
      m_sensor->AddInducedCharge(+1, p.path[0].x, p.path[0].y, p.path[0].z, 
                                 p.path.back().x, p.path.back().y, p.path.back().z);
    }
  }
  return true;
}

bool AvalancheMicroscopic::transportParticleStack(
  const bool aval,
  std::vector<std::pair<Point, Particle> > &particles,
  std::vector<std::pair<Point, Particle> > &newParticles,
  const bool primaryGeneration,
  const bool signal,
  const bool useBfield,
  const bool sc) {
  newParticles.clear();
  // Loop over the particles in the avalanche.
  for (const auto& particle : particles) {
    if (particle.second == Particle::Ion) {
      ++m_nIons;
      continue;
    }
    if (aval && m_sizeCut > 0 && m_nElectrons >= (int)m_sizeCut) { 
        newParticles.clear();
        break;
    }
    const bool isHole = (particle.second == Particle::Hole);
    std::vector<Point> path;
    std::vector<double> ts;
    std::vector<std::array<double, 3> > xs;
    int status = 0;
    if (sc) {
      status = TransportElectronSc(particle.first, isHole, aval,
                                    signal, primaryGeneration, ts, xs, path,
                                    newParticles);
    } else if (useBfield) {
      status = TransportElectronBfield(particle.first, isHole, aval,
                                        signal, primaryGeneration, ts, xs, path,
                                        newParticles);
    } else {
      status = TransportElectron(particle.first, isHole, aval, signal,
                                  primaryGeneration,
                                  ts, xs, path, newParticles);
    }
    double pathLength = 0.;
    if (m_computePathLength && xs.size() > 1) {
      const size_t ns = xs.size();
      for (size_t i = 0; i < ns - 1; ++i) {
        pathLength += Mag(xs[i + 1][0] - xs[i][0], 
                          xs[i + 1][1] - xs[i][1],
                          xs[i + 1][2] - xs[i][2]);
      }
    }
    if (isHole) {
      if (m_storeEndpoints || m_doInducedCharge) {
        Electron hole;
        hole.status = status;
        hole.path = std::move(path);
        hole.pathLength = pathLength;
        m_holes.push_back(std::move(hole));
      }
      if (status != StatusAttached) ++m_nHoles;
    } else {
      if (m_storeEndpoints || m_doInducedCharge) {
        Electron electron;
        electron.status = status;
        electron.path = std::move(path);
        electron.pathLength = pathLength;
        m_electrons.push_back(std::move(electron));
      }
      if (status != StatusAttached) ++m_nElectrons;
    }
    if (signal) {
      const double q = isHole ? 1. : -1.;
      if (m_useWeightingPotential) {
        m_sensor->AddSignalWeightingPotential(q, ts, xs);
      } else {
        m_sensor->AddSignalWeightingField(q, ts, xs, 
                                          m_integrateWeightingField);
      }
    }
  }
  return true;
}

int AvalancheMicroscopic::TransportElectron(const Point& p0,
  const bool hole, const bool aval, const bool signal,
  const bool primaryGeneration,
  std::vector<double>& ts, std::vector<std::array<double, 3> >& xs,
  std::vector<Point>& path, 
  std::vector<std::pair<Point, Particle> >& newParticles) {

  double x = p0.x;
  double y = p0.y;
  double z = p0.z;
  double t = p0.t;
  double en = p0.energy;
  int band = p0.band;
  double kx = p0.kx;
  double ky = p0.ky;
  double kz = p0.kz;
  bool firstPrimaryMultiplicationPending = primaryGeneration;
  const bool storePath = m_storeEndpoints || m_doInducedCharge;
  if (storePath) path.push_back(p0);
  if (signal || m_computePathLength) {
    ts.push_back(t);
    xs.push_back({x, y, z});
  }
  size_t did = 0;
  if (m_viewer) {
    if (hole) {
      did = m_viewer->NewDriftLine(Particle::Hole, 1, x, y, z); 
    } else { 
      did = m_viewer->NewDriftLine(Particle::Electron, 1, x, y, z);
    }
  }

  // Numerical prefactors in equation of motion
  const double c1 = SpeedOfLight * sqrt(2. / ElectronMass);
  const double c2 = 0.25 * c1 * c1;

  // Get the local electric field and medium.
  double ex = 0., ey = 0., ez = 0.;
  Medium* medium = nullptr;
  int status = 0;
  m_sensor->ElectricField(x, y, z, ex, ey, ez, medium, status);
  // Sign change for electrons.
  if (!hole) {
    ex = -ex;
    ey = -ey;
    ez = -ez;
  }
  if (m_debug) {
    std::cout << "    Drift line starts at (" 
              << x << ", " << y << ", " << z << ").\n"
              << "    Status: " << status << "\n";
    if (medium) std::cout << "    Medium: " << medium->GetName() << "\n";
  }

  if (status != 0 || !medium || !medium->IsDriftable() || 
      !medium->IsMicroscopic()) {
    if (m_debug) std::cout << "    Not in a valid medium.\n";
    return StatusLeftDriftMedium;
  }
  // Get the id number of the drift medium.
  auto mid = medium->GetId();
  // Get the null-collision rate.
  double fLim = medium->GetElectronNullCollisionRate(band);
  if (fLim <= 0.) {
    std::cerr << m_className 
              << "::TransportElectron: Got null-collision rate <= 0.\n";
    return StatusCalculationAbandoned;
  }
  if (m_debug) std::cout << "    Null collision rate: " << fLim << "\n";
  double tLim = 1. / fLim;

  std::vector<Medium::Secondary> secondaries;
  // Keep track of the previous coordinates for distance histogramming.
  double xLast = x;
  double yLast = y;
  double zLast = z;
  auto hEnergy = hole ? m_histHoleEnergy : m_histElectronEnergy;
  // Trace the electron/hole.
  size_t nColl = 0;
  size_t nCollPlot = 0;
  while (1) {
    // Make sure the kinetic energy exceeds the transport cut.
    if (en < m_deltaCut) {
      if (m_debug) std::cout << "    Kinetic energy below transport cut.\n";
      status = StatusBelowTransportCut;
      break;
    }

    // Fill the energy distribution histogram.
    if (hEnergy) hEnergy->Fill(en);

    // Make sure the particle is within the specified time window.
    if (m_hasTimeWindow && (t < m_tMin || t > m_tMax)) {
      if (m_debug) std::cout << "    Outside the time window.\n";
      status = StatusOutsideTimeWindow;
      break;
    }

    if (medium->GetId() != mid) {
      // Medium has changed.
      if (!medium->IsMicroscopic()) {
        // Electron/hole has left the microscopic drift medium.
        if (m_debug) std::cout << "    Not in a microscopic medium.\n";
        status = StatusLeftDriftMedium;
        break;
      }
      mid = medium->GetId();
      // Update the null-collision rate.
      fLim = medium->GetElectronNullCollisionRate(band);
      if (fLim <= 0.) {
        std::cerr << m_className 
                  << "::TransportElectron: Got null-collision rate <= 0.\n";
        status = StatusCalculationAbandoned;
        break;
      }
      if (m_debug) std::cout << "    Null collision rate: " << fLim << "\n";
      tLim = 1. / fLim;
    }

    // Calculate the initial velocity vector.
    const double vmag = c1 * sqrt(en);
    double vx = vmag * kx;
    double vy = vmag * ky;
    double vz = vmag * kz;
    const double a1 = vx * ex + vy * ey + vz * ez;
    const double a2 = c2 * (ex * ex + ey * ey + ez * ez);

    if (m_userHandleStep) {
      m_userHandleStep(x, y, z, t, en, kx, ky, kz, hole);
    }

    // Variables for the RKN method.
    double r0[3], vr[3];
    std::vector<std::array<double, 3> > rknIntPoints = {};
    int nsteps = 2;
    
    // Energy after the step.
    double en1 = en;
    // Determine the timestep.
    double dt = 0.;
    bool isNullCollision = true;
    tLim *= m_nullCollScale;
    while (isNullCollision) {
      // Sample the flight time.
      const double r = RndmUniformPos();
      dt += -log(r) * tLim;
      // Calculate the energy after the proposed step.
      if (m_rknSteps) {
        if (m_debug) {
          std::cout << "\n=============================\n\n"
                    << "RKN: (x,y,z) = (" << x << ", " << y << ", " << z << ")\n";
        }
      
        double h = dt / nsteps;
        nsteps = 0;
        double timeholder = 0.;
        bool loopholder = true;
        double ex0 = ex, ey0 = ey, ez0 = ez;
        r0[0] = x;
        r0[1] = y;
        r0[2] = z;
        vr[0] = vx;
        vr[1] = vy;
        vr[2] = vz;
        const double c3 = -2 * c2;
        std::array<double, 3> k1 = {c3 * ex0, c3 * ey0, c3 * ez0};
        rknIntPoints.clear();
        while (loopholder) {
          // Check if it is the last step of the time interval dt
          if (0 > dt - timeholder - h) {
            loopholder = false;
            h = dt - timeholder;
            if (m_debug) std::cout  << "\n";
          }
          
          timeholder += h;
          if (m_debug) {
            std::cout << "RKN: Time keeper = " << timeholder
                      << " of the " << dt << " ns.\n";
          }
          const double h2 = h * h;
          Medium* med0 = nullptr;
          int stat0 = 0;
          m_sensor->ElectricField(r0[0] + h * vr[0] * 0.5 + 0.125 * h2 * k1[0],
                                  r0[1] + h * vr[1] * 0.5 + 0.125 * h2 * k1[1],
                                  r0[2] + h * vr[2] * 0.5 + 0.125 * h2 * k1[2],
                                  ex0, ey0, ez0, med0, stat0);
          
          std::array<double, 3> k2 = {c3 * ex0, c3 * ey0, c3 * ez0}; // k3 = k2
          m_sensor->ElectricField(r0[0] + h * vr[0] + h2 * k2[0] * 0.5,
                                  r0[1] + h * vr[1] + h2 * k2[1] * 0.5,
                                  r0[2] + h * vr[2] + h2 * k2[2] * 0.5,
                                  ex0, ey0, ez0, med0, stat0);
          std::array<double, 3> k4 = {c3 * ex0, c3 * ey0, c3 * ez0};
          
          // Check error tolerance
          const double steperror = h2 * (
            Mag(k1[0], k1[1], k1[2]) - 2 * Mag(k2[0], k2[1], k2[2]) + 
            Mag(k4[0], k4[1], k4[2]));
          
          if (m_debug) std::cout << "RKN: steperror = " << steperror << ".\n";
          
          if (std::abs(steperror) < 4 * m_rknsteperrortol) {
            for (int j = 0; j <= 2; j++) {
              r0[j] += h * vr[j] + (k1[j] + k2[j] + k2[j]) * h2 / 6.;
              vr[j] += (k1[j] + 4 * k2[j] + k4[j]) * h / 6.;
            }
            
            const double hholder = h;
            h *= pow((m_rknsteperrortol / std::abs(steperror)), 0.25);
            
            if ((0.25 * hholder <= h && 4 * hholder >= h) ||
                steperror < m_rknsteperrortol * 1.e-10 || h < m_rknMinh) {
              h = hholder;
            }
            
            // Final point of current stage is first point of the next
            k1.swap(k4);
            
            if (stat0 != 0) {
              if (m_debug) {
                std::cout << "RKN: Outside drift medium! Breaking loop.\n";
              }
              dt = timeholder + h;
              break;
            }
            
            // Plot intermediate points
            if (m_viewer) rknIntPoints.push_back({r0[0], r0[1], r0[2]});
            nsteps++;
            
          } else {
            timeholder -= h;
            // Adjust step size
            h *= pow((m_rknsteperrortol / std::abs(steperror)), 0.25);
          }
          if (m_debug) std::cout << "RKN: h = " << h << "\n";
        }
        
        en1 = std::max(
                       (vr[0] * vr[0] + vr[1] * vr[1] + vr[2] * vr[2]) / (c1 * c1),
                       Small);
      } else {
        en1 = std::max(en + (a1 + a2 * dt) * dt, Small);
      }
      
      if (m_debug) {
        std::cout << "RKN: en1 = " << en1 << ","
                  << std::max(en + (a1 + a2 * dt) * dt, Small) << " eV.\n";
      }
      // Get the real collision rate at the updated energy.
      const double fReal = medium->GetElectronCollisionRate(en1, band);
      if (fReal <= 0.) {
        std::cerr << m_className << "::TransportElectron:\n"
                  << "    Got collision rate <= 0 at " << en1
                  << " eV (band " << band << ").\n";
        if (storePath) {
          path.emplace_back(MakePoint(x, y, z, t, en1, kx, ky, kz, band));
        }
        return StatusCalculationAbandoned;
      }
      if (fReal > fLim) {
        // Real collision rate is higher than null-collision rate.
        dt += log(r) * tLim;
        // Increase the null collision rate and try again.
        std::cerr << m_className << "::TransportElectron: "
                  << "Increasing null-collision rate by 5%.\n";
        fLim *= 1.05;
        tLim = 1. / fLim;
        continue;
      }
      // Check for real or null collision.
      if (RndmUniform() <= fReal * tLim) isNullCollision = false;
      if (m_useNullCollisionSteps) break;
    }

    // Increase the collision counters.
    ++nColl;
    ++nCollPlot;

    double x1, y1, z1;
    double kx1, ky1, kz1;
    if (m_rknSteps) {
      // Update the direction.
      const double a1 = 1/ (c1 * sqrt(en1));
      
      kx1 = vr[0] * a1;
      ky1 = vr[1] * a1;
      kz1 = vr[2] * a1;
      
      vx = vr[0];
      vy = vr[1];
      vz = vr[2];
      
      // Update the step in coordinate space.
      x1 = r0[0];
      y1 = r0[1];
      z1 = r0[2];
    } else {
      // Calculate the direction at the instant before the collision.
      const double b1 = sqrt(en / en1);
      const double b2 = 0.5 * c1 * dt / sqrt(en1);
      kx1 = kx * b1 + ex * b2;
      ky1 = ky * b1 + ey * b2;
      kz1 = kz * b1 + ez * b2;
      
      // Calculate the step in coordinate space.
      const double b3 = dt * dt * c2;
      x1 = x + vx * dt + ex * b3;
      y1 = y + vy * dt + ey * b3;
      z1 = z + vz * dt + ez * b3;
    }
    double t1 = t + dt;
    
    // Get the electric field and medium at the proposed new position.
    m_sensor->ElectricField(x1, y1, z1, ex, ey, ez, medium, status);
    
    if (!hole) {
      ex = -ex;
      ey = -ey;
      ez = -ez;
    }

    double xc = x, yc = y, zc = z, rc = 0.;
    // Is the particle still inside a drift medium/the drift area?
    if (status != 0) {
      // Try to terminate the drift line close to the boundary (endpoint
      // outside the drift medium/drift area) using iterative bisection.
      Terminate(x, y, z, t, x1, y1, z1, t1);
      if (m_debug) std::cout << "    Left the drift medium.\n";
      status = StatusLeftDriftMedium;
    } else if (!m_sensor->IsInArea(x1, y1, z1)) {
      Terminate(x, y, z, t, x1, y1, z1, t1);
      if (m_debug) std::cout << "    Left the drift area.\n";
      status = StatusLeftDriftArea;
    } else if (m_sensor->CrossedWire(x, y, z, x1, y1, z1, 
                                     xc, yc, zc, false, rc)) {
      t1 = t + dt * Mag(xc - x, yc - y, zc - z) /
                    Mag(x1 - x, y1 - y, z1 - z);
      x1 = xc;
      y1 = yc;
      z1 = zc;
      if (m_debug) std::cout << "    Hit a wire.\n";
      status = StatusLeftDriftMedium;
    } else if (m_sensor->CrossedPlane(x, y, z, x1, y1, z1, xc, yc, zc)) {
      t1 = t + dt * Mag(xc - x, yc - y, zc - z) /
                    Mag(x1 - x, y1 - y, z1 - z);
      x1 = xc;
      y1 = yc;
      z1 = zc;
      if (m_debug) std::cout << "    Hit a plane.\n";
      status = StatusHitPlane;
    }

    if (signal || m_computePathLength) {
      ts.push_back(t1);
      xs.push_back({x1, y1, z1});
    }
    // Update the coordinates.
    x = x1;
    y = y1;
    z = z1;
    t = t1;

    if (status != 0) {
      en = en1;
      kx = kx1;
      ky = ky1;
      kz = kz1;
      break;
    }

    if (isNullCollision) {
      en = en1;
      kx = kx1;
      ky = ky1;
      kz = kz1;
      continue;
    }

    // Get the collision type and parameters.
    int cstype = 0;
    int level = 0;
    secondaries.clear();
    medium->ElectronCollision(en1, cstype, level, en, kx1, ky1, kz1,
                              secondaries, band);
    if (m_debug) std::cout << "    Collision type " << cstype << ".\n";
    // If activated, histogram the distance with respect to the
    // last collision.
    if (m_histDistance) {
      FillDistanceHistogram(cstype, x, y, z, xLast, yLast, zLast);
    }

    if (m_userHandleCollision) {
      m_userHandleCollision(x, y, z, t, cstype, level, medium, en1, en, kx,
                            ky, kz, kx1, ky1, kz1);
    }
    switch (cstype) {
      // Elastic collision
      case ElectronCollisionTypeElastic:
        break;
      // Ionising collision
      case ElectronCollisionTypeIonisation:
        if (m_userHandleIonisation) {
          m_userHandleIonisation(x, y, z, t, cstype, level, medium);
        }
        if (firstPrimaryMultiplicationPending &&
            m_firstPrimaryMultiplicationScale < 1. &&
            RndmUniform() > m_firstPrimaryMultiplicationScale) {
          break;
        }
        if (firstPrimaryMultiplicationPending) {
          firstPrimaryMultiplicationPending = false;
        }
        for (const auto& secondary : secondaries) {
          if (secondary.type == Particle::Electron) {
            const double esec = std::max(secondary.energy, Small);
            if (m_histSecondary) m_histSecondary->Fill(esec);
            // Add the secondary electron to the stack.
            newParticles.emplace_back(std::make_pair(
              MakePoint(x, y, z, t, esec), Particle::Electron));
          } else if (secondary.type == Particle::Hole) {
            const double esec = std::max(secondary.energy, Small);
            // Add the secondary hole to the stack.
            newParticles.emplace_back(std::make_pair(
              MakePoint(x, y, z, t, esec), Particle::Hole));
          } else if (secondary.type == Particle::Ion) {
            newParticles.emplace_back(std::make_pair(
              MakePoint(x, y, z, t, 0.), Particle::Ion));
          }
        }
        break;
      // Attachment
      case ElectronCollisionTypeAttachment:
        if (m_userHandleAttachment) {
          m_userHandleAttachment(x, y, z, t, cstype, level, medium);
        }
        if (storePath) {
          path.emplace_back(MakePoint(x, y, z, t, en, kx1, ky1, kz1, band));
        }
        return StatusAttached;
        break;
      // Inelastic collision
      case ElectronCollisionTypeInelastic:
        if (m_userHandleInelastic) {
          m_userHandleInelastic(x, y, z, t, cstype, level, medium);
        }
        break;
      // Excitation
      case ElectronCollisionTypeExcitation:
        if (m_userHandleInelastic) {
          m_userHandleInelastic(x, y, z, t, cstype, level, medium);
        }
        // Loop over the particles produced in the deexcitation cascade.
        for (const auto& secondary : secondaries) {
          if (secondary.type == Particle::Electron) {
            if (firstPrimaryMultiplicationPending &&
                m_firstPrimaryMultiplicationScale < 1. &&
                RndmUniform() > m_firstPrimaryMultiplicationScale) {
              continue;
            }
            if (firstPrimaryMultiplicationPending) {
              firstPrimaryMultiplicationPending = false;
            }
            // Penning ionisation
            double xp = x, yp = y, zp = z;
            if (secondary.distance > Small) {
              // Randomise the point of creation.
              double dxp = 0., dyp = 0., dzp = 0.;
              RndmDirection(dxp, dyp, dzp);
              xp += secondary.distance * dxp;
              yp += secondary.distance * dyp;
              zp += secondary.distance * dzp;
            }
            // Get the electric field and medium at this location.
            Medium* med = nullptr;
            double fx = 0., fy = 0., fz = 0.;
            m_sensor->ElectricField(xp, yp, zp, fx, fy, fz, med, status);
            // Check if this location is inside a drift medium/area.
            if (status != 0 || !m_sensor->IsInArea(xp, yp, zp)) continue;
            // Make sure we haven't jumped across a wire.
            if (m_sensor->CrossedWire(x, y, z, xp, yp, zp, xc, yc, zc,
                                      false, rc)) {
              continue;
            }
            if (m_userHandleIonisation) {
              m_userHandleIonisation(xp, yp, zp, t, cstype, level, medium);
            }
            // Add the Penning electron to the list.
            const double tp = t + secondary.time;
            const double ep = std::max(secondary.energy, Small);
            newParticles.emplace_back(std::make_pair(
              MakePoint(xp, yp, zp, tp, ep), Particle::Electron));
            newParticles.emplace_back(std::make_pair(
              MakePoint(xp, yp, zp, tp, 0.), Particle::Ion));
          } else if (secondary.type == Particle::Photon && m_usePhotons &&
                     secondary.energy > m_gammaCut) {
            // Radiative de-excitation
            if (aval) TransportPhoton(x, y, z, t + secondary.time, 
                                      secondary.energy, newParticles);
          }
        }
        break;
      // Super-elastic collision
      case ElectronCollisionTypeSuperelastic:
      // Virtual/null collision
      case ElectronCollisionTypeVirtual:
      // Acoustic phonon scattering (intravalley)
      case ElectronCollisionTypeAcousticPhonon:
      // Optical phonon scattering (intravalley)
      case ElectronCollisionTypeOpticalPhonon:
      // Intervalley scattering (phonon assisted)
      case ElectronCollisionTypeIntervalleyG:
      case ElectronCollisionTypeIntervalleyF:
      case ElectronCollisionTypeInterbandXL:
      case ElectronCollisionTypeInterbandXG:
      case ElectronCollisionTypeInterbandLG:
      // Coulomb scattering
      case ElectronCollisionTypeImpurity:
        break;
      default:
        std::cerr << m_className 
                  << "::TransportElectron: Unknown collision type.\n";
        break;
    }
    if (m_viewer) {
      if (m_rknSteps) {
        for (const auto& pt : rknIntPoints) {
          PlotCollision(cstype, did, pt[0], pt[1], pt[2], nCollPlot);
        }
      } else {
        PlotCollision(cstype, did, x, y, z, nCollPlot);
      }
    }

    // Update the direction vector.
    kx = kx1;
    ky = ky1;
    kz = kz1;

    // Normalise the direction vector.
    if (nColl % 100 == 0) Normalise(kx, ky, kz);

    // Add a new point to the drift line (if enabled).
    if (storePath && m_storeDriftLines && nColl >= m_nCollSkip) {
      path.emplace_back(MakePoint(x, y, z, t, en, kx, ky, kz, band));
      nColl = 0;
    }
    if (m_viewer && nCollPlot >= m_nCollPlot) {
      m_viewer->AddDriftLinePoint(did, x, y, z);
      nCollPlot = 0;
    }
  }

  if (storePath && nColl > 0) {
    path.emplace_back(MakePoint(x, y, z, t, en, kx, ky, kz, band));
  }
  if (m_viewer && nCollPlot > 0) {
    m_viewer->AddDriftLinePoint(did, x, y, z);
  }
  if (m_debug) {
    std::cout << "    Drift line stops at (" 
              << x << ", " << y << ", " << z << ").\n";
  } 
  return status;
}

int AvalancheMicroscopic::TransportElectronBfield(const Point& p0,
  const bool hole, const bool aval, const bool signal,
  const bool primaryGeneration,
  std::vector<double>& ts, std::vector<std::array<double, 3> >& xs,
  std::vector<Point>& path, 
  std::vector<std::pair<Point, Particle> >& newParticles) {

  double x = p0.x;
  double y = p0.y;
  double z = p0.z;
  double t = p0.t;
  double en = p0.energy;
  int band = p0.band;
  double kx = p0.kx;
  double ky = p0.ky;
  double kz = p0.kz;
  bool firstPrimaryMultiplicationPending = primaryGeneration;
  const bool storePath = m_storeEndpoints || m_doInducedCharge;
  if (storePath) path.push_back(p0);
  if (signal || m_computePathLength) {
    ts.push_back(t);
    xs.push_back({x, y, z});
  }
  size_t did = 0;
  if (m_viewer) {
    if (hole) {
      did = m_viewer->NewDriftLine(Particle::Hole, 1, x, y, z); 
    } else { 
      did = m_viewer->NewDriftLine(Particle::Electron, 1, x, y, z);
    }
  }

  // Numerical prefactors in equation of motion
  const double c1 = SpeedOfLight * sqrt(2. / ElectronMass);
  const double c2 = 0.25 * c1 * c1;

  // Get the local electric field and medium.
  double ex = 0., ey = 0., ez = 0.;
  Medium* medium = nullptr;
  int status = 0;
  m_sensor->ElectricField(x, y, z, ex, ey, ez, medium, status);
  // Sign change for electrons.
  if (!hole) {
    ex = -ex;
    ey = -ey;
    ez = -ez;
  }
  if (m_debug) {
    std::cout << "    Drift line starts at (" 
              << x << ", " << y << ", " << z << ").\n"
              << "    Status: " << status << "\n";
    if (medium) std::cout << "    Medium: " << medium->GetName() << "\n";
  }

  if (status != 0 || !medium || !medium->IsDriftable() || 
      !medium->IsMicroscopic()) {
    if (m_debug) std::cout << "    Not in a valid medium.\n";
    return StatusLeftDriftMedium;
  }
  // Get the id number of the drift medium.
  auto mid = medium->GetId();
  // Get the null-collision rate.
  double fLim = medium->GetElectronNullCollisionRate(band);
  if (fLim <= 0.) {
    std::cerr << m_className 
              << "::TransportElectron: Got null-collision rate <= 0.\n";
    return StatusCalculationAbandoned;
  }
  double tLim = 1. / fLim;

  std::array<std::array<double, 3>, 3> rot;
  // Get the local magnetic field.
  double bx = 0., by = 0., bz = 0.;
  int st = 0;
  m_sensor->MagneticField(x, y, z, bx, by, bz, st);
  const double scale = hole ? Tesla2Internal : -Tesla2Internal;
  bx *= scale;
  by *= scale;
  bz *= scale;
  double bmag = Mag(bx, by, bz);
  // Calculate the rotation matrix to a local coordinate system
  // with B along x and E in the x-z plane.
  RotationMatrix(bx, by, bz, bmag, ex, ey, ez, rot);
  // Calculate the cyclotron frequency.
  double omega = OmegaCyclotronOverB * bmag;
  // Calculate the electric field in the local frame.
  ToLocal(rot, ex, ey, ez, ex, ey, ez);
  // Ratio of transverse electric field component and magnetic field.
  double ezovb = bmag > Small ? ez / bmag : 0.;

  std::vector<Medium::Secondary> secondaries;
  // Keep track of the previous coordinates for distance histogramming.
  double xLast = x;
  double yLast = y;
  double zLast = z;
  auto hEnergy = hole ? m_histHoleEnergy : m_histElectronEnergy;
  // Trace the electron/hole.
  size_t nColl = 0;
  size_t nCollPlot = 0;
  while (1) {
    // Make sure the kinetic energy exceeds the transport cut.
    if (en < m_deltaCut) {
      if (m_debug) std::cout << "    Kinetic energy below transport cut.\n";
      status = StatusBelowTransportCut;
      break;
    }

    // Fill the energy distribution histogram.
    if (hEnergy) hEnergy->Fill(en);

    // Make sure the particle is within the specified time window.
    if (m_hasTimeWindow && (t < m_tMin || t > m_tMax)) {
      if (m_debug) std::cout << "    Outside the time window.\n";
      status = StatusOutsideTimeWindow;
      break;
    }

    if (medium->GetId() != mid) {
      // Medium has changed.
      if (!medium->IsMicroscopic()) {
        // Electron/hole has left the microscopic drift medium.
        if (m_debug) std::cout << "    Not in a microscopic medium.\n";
        status = StatusLeftDriftMedium;
        break;
      }
      mid = medium->GetId();
      // Update the null-collision rate.
      fLim = medium->GetElectronNullCollisionRate(band);
      if (fLim <= 0.) {
        std::cerr << m_className 
                  << "::TransportElectron: Got null-collision rate <= 0.\n";
        status = StatusCalculationAbandoned;
        break;
      }
      tLim = 1. / fLim;
    }

    // Calculate the initial velocity vector in the local frame.
    const double vmag = c1 * sqrt(en);
    double vx = 0., vy = 0., vz = 0.;
    ToLocal(rot, vmag * kx, vmag * ky, vmag * kz, vx, vy, vz);
    double a1 = vx * ex;
    double a2 = c2 * ex * ex;
    if (omega > Small) {
      vy -= ezovb;
    } else {
      a1 += vz * ez;
      a2 += c2 * ez * ez;
    }

    if (m_userHandleStep) {
      m_userHandleStep(x, y, z, t, en, kx, ky, kz, hole);
    }

    // Energy after the step.
    double en1 = en;
    // Determine the timestep.
    double dt = 0.;
    // Parameters for B-field stepping.
    double cphi = 1., sphi = 0.;
    double a3 = 0., a4 = 0.;
    bool isNullCollision = true;
    while (isNullCollision) {
      // Sample the flight time.
      const double r = RndmUniformPos();
      dt += -log(r) * tLim;
      // Calculate the energy after the proposed step.
      en1 = en + (a1 + a2 * dt) * dt;
      if (omega > Small) {
        const double phi = omega * dt;
        cphi = cos(phi);
        sphi = sin(phi);
        a3 = sphi / omega;
        a4 = (1. - cphi) / omega;
        en1 += ez * (vz * a3 - vy * a4);
      }
      en1 = std::max(en1, Small);
      // Get the real collision rate at the updated energy.
      const double fReal = medium->GetElectronCollisionRate(en1, band);
      if (fReal <= 0.) {
        std::cerr << m_className << "::TransportElectron:\n"
                  << "    Got collision rate <= 0 at " << en1
                  << " eV (band " << band << ").\n";
        if (storePath) {
          path.emplace_back(MakePoint(x, y, z, t, en1, kx, ky, kz, band));
        }
        return StatusCalculationAbandoned;
      }
      if (fReal > fLim) {
        // Real collision rate is higher than null-collision rate.
        dt += log(r) * tLim;
        // Increase the null collision rate and try again.
        std::cerr << m_className << "::TransportElectron: "
                  << "Increasing null-collision rate by 5%.\n";
        fLim *= 1.05;
        tLim = 1. / fLim;
        continue;
      }
      // Check for real or null collision.
      if (RndmUniform() <= fReal * tLim) isNullCollision = false;
      if (m_useNullCollisionSteps) break;
    }

    // Increase the collision counters.
    ++nColl;
    ++nCollPlot;

    // Calculate the direction at the instant before the collision
    // and the proposed new position.
    double kx1 = 0., ky1 = 0., kz1 = 0.;
    double dx = 0., dy = 0., dz = 0.;
    // Calculate the new velocity.
    double vx1 = vx + 2. * c2 * ex * dt;
    double vy1 = vy * cphi + vz * sphi + ezovb;
    double vz1 = vz * cphi - vy * sphi;
    if (omega < Small) vz1 += 2. * c2 * ez * dt;
    // Rotate back to the global frame and normalise.
    ToGlobal(rot, vx1, vy1, vz1, kx1, ky1, kz1);
    Normalise(kx1, ky1, kz1);
    // Calculate the step in coordinate space.
    dx = vx * dt + c2 * ex * dt * dt;
    if (omega > Small) {
      dy = vy * a3 + vz * a4 + ezovb * dt;
      dz = vz * a3 - vy * a4;
    } else {
      dy = vy * dt;
      dz = vz * dt + c2 * ez * dt * dt;
    }
    // Rotate back to the global frame.
    ToGlobal(rot, dx, dy, dz, dx, dy, dz);

    double x1 = x + dx;
    double y1 = y + dy;
    double z1 = z + dz;
    double t1 = t + dt;
    // Get the electric field and medium at the proposed new position.
    m_sensor->ElectricField(x1, y1, z1, ex, ey, ez, medium, status);
    if (!hole) {
      ex = -ex;
      ey = -ey;
      ez = -ez;
    }

    double xc = x, yc = y, zc = z, rc = 0.;
    // Is the particle still inside a drift medium/the drift area?
    if (status != 0) {
      // Try to terminate the drift line close to the boundary (endpoint
      // outside the drift medium/drift area) using iterative bisection.
      Terminate(x, y, z, t, x1, y1, z1, t1);
      if (m_debug) std::cout << "    Left the drift medium.\n";
      status = StatusLeftDriftMedium;
    } else if (!m_sensor->IsInArea(x1, y1, z1)) {
      Terminate(x, y, z, t, x1, y1, z1, t1);
      if (m_debug) std::cout << "    Left the drift area.\n";
      status = StatusLeftDriftArea;
    } else if (m_sensor->CrossedWire(x, y, z, x1, y1, z1, 
                                     xc, yc, zc, false, rc)) {
      t1 = t + dt * Mag(xc - x, yc - y, zc - z) / Mag(dx, dy, dz);
      x1 = xc;
      y1 = yc;
      z1 = zc;
      if (m_debug) std::cout << "    Hit a wire.\n";
      status = StatusLeftDriftMedium;
    } else if (m_sensor->CrossedPlane(x, y, z, x1, y1, z1, xc, yc, zc)) {
      t1 = t + dt * Mag(xc - x, yc - y, zc - z) / Mag(dx, dy, dz);
      x1 = xc;
      y1 = yc;
      z1 = zc;
      if (m_debug) std::cout << "    Hit a plane.\n";
      status = StatusHitPlane;
    }

    if (signal || m_computePathLength) {
      ts.push_back(t1);
      xs.push_back({x1, y1, z1});
    }
    // Update the coordinates.
    x = x1;
    y = y1;
    z = z1;
    t = t1;

    if (status != 0) {
      en = en1;
      kx = kx1;
      ky = ky1;
      kz = kz1;
      break;
    }

    // Get the magnetic field at the new location.
    m_sensor->MagneticField(x, y, z, bx, by, bz, st);
    bx *= scale;
    by *= scale;
    bz *= scale;
    bmag = Mag(bx, by, bz);
    // Update the rotation matrix.
    RotationMatrix(bx, by, bz, bmag, ex, ey, ez, rot);
    omega = OmegaCyclotronOverB * bmag;
    // Calculate the electric field in the local frame.
    ToLocal(rot, ex, ey, ez, ex, ey, ez);
    ezovb = bmag > Small ? ez / bmag : 0.;

    if (isNullCollision) {
      en = en1;
      kx = kx1;
      ky = ky1;
      kz = kz1;
      continue;
    }

    // Get the collision type and parameters.
    int cstype = 0;
    int level = 0;
    secondaries.clear();
    medium->ElectronCollision(en1, cstype, level, en, kx1, ky1, kz1,
                              secondaries, band);

    if (m_debug) std::cout << "    Collision type " << cstype << ".\n";
    // If activated, histogram the distance with respect to the
    // last collision.
    if (m_histDistance) {
      FillDistanceHistogram(cstype, x, y, z, xLast, yLast, zLast);
    }

    if (m_userHandleCollision) {
      m_userHandleCollision(x, y, z, t, cstype, level, medium, en1, en, kx,
                            ky, kz, kx1, ky1, kz1);
    }
    switch (cstype) {
      // Elastic collision
      case ElectronCollisionTypeElastic:
        break;
      // Ionising collision
      case ElectronCollisionTypeIonisation:
        if (m_userHandleIonisation) {
          m_userHandleIonisation(x, y, z, t, cstype, level, medium);
        }
        if (firstPrimaryMultiplicationPending &&
            m_firstPrimaryMultiplicationScale < 1. &&
            RndmUniform() > m_firstPrimaryMultiplicationScale) {
          break;
        }
        if (firstPrimaryMultiplicationPending) {
          firstPrimaryMultiplicationPending = false;
        }
        for (const auto& secondary : secondaries) {
          if (secondary.type == Particle::Electron) {
            const double esec = std::max(secondary.energy, Small);
            if (m_histSecondary) m_histSecondary->Fill(esec);
            // Add the secondary electron to the stack.
            newParticles.emplace_back(std::make_pair(
              MakePoint(x, y, z, t, esec), Particle::Electron));
          } else if (secondary.type == Particle::Hole) {
            const double esec = std::max(secondary.energy, Small);
            // Add the secondary hole to the stack.
            newParticles.emplace_back(std::make_pair(
              MakePoint(x, y, z, t, esec), Particle::Hole));
          } else if (secondary.type == Particle::Ion) {
            newParticles.emplace_back(std::make_pair(
              MakePoint(x, y, z, t, 0.), Particle::Ion));
          }
        }
        break;
      // Attachment
      case ElectronCollisionTypeAttachment:
        if (m_userHandleAttachment) {
          m_userHandleAttachment(x, y, z, t, cstype, level, medium);
        }
        if (storePath) {
          path.emplace_back(MakePoint(x, y, z, t, en, kx1, ky1, kz1, band));
        }
        return StatusAttached;
        break;
      // Inelastic collision
      case ElectronCollisionTypeInelastic:
        if (m_userHandleInelastic) {
          m_userHandleInelastic(x, y, z, t, cstype, level, medium);
        }
        break;
      // Excitation
      case ElectronCollisionTypeExcitation:
        if (m_userHandleInelastic) {
          m_userHandleInelastic(x, y, z, t, cstype, level, medium);
        }
        // Loop over the particles produced in the deexcitation cascade.
        for (const auto& secondary : secondaries) {
          if (secondary.type == Particle::Electron) {
            if (firstPrimaryMultiplicationPending &&
                m_firstPrimaryMultiplicationScale < 1. &&
                RndmUniform() > m_firstPrimaryMultiplicationScale) {
              continue;
            }
            if (firstPrimaryMultiplicationPending) {
              firstPrimaryMultiplicationPending = false;
            }
            // Penning ionisation
            double xp = x, yp = y, zp = z;
            if (secondary.distance > Small) {
              // Randomise the point of creation.
              double dxp = 0., dyp = 0., dzp = 0.;
              RndmDirection(dxp, dyp, dzp);
              xp += secondary.distance * dxp;
              yp += secondary.distance * dyp;
              zp += secondary.distance * dzp;
            }
            // Get the electric field and medium at this location.
            Medium* med = nullptr;
            double fx = 0., fy = 0., fz = 0.;
            m_sensor->ElectricField(xp, yp, zp, fx, fy, fz, med, status);
            // Check if this location is inside a drift medium/area.
            if (status != 0 || !m_sensor->IsInArea(xp, yp, zp)) continue;
            // Make sure we haven't jumped across a wire.
            if (m_sensor->CrossedWire(x, y, z, xp, yp, zp, xc, yc, zc,
                                      false, rc)) {
              continue;
            }
            if (m_userHandleIonisation) {
              m_userHandleIonisation(xp, yp, zp, t, cstype, level, medium);
            }
            // Add the Penning electron to the list.
            const double tp = t + secondary.time;
            const double ep = std::max(secondary.energy, Small);
            newParticles.emplace_back(std::make_pair(
              MakePoint(xp, yp, zp, tp, ep), Particle::Electron));
            newParticles.emplace_back(std::make_pair(
              MakePoint(xp, yp, zp, tp, 0.), Particle::Ion));
          } else if (secondary.type == Particle::Photon && m_usePhotons &&
                     secondary.energy > m_gammaCut) {
            // Radiative de-excitation
            if (aval) TransportPhoton(x, y, z, t + secondary.time, 
                                      secondary.energy, newParticles);
          }
        }
        break;
      // Super-elastic collision
      case ElectronCollisionTypeSuperelastic:
      // Virtual/null collision
      case ElectronCollisionTypeVirtual:
      // Acoustic phonon scattering (intravalley)
      case ElectronCollisionTypeAcousticPhonon:
      // Optical phonon scattering (intravalley)
      case ElectronCollisionTypeOpticalPhonon:
      // Intervalley scattering (phonon assisted)
      case ElectronCollisionTypeIntervalleyG:
      case ElectronCollisionTypeIntervalleyF:
      case ElectronCollisionTypeInterbandXL:
      case ElectronCollisionTypeInterbandXG:
      case ElectronCollisionTypeInterbandLG:
      // Coulomb scattering
      case ElectronCollisionTypeImpurity:
        break;
      default:
        std::cerr << m_className 
                  << "::TransportElectron: Unknown collision type.\n";
        break;
    }
    if (m_viewer) PlotCollision(cstype, did, x, y, z, nCollPlot);

    // Update the direction vector.
    kx = kx1;
    ky = ky1;
    kz = kz1;

    // Normalise the direction vector.
    if (nColl % 100 == 0) Normalise(kx, ky, kz);

    // Add a new point to the drift line (if enabled).
    if (storePath && m_storeDriftLines && nColl >= m_nCollSkip) {
      path.emplace_back(MakePoint(x, y, z, t, en, kx, ky, kz, band));
      nColl = 0;
    }
    if (m_viewer && nCollPlot >= m_nCollPlot) {
      m_viewer->AddDriftLinePoint(did, x, y, z);
      nCollPlot = 0;
    }
  }

  if (storePath && nColl > 0) {
    path.emplace_back(MakePoint(x, y, z, t, en, kx, ky, kz, band));
  }
  if (m_viewer && nCollPlot > 0) {
    m_viewer->AddDriftLinePoint(did, x, y, z);
  }
  if (m_debug) {
    std::cout << "    Drift line stops at (" 
              << x << ", " << y << ", " << z << ").\n";
  } 
  return status;
}

int AvalancheMicroscopic::TransportElectronSc(const Point& p0,
  const bool hole, const bool aval, const bool signal,
  const bool primaryGeneration,
  std::vector<double>& ts, std::vector<std::array<double, 3> >& xs,
  std::vector<Point>& path, 
  std::vector<std::pair<Point, Particle> >& newParticles) {

  double x = p0.x;
  double y = p0.y;
  double z = p0.z;
  double t = p0.t;
  double en = p0.energy;
  int band = p0.band;
  double kx = p0.kx;
  double ky = p0.ky;
  double kz = p0.kz;
  bool firstPrimaryMultiplicationPending = primaryGeneration;
  const bool storePath = m_storeEndpoints || m_doInducedCharge;
  if (storePath) path.push_back(p0);
  if (signal || m_computePathLength) {
    ts.push_back(t);
    xs.push_back({x, y, z});
  }
  size_t did = 0;
  if (m_viewer) {
    if (hole) {
      did = m_viewer->NewDriftLine(Particle::Hole, 1, x, y, z); 
    } else { 
      did = m_viewer->NewDriftLine(Particle::Electron, 1, x, y, z);
    }
  }

  // Get the local electric field and medium.
  double ex = 0., ey = 0., ez = 0.;
  Medium* medium = nullptr;
  int status = 0;
  m_sensor->ElectricField(x, y, z, ex, ey, ez, medium, status);
  // Sign change for electrons.
  if (!hole) {
    ex = -ex;
    ey = -ey;
    ez = -ez;
  }
  if (m_debug) {
    std::cout << "    Drift line starts at (" 
              << x << ", " << y << ", " << z << ").\n"
              << "    Status: " << status << "\n";
    if (medium) std::cout << "    Medium: " << medium->GetName() << "\n";
  }

  if (status != 0 || !medium || !medium->IsDriftable() || 
      !medium->IsMicroscopic()) {
    if (m_debug) std::cout << "    Not in a valid medium.\n";
    return StatusLeftDriftMedium;
  }
  // Get the id number of the drift medium.
  auto mid = medium->GetId();
  // Get the null-collision rate.
  double fLim = medium->GetElectronNullCollisionRate(band);
  if (fLim <= 0.) {
    std::cerr << m_className 
              << "::TransportElectron: Got null-collision rate <= 0.\n";
    return StatusCalculationAbandoned;
  }
  double tLim = 1. / fLim;

  std::vector<Medium::Secondary> secondaries;
  // Keep track of the previous coordinates for distance histogramming.
  double xLast = x;
  double yLast = y;
  double zLast = z;
  auto hEnergy = hole ? m_histHoleEnergy : m_histElectronEnergy;
  // Trace the electron/hole.
  size_t nColl = 0;
  size_t nCollPlot = 0;
  while (1) {
    // Make sure the kinetic energy exceeds the transport cut.
    if (en < m_deltaCut) {
      if (m_debug) std::cout << "    Kinetic energy below transport cut.\n";
      status = StatusBelowTransportCut;
      break;
    }

    // Fill the energy distribution histogram.
    if (hEnergy) hEnergy->Fill(en);

    // Make sure the particle is within the specified time window.
    if (m_hasTimeWindow && (t < m_tMin || t > m_tMax)) {
      if (m_debug) std::cout << "    Outside the time window.\n";
      status = StatusOutsideTimeWindow;
      break;
    }

    if (medium->GetId() != mid) {
      // Medium has changed.
      if (!medium->IsMicroscopic()) {
        // Electron/hole has left the microscopic drift medium.
        if (m_debug) std::cout << "    Not in a microscopic medium.\n";
        status = StatusLeftDriftMedium;
        break;
      }
      mid = medium->GetId();
      // TODO!
      // sc = (medium->IsSemiconductor() && m_useBandStructure);
      // Update the null-collision rate.
      fLim = medium->GetElectronNullCollisionRate(band);
      if (fLim <= 0.) {
        std::cerr << m_className 
                  << "::TransportElectron: Got null-collision rate <= 0.\n";
        status = StatusCalculationAbandoned;
        break;
      }
      tLim = 1. / fLim;
    }

    // Initial velocity.
    double vx = 0., vy = 0., vz = 0.;
    en = medium->GetElectronEnergy(kx, ky, kz, vx, vy, vz, band);

    if (m_userHandleStep) {
      m_userHandleStep(x, y, z, t, en, kx, ky, kz, hole);
    }

    // Energy after the step.
    double en1 = en;
    // Determine the timestep.
    double dt = 0.;
    bool isNullCollision = true;
    while (isNullCollision) {
      // Sample the flight time.
      const double r = RndmUniformPos();
      dt += -log(r) * tLim;
      // Calculate the energy after the proposed step.
      const double cdt = dt * SpeedOfLight;
      const double kx1 = kx + ex * cdt;
      const double ky1 = ky + ey * cdt;
      const double kz1 = kz + ez * cdt;
      double vx1 = 0., vy1 = 0., vz1 = 0.;
      en1 = medium->GetElectronEnergy(kx1, ky1, kz1, vx1, vy1, vz1, band);
      en1 = std::max(en1, Small);
      // Get the real collision rate at the updated energy.
      const double fReal = medium->GetElectronCollisionRate(en1, band);
      if (fReal <= 0.) {
        std::cerr << m_className << "::TransportElectron:\n"
                  << "    Got collision rate <= 0 at " << en1
                  << " eV (band " << band << ").\n";
        if (storePath) {
          path.emplace_back(MakePoint(x, y, z, t, en1, kx, ky, kz, band));
        }
        return StatusCalculationAbandoned;
      }
      if (fReal > fLim) {
        // Real collision rate is higher than null-collision rate.
        dt += log(r) * tLim;
        // Increase the null collision rate and try again.
        std::cerr << m_className << "::TransportElectron: "
                  << "Increasing null-collision rate by 5% (band "
                  << band << ").\n";
        fLim *= 1.05;
        tLim = 1. / fLim;
        continue;
      }
      // Check for real or null collision.
      if (RndmUniform() <= fReal * tLim) isNullCollision = false;
      if (m_useNullCollisionSteps) break;
    }

    // Increase the collision counters.
    ++nColl;
    ++nCollPlot;

    // Calculate the direction at the instant before the collision
    // and the proposed new position.
    double kx1 = 0., ky1 = 0., kz1 = 0.;
    double dx = 0., dy = 0., dz = 0.;
    // Update the wave-vector.
    const double cdt = dt * SpeedOfLight;
    kx1 = kx + ex * cdt;
    ky1 = ky + ey * cdt;
    kz1 = kz + ez * cdt;
    double vx1 = 0., vy1 = 0, vz1 = 0.;
    en1 = medium->GetElectronEnergy(kx1, ky1, kz1, vx1, vy1, vz1, band);
    dx = 0.5 * (vx + vx1) * dt;
    dy = 0.5 * (vy + vy1) * dt;
    dz = 0.5 * (vz + vz1) * dt;

    double x1 = x + dx;
    double y1 = y + dy;
    double z1 = z + dz;
    double t1 = t + dt;
    // Get the electric field and medium at the proposed new position.
    m_sensor->ElectricField(x1, y1, z1, ex, ey, ez, medium, status);
    if (!hole) {
      ex = -ex;
      ey = -ey;
      ez = -ez;
    }

    double xc = x, yc = y, zc = z, rc = 0.;
    // Is the particle still inside a drift medium/the drift area?
    if (status != 0) {
      // Try to terminate the drift line close to the boundary (endpoint
      // outside the drift medium/drift area) using iterative bisection.
      Terminate(x, y, z, t, x1, y1, z1, t1);
      if (m_debug) std::cout << "    Left the drift medium.\n";
      status = StatusLeftDriftMedium;
    } else if (!m_sensor->IsInArea(x1, y1, z1)) {
      Terminate(x, y, z, t, x1, y1, z1, t1);
      if (m_debug) std::cout << "    Left the drift area.\n";
      status = StatusLeftDriftArea;
    } else if (m_sensor->CrossedWire(x, y, z, x1, y1, z1, 
                                     xc, yc, zc, false, rc)) {
      t1 = t + dt * Mag(xc - x, yc - y, zc - z) / Mag(dx, dy, dz);
      x1 = xc;
      y1 = yc;
      z1 = zc;
      if (m_debug) std::cout << "    Hit a wire.\n";
      status = StatusLeftDriftMedium;
    } else if (m_sensor->CrossedPlane(x, y, z, x1, y1, z1, xc, yc, zc)) {
      t1 = t + dt * Mag(xc - x, yc - y, zc - z) / Mag(dx, dy, dz);
      x1 = xc;
      y1 = yc;
      z1 = zc;
      if (m_debug) std::cout << "    Hit a plane.\n";
      status = StatusHitPlane;
    }

    if (signal || m_computePathLength) {
      ts.push_back(t1);
      xs.push_back({x1, y1, z1});
    }
    // Update the coordinates.
    x = x1;
    y = y1;
    z = z1;
    t = t1;

    if (status != 0) {
      en = en1;
      kx = kx1;
      ky = ky1;
      kz = kz1;
      break;
    }

    if (isNullCollision) {
      en = en1;
      kx = kx1;
      ky = ky1;
      kz = kz1;
      continue;
    }

    // Get the collision type and parameters.
    int cstype = 0;
    int level = 0;
    secondaries.clear();
    medium->ElectronCollision(en1, cstype, level, en, kx1, ky1, kz1,
                              secondaries, band);

    if (m_debug) std::cout << "    Collision type " << cstype << ".\n";
    // If activated, histogram the distance with respect to the
    // last collision.
    if (m_histDistance) {
      FillDistanceHistogram(cstype, x, y, z, xLast, yLast, zLast);
    }

    if (m_userHandleCollision) {
      m_userHandleCollision(x, y, z, t, cstype, level, medium, en1, en, kx,
                            ky, kz, kx1, ky1, kz1);
    }
    switch (cstype) {
      // Elastic collision
      case ElectronCollisionTypeElastic:
        break;
      // Ionising collision
      case ElectronCollisionTypeIonisation:
        if (m_userHandleIonisation) {
          m_userHandleIonisation(x, y, z, t, cstype, level, medium);
        }
        if (firstPrimaryMultiplicationPending &&
            m_firstPrimaryMultiplicationScale < 1. &&
            RndmUniform() > m_firstPrimaryMultiplicationScale) {
          break;
        }
        if (firstPrimaryMultiplicationPending) {
          firstPrimaryMultiplicationPending = false;
        }
        for (const auto& secondary : secondaries) {
          if (secondary.type == Particle::Electron) {
            const double esec = std::max(secondary.energy, Small);
            if (m_histSecondary) m_histSecondary->Fill(esec);
            // Add the secondary electron to the stack.
            double kxs = 0., kys = 0., kzs = 0.;
            int bs = -1;
            medium->GetElectronMomentum(esec, kxs, kys, kzs, bs);
            newParticles.emplace_back(std::make_pair(
              MakePoint(x, y, z, t, esec, kxs, kys, kzs, bs), Particle::Electron));
          } else if (secondary.type == Particle::Hole) {
            const double esec = std::max(secondary.energy, Small);
            // Add the secondary hole to the stack.
            double kxs = 0., kys = 0., kzs = 0.;
            int bs = -1;
            medium->GetElectronMomentum(esec, kxs, kys, kzs, bs);
            newParticles.emplace_back(std::make_pair(
              MakePoint(x, y, z, t, esec, kxs, kys, kzs, bs), Particle::Hole));
          } else if (secondary.type == Particle::Ion) {
            newParticles.emplace_back(std::make_pair(
              MakePoint(x, y, z, t, 0.), Particle::Ion));
          }
        }
        break;
      // Attachment
      case ElectronCollisionTypeAttachment:
        if (m_userHandleAttachment) {
          m_userHandleAttachment(x, y, z, t, cstype, level, medium);
        }
        if (hole) {
          --m_nHoles;
        } else {
          --m_nElectrons;
        }
        if (storePath) {
          path.emplace_back(MakePoint(x, y, z, t, en, kx1, ky1, kz1, band));
        }
        return StatusAttached;
        break;
      // Inelastic collision
      case ElectronCollisionTypeInelastic:
        if (m_userHandleInelastic) {
          m_userHandleInelastic(x, y, z, t, cstype, level, medium);
        }
        break;
      // Excitation
      case ElectronCollisionTypeExcitation:
        if (m_userHandleInelastic) {
          m_userHandleInelastic(x, y, z, t, cstype, level, medium);
        }
        // Loop over the particles produced in the deexcitation cascade.
        for (const auto& secondary : secondaries) {
          if (secondary.type == Particle::Electron) {
            if (firstPrimaryMultiplicationPending &&
                m_firstPrimaryMultiplicationScale < 1. &&
                RndmUniform() > m_firstPrimaryMultiplicationScale) {
              continue;
            }
            if (firstPrimaryMultiplicationPending) {
              firstPrimaryMultiplicationPending = false;
            }
            // Penning ionisation
            double xp = x, yp = y, zp = z;
            if (secondary.distance > Small) {
              // Randomise the point of creation.
              double dxp = 0., dyp = 0., dzp = 0.;
              RndmDirection(dxp, dyp, dzp);
              xp += secondary.distance * dxp;
              yp += secondary.distance * dyp;
              zp += secondary.distance * dzp;
            }
            // Get the electric field and medium at this location.
            Medium* med = nullptr;
            double fx = 0., fy = 0., fz = 0.;
            m_sensor->ElectricField(xp, yp, zp, fx, fy, fz, med, status);
            // Check if this location is inside a drift medium/area.
            if (status != 0 || !m_sensor->IsInArea(xp, yp, zp)) continue;
            // Make sure we haven't jumped across a wire.
            if (m_sensor->CrossedWire(x, y, z, xp, yp, zp, xc, yc, zc,
                                      false, rc)) {
              continue;
            }
            if (m_userHandleIonisation) {
              m_userHandleIonisation(xp, yp, zp, t, cstype, level, medium);
            }
            // Add the Penning electron to the list.
            const double tp = t + secondary.time;
            const double ep = std::max(secondary.energy, Small);
            newParticles.emplace_back(std::make_pair(
              MakePoint(xp, yp, zp, tp, ep), Particle::Electron));
            newParticles.emplace_back(std::make_pair(
              MakePoint(xp, yp, zp, tp, 0.), Particle::Ion));
          } else if (secondary.type == Particle::Photon && m_usePhotons &&
                     secondary.energy > m_gammaCut) {
            // Radiative de-excitation
            if (aval) TransportPhoton(x, y, z, t + secondary.time, 
                                      secondary.energy, newParticles);
          }
        }
        break;
      // Super-elastic collision
      case ElectronCollisionTypeSuperelastic:
      // Virtual/null collision
      case ElectronCollisionTypeVirtual:
      // Acoustic phonon scattering (intravalley)
      case ElectronCollisionTypeAcousticPhonon:
      // Optical phonon scattering (intravalley)
      case ElectronCollisionTypeOpticalPhonon:
      // Intervalley scattering (phonon assisted)
      case ElectronCollisionTypeIntervalleyG:
      case ElectronCollisionTypeIntervalleyF:
      case ElectronCollisionTypeInterbandXL:
      case ElectronCollisionTypeInterbandXG:
      case ElectronCollisionTypeInterbandLG:
      // Coulomb scattering
      case ElectronCollisionTypeImpurity:
        break;
      default:
        std::cerr << m_className 
                  << "::TransportElectron: Unknown collision type.\n";
        break;
    }
    if (m_viewer) PlotCollision(cstype, did, x, y, z, nCollPlot);

    // Update the direction vector.
    kx = kx1;
    ky = ky1;
    kz = kz1;

    // Add a new point to the drift line (if enabled).
    if (storePath && m_storeDriftLines && nColl >= m_nCollSkip) {
      path.emplace_back(MakePoint(x, y, z, t, en, kx, ky, kz, band));
      nColl = 0;
    }
    if (m_viewer && nCollPlot >= m_nCollPlot) {
      m_viewer->AddDriftLinePoint(did, x, y, z);
      nCollPlot = 0;
    }
  }

  if (storePath && nColl > 0) {
    path.emplace_back(MakePoint(x, y, z, t, en, kx, ky, kz, band));
  }
  if (m_viewer && nCollPlot > 0) {
    m_viewer->AddDriftLinePoint(did, x, y, z);
  }
  if (m_debug) {
    std::cout << "    Drift line stops at (" 
              << x << ", " << y << ", " << z << ").\n";
  } 
  return status;
}

void AvalancheMicroscopic::PlotCollision(const int cstype, const size_t did,
    const double x, const double y, const double z,
    size_t& nCollPlot) const {
  if (!m_viewer) return;
  if (cstype == ElectronCollisionTypeIonisation) {
    if (m_plotIonisations) {
      m_viewer->AddIonisation(x, y, z);
      m_viewer->AddDriftLinePoint(did, x, y, z);
      nCollPlot = 0;
    }
  } else if (cstype == ElectronCollisionTypeExcitation) {
    if (m_plotExcitations) {
      m_viewer->AddExcitation(x, y, z);
      m_viewer->AddDriftLinePoint(did, x, y, z);
      nCollPlot = 0;
    }
  } else if (cstype == ElectronCollisionTypeAttachment) {
    if (m_plotAttachments) {
      m_viewer->AddAttachment(x, y, z);
      m_viewer->AddDriftLinePoint(did, x, y, z);
    }
  }
}

void AvalancheMicroscopic::FillDistanceHistogram(const int cstype,
    const double x, const double y, const double z,
    double& xLast, double& yLast, double& zLast) const {
 
  for (const auto& htype : m_distanceHistogramType) {
    if (htype != cstype) continue;
    if (m_debug) std::cout << "    Filling distance histogram.\n";
    switch (m_distanceOption) {
      case 'x':
        m_histDistance->Fill(xLast - x);
        break;
      case 'y':
        m_histDistance->Fill(yLast - y);
        break;
      case 'z':
        m_histDistance->Fill(zLast - z);
        break;
      case 'r':
        m_histDistance->Fill(Mag(xLast - x, yLast - y, zLast - z));
        break;
    }
    xLast = x;
    yLast = y;
    zLast = z;
    return;
  }
}

void AvalancheMicroscopic::TransportPhoton(
    const double x0, const double y0, const double z0, const double t0,
    const double e0, 
    std::vector<std::pair<Point, Particle> > & newParticles) {
  // Make sure that the sensor is defined.
  if (!m_sensor) {
    std::cerr << m_className << "::TransportPhoton: Sensor is not defined.\n";
    return;
  }

  // Make sure that the starting point is inside the active area.
  if (!m_sensor->IsInArea(x0, y0, z0)) {
    std::cerr << m_className << "::TransportPhoton:\n"
              << "    No valid field at initial position.\n";
    return;
  }

  // Make sure that the starting point is inside a medium.
  Medium* medium = m_sensor->GetMedium(x0, y0, z0);
  if (!medium) {
    std::cerr << m_className << "::TransportPhoton:\n"
              << "    No medium at initial position.\n";
    return;
  }

  // Make sure that the medium is "driftable" and microscopic.
  if (!medium->IsDriftable() || !medium->IsMicroscopic()) {
    std::cerr << m_className << "::TransportPhoton:\n"
              << "    Medium at initial position does not provide "
              << " microscopic tracking data.\n";
    return;
  }

  // Get the id number of the drift medium.
  int id = medium->GetId();

  // Position
  double x = x0, y = y0, z = z0;
  double t = t0;
  // Initial direction (randomised).
  double dx = 0., dy = 0., dz = 0.;
  RndmDirection(dx, dy, dz);
  // Energy
  double e = e0;

  // Photon collision rate
  double f = medium->GetPhotonCollisionRate(e);
  if (f <= 0.) return;
  // Timestep
  double dt = -log(RndmUniformPos()) / f;
  t += dt;
  dt *= SpeedOfLight;
  x += dt * dx;
  y += dt * dy;
  z += dt * dz;

  // Check if the photon is still inside a medium.
  medium = m_sensor->GetMedium(x, y, z);
  if (!medium || medium->GetId() != id) {
    // Try to terminate the photon track close to the boundary
    // by means of iterative bisection.
    dx *= dt;
    dy *= dt;
    dz *= dt;
    x -= dx;
    y -= dy;
    z -= dz;
    double delta = Mag(dx, dy, dz);
    if (delta > 0) {
      dx /= delta;
      dy /= delta;
      dz /= delta;
    }
    // Mid-point
    double xM = x, yM = y, zM = z;
    while (delta > BoundaryDistance) {
      delta *= 0.5;
      dt *= 0.5;
      xM = x + delta * dx;
      yM = y + delta * dy;
      zM = z + delta * dz;
      // Check if the mid-point is inside the drift medium.
      medium = m_sensor->GetMedium(xM, yM, zM);
      if (medium && medium->GetId() == id) {
        x = xM;
        y = yM;
        z = zM;
        t += dt;
      }
    }
    photon newPhoton;
    newPhoton.x0 = x0;
    newPhoton.y0 = y0;
    newPhoton.z0 = z0;
    newPhoton.x1 = x;
    newPhoton.y1 = y;
    newPhoton.z1 = z;
    newPhoton.energy = e0;
    newPhoton.status = StatusLeftDriftMedium;
    m_photons.push_back(std::move(newPhoton));
    if (m_viewer) m_viewer->AddPhoton(x0, y0, z0, x, y, z);
    return;
  }

  int type, level;
  double e1;
  double ctheta = 0.;
  std::vector<Medium::Secondary> secondaries;
  if (!medium->PhotonCollision(e, type, level, e1, ctheta, secondaries)) {
    return;
  }
  if (type == PhotonCollisionTypeIonisation) {
    for (const auto& secondary : secondaries) {
      if (secondary.type == Particle::Electron) {
        // Add the secondary electron (random direction).
        newParticles.emplace_back(std::make_pair(
          MakePoint(x, y, z, t, std::max(secondary.energy, Small)), 
          Particle::Electron));
      } else if (secondary.type == Particle::Ion) {
        newParticles.emplace_back(std::make_pair(
          MakePoint(x, y, z, t, 0.), Particle::Ion));
      }
    }
  } else if (type == PhotonCollisionTypeExcitation) {
    std::vector<double> tPhotons;
    std::vector<double> ePhotons;
    for (const auto& secondary : secondaries) {
      if (secondary.type == Particle::Electron) {
        // Ionisation.
        const double esec = std::max(secondary.energy, Small);
        newParticles.emplace_back(std::make_pair(
          MakePoint(x, y, z, t + secondary.time, esec), 
          Particle::Electron));
      } else if (secondary.type == Particle::Photon && m_usePhotons) {
        // Radiative de-excitation
        if (secondary.energy > m_gammaCut) {
          tPhotons.push_back(t + secondary.time);
          ePhotons.push_back(secondary.energy);
        }
      }
    }
    // Transport the photons (if any).
    const size_t nSizePhotons = tPhotons.size();
    for (size_t k = 0; k < nSizePhotons; ++k) {
      TransportPhoton(x, y, z, tPhotons[k], ePhotons[k], newParticles);
    }
  }

  photon newPhoton;
  newPhoton.x0 = x0;
  newPhoton.y0 = y0;
  newPhoton.z0 = z0;
  newPhoton.x1 = x;
  newPhoton.y1 = y;
  newPhoton.z1 = z;
  newPhoton.energy = e0;
  newPhoton.status = -2;
  if (m_viewer) m_viewer->AddPhoton(x0, y0, z0, x, y, z);
  m_photons.push_back(std::move(newPhoton));
}

void AvalancheMicroscopic::Terminate(double x0, double y0, double z0, double t0,
                                     double& x1, double& y1, double& z1,
                                     double& t1) const {
  const double dx = x1 - x0;
  const double dy = y1 - y0;
  const double dz = z1 - z0;
  double d = Mag(dx, dy, dz);
  while (d > BoundaryDistance) {
    d *= 0.5;
    const double xm = 0.5 * (x0 + x1);
    const double ym = 0.5 * (y0 + y1);
    const double zm = 0.5 * (z0 + z1);
    const double tm = 0.5 * (t0 + t1);
    // Check if the mid-point is inside the drift medium.
    if (m_sensor->IsInside(xm, ym, zm) && m_sensor->IsInArea(xm, ym, zm)) {
      x0 = xm;
      y0 = ym;
      z0 = zm;
      t0 = tm;
    } else {
      x1 = xm;
      y1 = ym;
      z1 = zm;
      t1 = tm;
    }
  }
  // TODO
  bool outsideMedium = true;
  while (outsideMedium) {
    d *= 0.5;
    const double xm = 0.5 * (x0 + x1);
    const double ym = 0.5 * (y0 + y1);
    const double zm = 0.5 * (z0 + z1);
    const double tm = 0.5 * (t0 + t1);
    // Check if the mid-point is inside the drift medium.
    if (m_sensor->IsInside(xm, ym, zm) && m_sensor->IsInArea(xm, ym, zm)) {
      outsideMedium = false;
    }
    x1 = xm;
    y1 = ym;
    z1 = zm;
    t1 = tm;
  }
}

void AvalancheMicroscopic::SetRunModeOptions(MPRunMode mode, int device) {
  m_runMode = mode;
  m_cudaDevice = device;
}
}  // namespace Garfield
