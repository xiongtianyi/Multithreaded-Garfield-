#include "Garfield/AvalancheGridSpaceCharge.hh"

#include <iostream>
#include <sstream>

#include "Garfield/Random.hh"

namespace {

double Mag(const double x, const double y) { return std::sqrt(x * x + y * y); }

// Get size of avalanche when going from x to x+dx in Monte Carlo fashion
void GetAvalancheSizeFromStep(double dx, const long nElectronIn,
                              const double alpha, const double eta,
                              long &nElectronOut, double &nPosIonOut,
                              double &nNegIonOut) {
  // Monte Carlo Avalanche gain per travelled distance dx (cm)
  nElectronOut = 0;
  nPosIonOut = 0;
  nNegIonOut = 0;

  if (std::abs(alpha - eta) < 1.e-8 && alpha > 1.e-8) {
    // alpha == eta
    if (nElectronIn < 1000L) {
      // Condition to which the random number will be compared.
      // If the number is smaller than the condition, nothing happens.
      // Otherwise, the single electron will be attached or retrieve
      // additional electrons from the gas.
      const double prob = alpha * dx / (1 + alpha * dx);
      // Running over all electrons in the avalanche.
      for (long i = 0; i < nElectronIn; i++) {
        // Draw a random number from the uniform distribution (0,1).
        const double s = Garfield::RndmUniformPos();
        // We (wrongly) assume if s >= prob only pos ions are created
        // else 1 neg ion.
        if (s >= prob) {
          // deviation/improvement wrt Lippmann?
          // HS: compute log(prob) outside of the loop.
          nElectronOut += (long)(log((1 - s) * (1 + alpha * dx)) / log(prob));
        } else {
          nNegIonOut += 1;
        }
      }
      // charge conservation
      nPosIonOut = (nElectronOut - nElectronIn) + nNegIonOut;

    } else {
      // Central limit theorem.
      const double sigma = sqrt(2 * alpha * dx * nElectronIn);
      nElectronOut = (long)Garfield::RndmGaussian(nElectronIn, sigma);

      // boundary conditions (alpha dx ElectronIn = dPosOut),
      //  the procedure guarantees positive values
      //  and net Nion = nPos - nNeg = nOut - nIn
      if (nElectronOut <= 0) nElectronOut = 0;  //< unphysical
      if (nElectronOut >= nElectronIn) {
        nNegIonOut = (std::exp(eta * dx) - 1) * nElectronIn;     //< >= 0
        nPosIonOut = (nElectronOut - nElectronIn) + nNegIonOut;  //< >= 0
      } else {
        nPosIonOut = (std::exp(alpha * dx) - 1) * nElectronIn;   //< >= 0
        nNegIonOut = nPosIonOut - (nElectronOut - nElectronIn);  //< >= 0
      }
    }
  } else if (alpha < 1.e-8 && eta > 1.e-8) {
    // alpha == 0, only attachment possible
    if (nElectronIn < 1000L) {
      const double prob = exp(-eta * dx);
      for (long i = 0; i < nElectronIn; i++) {
        // Draw a random number from the uniform distribution (0,1).
        const double s = Garfield::RndmUniformPos();
        if (s >= prob) {
          nNegIonOut += 1;
        } else {
          nElectronOut += 1;
        }
      }
    } else {
      // Central limit theorem.
      // HS: compute exp(-eta * dx) only once?
      const double sigma =
          std::sqrt(nElectronIn * exp(-2 * eta * dx) * (exp(-eta * dx) - 1));
      nElectronOut =
          (long)Garfield::RndmGaussian(nElectronIn * exp(-eta * dx), sigma);

      // boundary conditions
      if (nElectronOut <= 0) nElectronOut = 0;  //< unphysical
      if (nElectronOut > nElectronIn)
        nElectronOut = nElectronIn;  //< unphysical with alpha = 0

      // charge conservation
      nNegIonOut = -(nElectronOut - nElectronIn);
    }
  } else {
    // alpha != 0 =! eta
    const double k = eta / alpha;
    const double ndx = exp((alpha - eta) * dx);

    if (nElectronIn < 1000L) {
      // Condition to which the random number will be compared.
      // If the number is smaller than the condition, nothing happens.
      // Otherwise, the single electron will be attached or retrieve
      // additional electrons from the gas.
      const double prob = k * (ndx - 1) / (ndx - k);
      // Running over all electrons in the avalanche.
      for (long i = 0; i < nElectronIn; i++) {
        // Draw a random number from the uniform distribution (0,1).
        const double s = Garfield::RndmUniformPos();
        if (s >= prob) {
          // deviation/improvement wrt Lippmann?
          // HS: compute the denominator outside of the loop.
          nElectronOut +=
              (long)(1 + log((ndx - k) * (1 - s) / (ndx * (1 - k))) /
                             log(1 - (1 - k) / (ndx - k)));
        } else {
          nNegIonOut += 1;
        }
      }
      // charge conservation
      nPosIonOut = (nElectronOut - nElectronIn) + nNegIonOut;

    } else {
      // Central limit theorem.
      const double sigma =
          sqrt(nElectronIn * (1 + k) * ndx * (ndx - 1) / (1 - k));
      nElectronOut = (long)Garfield::RndmGaussian(nElectronIn * ndx, sigma);

      // boundary conditions (alpha dx ElectronIn = dPosOut),
      //  the procedure guarantees positive values
      //  and netto Nion = nPos - nNeg = nOut - nIn
      if (nElectronOut <= 0) nElectronOut = 0;  //< unphysical

      // either the above has not been executed or nPosIonOut was not positive
      if (nElectronOut >= nElectronIn) {
        nNegIonOut = (std::exp(eta * dx) - 1) * nElectronIn;  //< >= 0
        // nNegIonOut = eta / (alpha - eta) * (nElectronOut - nElectronIn);
        nPosIonOut = (nElectronOut - nElectronIn) + nNegIonOut;  //< >= 0
      } else {
        nPosIonOut = (std::exp(alpha * dx) - 1) * nElectronIn;  //< >= 0
        // nPosIonOut = alpha / (alpha - eta) * (nElectronOut - nElectronIn);
        nNegIonOut = nPosIonOut - (nElectronOut - nElectronIn);  //< >= 0
      }
    }
  }
}

// Get mean size of avalanche when going from x to x+dx
void GetMeanAvalancheSizeFromStep(double dx, const long nElectronIn,
                                  const double alpha, const double eta,
                                  long &nElectronOut, double &nPosIonOut,
                                  double &nNegIonOut) {
  // Mean size gain
  nPosIonOut = 0;
  nNegIonOut = 0;
  const double ndx = exp((alpha - eta) * dx);
  // for electrons
  nElectronOut = nElectronIn * ndx;
  if (nElectronOut <= 0) nElectronOut = 0;  //< unphysical

  // either the above has not been executed or nPosIonOut was not positive
  if (nElectronOut >= nElectronIn) {
    nNegIonOut = (std::exp(eta * dx) - 1) * nElectronIn;     //< >= 0
    nPosIonOut = (nElectronOut - nElectronIn) + nNegIonOut;  //< >= 0
  } else {
    nPosIonOut = (std::exp(alpha * dx) - 1) * nElectronIn;   //< >= 0
    nNegIonOut = nPosIonOut - (nElectronOut - nElectronIn);  //< >= 0
  }
}

}  // namespace

namespace Garfield {

// Public:
AvalancheGridSpaceCharge::AvalancheGridSpaceCharge() {
  m_vEElliptic.reserve(20000);
  m_vKElliptic.reserve(20000);
  m_vXElliptic.reserve(20000);
  m_zGrid.reserve(5000);
  m_rGrid.reserve(1000);

  // Import the elliptic integral values
  const std::string path = std::getenv("GARFIELD_INSTALL");
  ImportEllipticIntegralValues(path +
                               "/share/Garfield/Data/elliptic_integrals.txt");
}

void AvalancheGridSpaceCharge::Reset() {
  m_time = 0.;
  m_time0 = 0.;
  m_dt = 0.;
  m_nTotElectron = 0;
  m_nTotPosIons = 0;
  m_run = true;

  m_vCoNGasLayer.resize(0);
  m_vElectrons.resize(0);
  m_vNElectronEvolution.resize(0);
  m_grid.resize(0);
  m_vYPointInGasGap.resize(0);
  m_vIndexGasGaps = {0};
  m_ezBkg = {0};
  m_vSaturatedGaps.resize(0);

  m_bDriftAvalanche = false;
  m_bImportAvalanche = false;
  m_bPreparedImportAvalanche = false;
  m_bFieldK = false;

  std::cout << m_className << "::Reset: Instance reset, ready to use again.\n";
}

void AvalancheGridSpaceCharge::Set2dGrid(const double zmin, const double zmax,
                                         const int zsteps, const double rmax,
                                         const int rsteps) {
  m_isgridset = true;

  if (zmin >= zmax || zsteps <= 0 || 0 >= rmax || rsteps <= 0) {
    std::cerr << m_className
              << "::Set2dGrid: Error. Grid is not properly defined.\n";
    return;
  }

  // set z grid
  m_zSteps = zsteps;
  m_zStepSize = (zmax - zmin) / zsteps;
  // m_zGrid.resize(zsteps + 1);
  for (int i = 0; i < zsteps + 1; i++) {  //< put one more to include the last point on the grid
    m_zGrid.push_back(zmin + i * m_zStepSize);
  }

  // set r grid
  m_rSteps = rsteps;
  m_rStepSize = rmax / rsteps;
  // m_rGrid.resize(rsteps + 1);
  for (int i = 0; i < rsteps + 1; i++) {  //< put one more to include the last point on the grid
    m_rGrid.push_back(0 + i * m_rStepSize);
  }

  if (m_bDebug) {
    std::cout << m_className << "::Set2dGrid: Grid created:\n"
              << "       z range = (" << zmin << "," << zmax << ").\n"
              << "       r range = (" << 0 << "," << rmax << ").\n";
  }
}

void AvalancheGridSpaceCharge::ImportElectronsFromAvalancheMicroscopic(
    Garfield::AvalancheMicroscopic *avmc) {
  if (!avmc) return;

  if (!m_bImportAvalanche) {
    m_bImportAvalanche = true;

    // resize the electrons according to # gap
    if (m_ParallelPlate) {
      m_ParallelPlate->IndexOfGasGaps(m_vIndexGasGaps);
    }
    m_vElectrons.resize(m_vIndexGasGaps.size());
  }

  for (const auto &electron : avmc->GetElectrons()) {
    // StatusOutsideTimeWindow (if electrons has been transported until reaching
    // status -17 ~ still active)
    if (electron.status != -17) {
      if (m_bDebug)
        std::cerr << m_className
                  << "::ImportElectronsFromAvalancheMicroscopic: Status is not "
                     "-17, continue.\n";
      continue;
    }
    int k = 0;
    if (m_ParallelPlate) {
      int ind;
      double eps;
      if (!m_ParallelPlate->getLayer(electron.path.back().y, ind, eps)) {
        std::cerr << m_className
                  << "::ImportElectronsFromAvalancheMicroscopic: Electron "
                     "outside component.\n";
        continue;
      }
      k = GetGasGapNumber(ind);
      if (k == -1) {
        std::cerr << m_className
                  << "::ImportElectronsFromAvalancheMicroscopic: Electron is "
                     "not in a gas gap, continue.\n";
        continue;
      }
    }
    // add electrons to a vector-list
    Point pt{};
    pt.x = electron.path.back().x;
    pt.y = electron.path.back().y;
    pt.z = electron.path.back().z;
    pt.t = electron.path.back().t;
    m_vElectrons[k].push_back(std::move(pt));

    if (m_bDebug)
      std::cout
          << m_className
          << "::ImportElectronsFromAvalancheMicroscopic: Electron added, y: "
          << electron.path.back().y << " and gas gap: " << k + 1 << "\n";
  }
}

void AvalancheGridSpaceCharge::AvalancheElectron(const double x, const double y,
                                                 const double z, const double t,
                                                 const int n) {
  int gasGap = 0;
  // check if avalanche electron in a gas gap
  if (m_ParallelPlate) {
    int ind;
    double eps = -1;
    if (!m_ParallelPlate->getLayer(y, ind, eps) && eps != 1.) {
      std::cerr << m_className
                << "AvalancheElectron: Electron is not in a gas gap.";
      return;
    }
    // determine indices of gas gaps
    m_ParallelPlate->IndexOfGasGaps(m_vIndexGasGaps);
    gasGap = GetGasGapNumber(ind);
  } else {
    // put the y-electron coord as reference
    m_vYPointInGasGap.push_back(y);
  }

  if (!m_bDriftAvalanche) {
    m_bDriftAvalanche = true;
  }

  if (m_time == 0 && m_time != t && m_bDebug)
    std::cerr
        << m_className
        << "::AvalancheElectron: Overwriting start time of avalanche for t "
           "= 0 to "
        << t << ".\n";

  m_time = t;
  m_time0 = t;

  // prepare the CoN
  for (int i = 0; i < (int)m_vIndexGasGaps.size(); i++) {
    if (i != gasGap) {
      // HS: not sure it's a good idea to initialize with NAN...
      m_vCoNGasLayer.push_back({NAN, NAN, NAN});
    } else if (i == gasGap) {
      m_vCoNGasLayer.push_back({x, y, z});
    }
  }

  if (m_vCoNGasLayer.size() == 0) {
    std::cerr << m_className
              << "::AvalancheElectron: Could not determine center.\n";
  }
  Prepare2dMesh();
  // HS: check!!
  if (SnapTo2dGrid(x, y, z, n, gasGap) && m_bDebug)
    std::cerr << m_className
              << "::AvalancheElectron: Electron added at (t, x, y, z) =  (" << t
              << ", " << x << ", " << y << ", " << z << ").\n";
}

void AvalancheGridSpaceCharge::AddExtraAvalancheElectron(double y, int n) {
  if (!m_bDriftAvalanche) {
    std::cerr << m_className
              << "::AddExtraAvalancheElectron: First use AvalancheElectron.\n";
    return;
  }

  // check if avalanche electron in a gas gap
  int gasGap = 0;
  if (m_ParallelPlate) {
    int ind;
    double eps = -1;
    if (!m_ParallelPlate->getLayer(y, ind, eps) && eps != 1.) {
      std::cerr << m_className
                << "AddExtraAvalancheElectron: Electron is not in a gas gap.";
      return;
    }
    gasGap = GetGasGapNumber(ind);
  }

  // nothing yet in this gas gap
  if (std::isnan(m_vCoNGasLayer[gasGap][0])) {
    m_vCoNGasLayer[gasGap] = {0, y, 0};
  }
  // HS: check!!
  if (SnapTo2dGrid(m_vCoNGasLayer[gasGap][0], y, m_vCoNGasLayer[gasGap][2], n,
                   gasGap) &&
      m_bDebug)
    std::cout << m_className << "::AddExtraAvalancheElectron: "
              << "Electron added at (t, x, y, z) =  (" << m_time << ", "
              << m_vCoNGasLayer[gasGap][0] << ", " << y << ", "
              << m_vCoNGasLayer[gasGap][2] << ").\n";
}

void AvalancheGridSpaceCharge::StartGridAvalanche(double dtime) {
  // avalanche the electrons until a certain delta-time OR there are no
  // electrons left in the gap
  if ((!m_bImportAvalanche && !m_bDriftAvalanche) || !m_sensor) return;

  // prepare the imported avalanche
  if (m_bImportAvalanche && !m_bPreparedImportAvalanche) {
    PrepareElectronsFromMicroscopicAvalanche();
    if (m_bDebug)
      std::cerr << m_className
                << "::StartGridAvalanche: Microscopic electrons successfully "
                   "prepared.\n";
    m_bPreparedImportAvalanche = true;
  }

  // check if electrons are on grid
  if (m_nTotElectron <= 0) {
    std::cerr << m_className << "::StartGridAvalanche: Cancelled "
              << m_nTotElectron << " electrons on grid.\n";
    return;
  }

  if (dtime == -1) {
    // transport until no electrons in gap or an error is reached
    while (true) {
      // Transport the nodes (returns false if 0 electrons in gap)
      if (!TransportTimeStep()) {
        break;
      }
    }
  } else {
    double tStart = m_time;

    while (m_time + m_dt - tStart < dtime) {
      // Transport the nodes (returns false if 0 electrons in gap)
      if (!TransportTimeStep()) {
        break;
      }
    }
  }

  if (!m_vNElectronEvolution.empty()) {
    // determine maximal size of electron maxSize at time maxTime.
    auto maxSize = std::max_element(m_vNElectronEvolution.begin(),
                                    m_vNElectronEvolution.end(),
                                    [](const std::pair<double, long> &p1,
                                       const std::pair<double, long> &p2) {
                                      return p1.second < p2.second;
                                    });
    double maxTime = m_time0 + m_dt *
            ((double)std::distance(m_vNElectronEvolution.begin(), maxSize) + 1);

    std::cout << m_className
              << "::StartGridAvalanche: Avalanche maximum size of "
              << maxSize->second << " electrons reached at " << maxTime
              << " ns.\n";

    std::cout << m_className
              << "::StartGridAvalanche: Final avalanche size (produced "
                 "positive charge) = "
              << m_nTotPosIons << " ended at t = " << m_time << " ns.\n";
  }
}

void AvalancheGridSpaceCharge::ExportGrid(const std::string &filename) {


  std::ofstream exportElectrons(filename + "_electrons.csv");
  if (!exportElectrons.is_open()) {
    std::cerr << "Error opening e- file.\n";
    return;
  }
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      exportElectrons << m_grid[iz][ir].nElectron << " ";
    }
    exportElectrons << "\n";
  }
  exportElectrons.close();

  std::ofstream exportPosIon(filename + "_posion.csv");
  if (!exportPosIon.is_open()) {
    std::cerr << "Error opening p+ file.\n";
    return;
  }
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      exportPosIon << std::floor(m_grid[iz][ir].nPosIon) << " ";
    }
    exportPosIon << "\n";
  }
  exportPosIon.close();

  std::ofstream exportNegIon(filename + "_negion.csv");
  if (!exportNegIon.is_open()) {
    std::cerr << "Error opening n- file.\n";
    return;
  }
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      exportNegIon << std::floor(m_grid[iz][ir].nNegIon) << " ";
    }
    exportNegIon << "\n";
  }
  exportNegIon.close();

  std::ofstream exportZField(filename + "_eFieldZ.csv");
  if (!exportZField.is_open()) {
    std::cerr << "Error opening E_z file.\n";
    return;
  }
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      GridNode *nd = &m_grid[iz][ir];
      int gasGap = nd->gasGapIndex;
      double EField = nd->eFieldZ + m_ezBkg[gasGap];
      exportZField << std::round(EField) << " ";
    }
    exportZField << "\n";
  }
  exportZField.close();

  std::ofstream exportRField(filename + "_eFieldR.csv");
  if (!exportRField.is_open()) {
    std::cerr << "Error opening E_r file.\n";
    return;
  }
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      double EField = m_grid[iz][ir].eFieldR;
      exportRField << std::round(EField) << " ";
    }
    exportRField << "\n";
  }
  exportRField.close();

  std::ofstream exportMagField(filename + "_MagField.csv");
  if (!exportMagField.is_open()) {
    std::cerr << "Error opening E_r file.\n";
    return;
  }
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      GridNode *nd = &m_grid[iz][ir];
      int gasGap = nd->gasGapIndex;
      double EField = Mag(nd->eFieldZ + m_ezBkg[gasGap], nd->eFieldR);
      exportMagField << std::round(EField) << " ";
    }
    exportMagField << "\n";
  }
  exportMagField.close();

  if (m_bDebug) {
    std::cout << m_className << "::ExportGrid: Grids exported.\n";
  }
}

void AvalancheGridSpaceCharge::ImportEllipticIntegralValues(
    const std::string &filename) {
  // reads values of the elliptic functions
  // it has a very special form to find the values rapidly (not given for a
  // completely non uniform grid) from x = 0 to 10 it is in steps of 1e-3. From
  // 10 to 1e4 in steps of 1. Then in steps of 1000 until 1e7.

  m_vXElliptic.resize(0);
  m_vKElliptic.resize(0);
  m_vEElliptic.resize(0);

  std::ifstream ellipticStream(filename);

  if (!ellipticStream) {
    std::cerr << m_className
              << "::ImportEllipticIntegralValues: Could not open file.\n";
  }

  for (std::string line; std::getline(ellipticStream, line);) {
    std::istringstream iss(line);
    double value = 0.;
    iss >> value;
    m_vXElliptic.push_back(value);
    iss >> value;
    m_vKElliptic.push_back(value);
    iss >> value;
    m_vEElliptic.push_back(value);
  }

  ellipticStream.close();
  m_bImportElliptic = true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Private Section:
bool AvalancheGridSpaceCharge::SnapTo2dGrid(const double x, const double y,
                                            const double z, const long n,
                                            const int gasLayer) {
  // Snap electron from AvalancheMicroscopic to the predefined grid
  if (!m_isgridset) {
    std::cerr << m_className << "::SnapTo2dGrid: Grid is not defined.\n";
    return false;
  }

  // HS: why make a copy?
  auto CoN = m_vCoNGasLayer[gasLayer];
  // y in micro is z in grid space-charge
  double r = std::sqrt((x - CoN[0]) * (x - CoN[0]) + (z - CoN[2]) * (z - CoN[2]));
  int iZ = (int)std::round((y - m_zGrid.front()) / m_zStepSize);
  int iR = (int)std::round(r / m_rStepSize);

  if (m_bDebug) {
    std::cout << m_className << "::SnapTo2dGrid: iz = " << iZ
              << ", ir = " << iR << ".\n";
  }

  if (iZ < 0 || iZ > m_zSteps || iR < 0 || iR > m_rSteps) {
    if (m_bDebug) {
      std::cerr << m_className
                << "::SnapTo2dGrid: Point is outside the grid.\n";
    }
    return false;
  }

  // Add point to the grid.
  m_grid[iZ][iR].time = m_time;

  // When snapping the electron to the grid the distance traveled can yield
  // additional electrons or get attached. (depends on if against E field or
  // along ...). e-field is along y (micro)
  double step = m_zGrid[iZ] - y;
  // determine if against (ok) or with e field (not ok):
  int against = (step > 0 && m_ezBkg[gasLayer] < 0) ||
                (step < 0 && m_ezBkg[gasLayer] > 0);

  // sanity check
  if (m_grid[iZ][iR].gasGapIndex != gasLayer) {
    std::cerr << m_className
              << "::SnapTo2dGrid: Gas layer index does not match.\n";
    return false;
  }

  if (!against) {
    m_grid[iZ][iR].nElectron += n;
    m_nTotElectron += n;
    if (m_bDebug)
      std::cerr << m_className
                << "::SnapTo2dGrid: snap along e-field, continue.\n";
    return true;
  }

  // make step positive
  long nEOut;
  double nPosOut, nNegOut;
  GetAvalancheSizeFromStep(
      std::abs(step), n, m_grid[iZ][iR].townsend,
      m_grid[iZ][iR].attachment, nEOut, nPosOut, nNegOut);
  if (nEOut == 0) {
    if (m_bDebug)
      std::cerr << m_className << "::SnapTo2dGrid: e- from " << n
                << " to 0 -> cancel.\n";
    return false;
  }

  m_grid[iZ][iR].nElectron += nEOut;
  m_grid[iZ][iR].nPosIon += nPosOut;
  m_grid[iZ][iR].nNegIon += nNegOut;
  m_nTotElectron += nEOut;
  m_nTotPosIons += (long)nPosOut;

  if (m_bDebug) {
    std::cout << m_className << "::SnapTo2dGrid: e- from " << n << " to "
              << nEOut << " p+: " << nPosOut << " n-: " << nNegOut << ".\n"
              << "    Snapped to (z, r) = (" << y << " -> "
              << m_zGrid[iZ] << ", " << r << " -> "
              << m_rGrid[iR] << ").\n";
  }
  return true;
}

void AvalancheGridSpaceCharge::Prepare2dMesh() {
  // check if sensor is defined
  if (!m_sensor) {
    std::cerr << m_className
              << "::Prepare2dMesh: Sensor is not defined. Abort.\n";
  }

  // get a point (Y global coordinate) in each gas gap
  int n = m_vIndexGasGaps.size();
  m_ezBkg.resize(n);

  if (m_ParallelPlate) {
    m_vYPointInGasGap.resize(n);
    for (int iz = 0; iz <= m_zSteps; iz++) {
      // determine layer index
      int layerIndex = 0;
      double eps = 0.;
      m_ParallelPlate->getLayer(m_zGrid[iz], layerIndex, eps);
      // determine gap number from m_iIndexGasGaps and layer index
      int k = GetGasGapNumber(layerIndex);
      if (k != -1 && k < n) m_vYPointInGasGap[k] = m_zGrid[iz];
    }
  }

  std::vector<double> alpha(n), eta(n), drift(n), dSigmaL(n), dSigmaT(n), wv(n),
      wr(n), alphaPT(n), etaPT(n);
  double e[3], v;
  int status;
  Medium *m = nullptr;
  // iterate through the gas gaps
  for (int k = 0; k < n; k++) {
    m_sensor->ElectricField(0, m_vYPointInGasGap[k], 0, e[0], e[1], e[2], v, m,
                            status);

    if (status != 0) {
      std::cerr
          << m_className
          << "::Prepare2dMesh: Cannot estimate background field for gas gap "
          << k + 1 << ".\n";
    }

    // one expects (ComponentParallelPlate) that the electric field is pointing
    // along y-axis
    //  i.e. Z-axis in our coordinate system.
    m_ezBkg[k] = e[1];
    GetSwarmParameters(std::abs(e[1]), alpha[k], eta[k], drift[k], dSigmaL[k],
                       dSigmaT[k], wv[k], wr[k], alphaPT[k], etaPT[k], k);

    // print-out to double-check the swarm parameters
    std::cout << m_className << "::Prepare2dMesh:\n"
              << "  Gas gap " << k + 1 << "\n"
              << "     Ez: " << m_ezBkg[k] << " (V/cm)\n"
              << "     alphaSST: " << alpha[k] << " (1/cm)\n"
              << "     alphaPT:  " << alphaPT[k] << " (1/cm)\n"
              << "     etaSST: " << eta[k] << " (1/cm)\n"
              << "     etaPT:  " << etaPT[k] << " (1/cm)\n"
              << "     drift velocity (Wv): " << drift[k] << " (cm/ns)\n"
              << "     Wr (!= Wv): " << wr[k] << " (cm/ns).\n";
  }

  // Set up mesh
  m_grid.resize(m_zSteps + 1);
  m_zGasGapBoundaries.resize(n);
  for (int iz = 0; iz <= m_zSteps; iz++) {
    m_grid[iz].resize(m_rSteps + 1);
    int k = 0;
    // Determine layer index.
    int layerIndex = 0;
    if (m_ParallelPlate) {
      double eps = 0.;
      m_ParallelPlate->getLayer(m_zGrid[iz], layerIndex, eps);
      // Determine gap number.
      k = GetGasGapNumber(layerIndex);
    }
    if (k != -1) {
      // store z index for gap k
      m_zGasGapBoundaries[k].push_back(iz);
    }
    for (int ir = 0; ir <= m_rSteps; ir++) {
      // Set layer index
      m_grid[iz][ir].layerIndex = layerIndex;
      m_grid[iz][ir].gasGapIndex = k;
      // Continue if nodes are not in gas (eps != 1 or k == -1)
      if (k == -1) {
        m_grid[iz][ir].isGasGap = false;
        continue;
      }
      // set swarm parameters & time
      m_grid[iz][ir].townsend = alpha[k];
      m_grid[iz][ir].attachment = eta[k];
      m_grid[iz][ir].velocity = drift[k];  //< magnitude! direction against E_field. E/|E|
      m_grid[iz][ir].dSigmaL = dSigmaL[k];
      m_grid[iz][ir].dSigmaT = dSigmaT[k];
      m_grid[iz][ir].Wv = wv[k];
      m_grid[iz][ir].Wr = wr[k];
      m_grid[iz][ir].townsendPT = alphaPT[k];
      m_grid[iz][ir].attachmentPT = etaPT[k];
      m_grid[iz][ir].time = m_time;
    }
  }

  // set anode in each gas gap
  for (int k = 0; k < n; k++) {
    int izMin = m_zGasGapBoundaries[k].front();
    int izMax = m_zGasGapBoundaries[k].back();
    int izAnode = (m_ezBkg[k] > 0) ? izMin : izMax;
    for (int ir = 0; ir <= m_rSteps; ir++) {
      m_grid[izAnode][ir].anode = true;
    }
  }

  // we define the time step as dz / max(Wr)
  m_dt = m_zStepSize / *std::max_element(wr.begin(), wr.end());

  if (m_bDebug) {
    std::cout << m_className << "::Prepare2dMesh: Time step per loop: "
              << m_dt << " ns.\n";
  }
}

void AvalancheGridSpaceCharge::PrepareElectronsFromMicroscopicAvalanche() {
  double tMicro = 0;
  long neTotal = 0;
  for (int k = 0; k < (int)m_vIndexGasGaps.size(); k++) {
    // calculate middle coord of electron cloud and add all electrons to the
    // grid/mesh per gas gap.
    auto neGasGap = m_vElectrons[k].size();

    // continue if no electron in the gas gap
    if (neGasGap <= 0) {
      // HS: not sure it's a good idea to initialize with NAN.
      m_vCoNGasLayer.push_back({NAN, NAN, NAN});
      continue;
    }
    neTotal += neGasGap;
    double xMicro = 0, yMicro = 0, zMicro = 0;

    // set time and center of electron number
    for (const auto &electron : m_vElectrons[k]) {
      xMicro += electron.x;
      yMicro += electron.y;
      zMicro += electron.z;
      tMicro += electron.t;
    }
    m_vCoNGasLayer.push_back({xMicro / (double)neGasGap,
                              yMicro / (double)neGasGap,
                              zMicro / (double)neGasGap});

    if (m_bDebug) {
      std::cout
          << m_className
          << "::PrepareElectronsFromMicroscopicAvalanche: mean center gas gap "
          << k + 1 << " X = (" << m_vCoNGasLayer[k][0] << ", "
          << m_vCoNGasLayer[k][1] << ", " << m_vCoNGasLayer[k][2] << ").\n";
    }
  }

  // set the times
  m_time = tMicro / (double)neTotal;
  m_time0 = tMicro / (double)neTotal;

  // prepare the mesh (needs m_vCoN)
  Prepare2dMesh();

  // place all electrons onto the grid
  for (int k = 0; k < (int)m_vIndexGasGaps.size(); k++)
    for (auto &electron : m_vElectrons[k]) {  //< does not proceed if no
                                              //electrons in respective gap
      // HS: check!!
      if (SnapTo2dGrid(electron.x, electron.y, electron.z, 1, k) && m_bDebug) {
        std::cout << m_className
                  << "::PrepareElectronsFromMicroscopicAvalanche: "
		  << "Electron added in gas gap " << k + 1 << "\n"
                  << "      at (x,y,z) =  (" << electron.x << "," << electron.y
                  << "," << electron.z << ").\n";
      }
    }
}

void AvalancheGridSpaceCharge::GetSwarmParameters(
    const double MagEField, double &alpha, double &eta, double &drift,
    double &dSigmaL, double &dSigmaT, double &wv, double &wr, double &alphaPT,
    double &etaPT, int gasGap) {
  if (m_bDebug && false)
    std::cerr << m_className
              << "::GetSwarmParameters: Getting parameters at "
                 "|E| = "
              << MagEField << ".\n";

  // medium from sensor
  Medium *m = m_sensor->GetMedium(0, m_vYPointInGasGap[gasGap], 0);

  // alpha
  m->ElectronTownsend(0., MagEField, 0., 0., 0., 0., alpha);
  // eta
  m->ElectronAttachment(0., MagEField, 0., 0., 0., 0., eta);

  // mag of velocity
  double vx, vy, vz;
  m->ElectronVelocity(0., MagEField, 0., 0., 0., 0., vx, vy, vz);
  drift = std::sqrt(vx * vx + vy * vy + vz * vz);  //< Wv in Magboltz
  // wv, wr <- take only Wr and 'drift' for Wv
  wr = 0.;
  if (!m_bUseTOF ||
      !m->ElectronVelocityFluxBulk(0., MagEField, 0., 0., 0., 0., wv, wr) ||
      wr < Small) {
    m_bWrAvailable = false;
    wr = drift;
  }
  // rates, if not available we take (alpha-eta)SST and ratio of alpha/eta =
  // Rion/Ratt Rion-Ratt = Reff (tagashira eq.)
  //  -> Rion converged to a rate with Wr and DL (using the one equation),
  //  alphaSST is either from SST and if not converged from magboltz itself.
  double rion = 0, ratt = 0;
  if (!m_bUseTOF ||
      !m->ElectronTOFIonisation(0., MagEField, 0., 0., 0., 0., rion) ||
      !m->ElectronTOFAttachment(0., MagEField, 0., 0., 0., 0., ratt)) {
    if (m_bDebug) {
      std::cerr << m_className
                << "::GetSwarmParameters: TOF Rates not available.\n";
    }

    m_bRatesAvailable = false;

    // Diffusionless approximation
    rion = alpha * wr;
    ratt = eta * wr;
  }
  // calculate alpha/eta PT
  alphaPT = rion / wr;
  etaPT = ratt / wr;

  // diffusion coefficients
  m->ElectronDiffusion(0., MagEField, 0., 0., 0., 0., dSigmaL, dSigmaT);

  // print (and information about units!)
  if (m_bDebug && false) {
    std::cout << m_className << "::GetSwarmParameters:\n"
              << "  Townsend = " << alpha << " [1/cm], Attachment = " << eta
              << " [1/cm], Flux Velocity = "
              << drift
              //                      << " [cm/ns], Wv = " << wv
              << " [cm/ns], Bulk velocity = " << wr << " [cm/ns].\n"
              << "  TOF Ionization rate = " << alphaPT
              << " [1/ns], TOF Attachment rate = " << etaPT << "[1/ns].\n"
              << "  Longitudinal Diffusion = " << dSigmaL << " [sqrt(cm)],"
              << " Transversal Diffusion = " << dSigmaT << " [sqrt(cm)].\n";
  }
}

bool AvalancheGridSpaceCharge::TransportTimeStep() {
  // Transport GridMesh one time-step with updated E-fields (Lippmann et al.
  // approach)
  if (!m_run) return false;

  // local helper variables
  long nElectronOut;
  double nPosIonOut, nNegIonOut;

  if (m_bDebug) {
    std::cout << m_className
              << "::TransportTimeStep: Start time: " << m_time << "\n";
  }

  // choose MC or Mean version depending on m_bMC; total electron > 1e5
  std::function<void(double, const long, const double, const double, long &,
                     double &, double &)>
      AvalancheGain = GetAvalancheSizeFromStep;
  if (!m_bMC && m_nTotElectron > 1e5) {
    AvalancheGain = GetMeanAvalancheSizeFromStep;
    // TODO: Diffusion
  }

  // update the nodes for the next run (SC-field and swarm parameters)
  // MRPC: SC-effect only within each gas gap and option="coulomb"
  if (m_bSpaceCharge && m_nTotElectron > 1e5) {
    for (int iz = 0; iz <= m_zSteps; iz++) {
      // continue if not in gas gap
      int gasGap = m_grid[iz][0].gasGapIndex;
      if (gasGap == -1) continue;
      for (int ir = 0; ir <= m_rSteps; ir++) {
        GridNode *nd = &m_grid[iz][ir];

        // reset local fields at node
        nd->eFieldZ = 0;
        nd->eFieldR = 0;

        // continue if: no electrons, an anode
        if ((double)nd->nElectron < 0.5 || nd->anode) continue;

        // update space charge field
        // calculate field at the current bin from all other bins containing
        // charge
        GetLocalField(iz, ir, nd->eFieldZ, nd->eFieldR, m_sFieldOption, gasGap);

        // check if local field reaches background field values.
        double MagEField = Mag(nd->eFieldZ + m_ezBkg[gasGap], nd->eFieldR);
        if (MagEField - std::abs(m_ezBkg[gasGap]) >=
                m_fStreamerK * std::abs(m_ezBkg[gasGap]) &&
            !m_bFieldK) {
          std::cout << m_className << ":TransportTimeStep:\n"
                    << "    Space-charge field reached "
                    << std::to_string(int(m_fStreamerK * 100))
                    << "% of background field in gas gap " << gasGap + 1
                    << "\n";
          // TODO: Total electrons in gas gap "gasGap"
          m_lElectronsK = m_vNElectronEvolution.back().second;
          m_bFieldK = true;
        }

        // calculate the swarm parameters
        // HS: why calculate MagEField again?
        MagEField = Mag(nd->eFieldZ + m_ezBkg[gasGap], nd->eFieldR);
        GetSwarmParameters(MagEField, nd->townsend, nd->attachment,
                           nd->velocity, nd->dSigmaL, nd->dSigmaT, nd->Wv,
                           nd->Wr, nd->townsendPT, nd->attachmentPT, gasGap);

        // get new step distance
        double step = std::abs(nd->Wr * m_dt);

        // adaptive time stepping (this routine takes the smallest dt needed for
        // the current sc-field)
        if (m_bAdaptiveTime && step >= 2. * m_zStepSize) {
          // reset dt
          double dtPrev = m_dt;
          m_dt = m_zStepSize / nd->Wr;

          if (m_bDebug) {
            std::cout << m_className << "::TransportTimeStep: Changed dt from "
                      << dtPrev << " to: " << m_dt << "\n"
                      << "      due to step size: " << step
                      << " bulk velocity: " << nd->Wr << "\n"
                      << "      electric field: " << MagEField
                      << " alpha: " << nd->townsendPT
                      << " eta: " << nd->attachmentPT << "\n"
                      << "      diffusion longitudinal/transversal: "
                      << nd->dSigmaL << " " << nd->dSigmaT << "\n";
            ExportGrid("TIME_STEP_ADAPTION_" + std::to_string(m_dt));
          }
        }
      }
    }
  }

  // finish if stop reached and set at 100 * K %
  if (m_bStopAtK && m_bFieldK) return false;

  // transport the electrons with the update sc-field, swarm parameter and dt
  for (int iz = 0; iz <= m_zSteps; iz++) {
    // continue if not in gas gap
    int gasGap = m_grid[iz][0].gasGapIndex;
    if (gasGap == -1) continue;

    for (int ir = 0; ir <= m_rSteps; ir++) {
      GridNode *nd = &m_grid[iz][ir];

      // continue if at anode or no electrons
      // HS: why static_cast<double>(nd->nElectron) < 0.5 instead of
      //     just nd->nElectron < 1?
      if (nd->anode || static_cast<double>(nd->nElectron) < 0.5) continue;

      // update step distance
      double step = std::abs(nd->Wr * m_dt);

      // calculate new avalanche size at X + step
      // HS: use a vector<bool> to keep track of which gaps are saturated?
      if (!m_bSpaceCharge &&
          (std::find(m_vSaturatedGaps.begin(), m_vSaturatedGaps.end(),
                     gasGap) != m_vSaturatedGaps.end()
               ? true
               : false)) {
        // Saturated case, don't evolve electrons in size
        nElectronOut = nd->nElectron;  
        nPosIonOut = 0,
        nNegIonOut = 0;  //< strictly this is completely wrong because
                         //SC-bremsung creates huge amounts of ions
      } else {
        AvalancheGain(step, nd->nElectron, nd->townsendPT, nd->attachmentPT,
                      nElectronOut, nPosIonOut, nNegIonOut);
      }
      m_nTotPosIons += std::round(nPosIonOut);

      // calculate steps against electric field i.e. correct sign.
      double MagEField = Mag(nd->eFieldZ + m_ezBkg[gasGap], nd->eFieldR);
      double stepZ = step * (-(nd->eFieldZ + m_ezBkg[gasGap]) / MagEField);
      double stepR = step * (-(nd->eFieldR) / MagEField);

      if (m_bDiffusion) {
        // correct the stepping from diffusion + charge distribution
        DiffuseTimeStep(step, nElectronOut, std::round(nPosIonOut),
                        std::round(nNegIonOut), iz, ir, gasGap);
      } else {
        // calculate steps and distribute charges (no diffusion)
        DistributeCharges(nElectronOut, std::round(nPosIonOut),
                          std::round(nNegIonOut), iz, ir, stepZ, stepR, gasGap);
      }

      if (m_sensor->GetNumberOfElectrodes() > 0) {
        // adding the movement to signal: (go to global cartesian coordinates
        // using phi = 0)
        // TODO: discretize phi and add signal for each with nElectron / M on
        // each element.
        // TODO: at anode the signal from diffusion is due to bounded plane not
        // netto 0 because diffusion
        //  tends more backwards, since forward they reach at earlier distance
        //  the boundary
        // HS: why? This means we have to compute cos(0), sin(0) every time.
        double x0, y0, z0;
        GetGlobalCoordinates(m_rGrid[ir], m_zGrid[iz], 0., x0, y0, z0, gasGap);

        // z-step outside gasGap domain, resize to stepZ = Anode - Current
        int izMin = m_zGasGapBoundaries[gasGap].front();
        int izMax = m_zGasGapBoundaries[gasGap].back();
        if (m_zGrid[iz] + stepZ < m_zGrid[izMin]) {
          stepZ = (m_zGrid[izMin] - m_zGrid[iz]);
        } else if (m_zGrid[iz] + stepZ > m_zGrid[izMax]) {
          stepZ = (m_zGrid[izMax] - m_zGrid[iz]);
        }

        // Induced current from flux drift velocity i.e. introduce weight factor
        double weight =
            nd->velocity / nd->Wr;  //< 1 if (Wv = velocity): Wr = flux
        double x1, y1, z1;
        GetGlobalCoordinates(m_rGrid[ir] + stepR, m_zGrid[iz] + stepZ, 0.,
                             x1, y1, z1, gasGap);
        m_sensor->AddSignalWeightingPotential(
            -weight, {nd->time, nd->time + m_dt},
            {{x0, y0, z0}, {x1, y1, z1}},
            {(double)nd->nElectron, (double)nElectronOut});
      }
    }
  }

  // propagate Grid time (after above for loops due to the adaptive time
  // stepping)
  m_time += m_dt;

  // sums electrons left in gap (count per gap)
  std::vector<long> eOnGrid(m_vIndexGasGaps.size(), 0);

  // update nodes with transported electrons
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      // get node
      // HS: why a pointer?
      GridNode *nd = &m_grid[iz][ir];
      int gasGap = nd->gasGapIndex;
      if (nd->anode && m_bStick) {
        // sticky anode: electron stay and holder electrons add it up
        nd->nElectron += nd->nElectronHolder;
      } else {
        // update node with electrons from holder
        nd->nElectron = nd->nElectronHolder;
      }
      // ions add up, also at the anode
      nd->nPosIon += nd->nPosIonHolder;
      nd->nNegIon += nd->nNegIonHolder;

      // reset node Holder
      nd->nElectronHolder = 0;
      nd->nPosIonHolder = 0;
      nd->nNegIonHolder = 0;

      // add electrons if they are not stuck
      if (!(nd->anode && m_bStick)) eOnGrid[gasGap] += nd->nElectron;

      // move node in time (even if no electrons in there)
      nd->time += m_dt;
    }
  }
  // add total electrons in gap to grid and to evolution vector
  m_nTotElectron = std::accumulate(eOnGrid.begin(), eOnGrid.end(), 0.);
  m_vNElectronEvolution.push_back(std::make_pair(m_time, m_nTotElectron));
  // continue run if electrons left on grid
  m_run = m_nTotElectron > 0;

  // determine saturated gaps at each time step
  // clear: anode-absorption activates avalanche to grow again
  m_vSaturatedGaps.resize(0);  
  for (int k = 0; k < (int)eOnGrid.size(); k++) {
    if (!m_bSpaceCharge && eOnGrid[k] > m_lNCrit) {
      m_vSaturatedGaps.push_back(k);
    }
    if (m_bDebug) {
      std::cout << m_className
                << "::TransportTimeStep: Electrons active on grid in gas gap "
                << k + 1 << ": " << eOnGrid[k] << "\n";
    }
  }

  return true;
}

void AvalancheGridSpaceCharge::DiffuseTimeStep(double dx, long nElectron,
                                               double nPosIon, double nNegIon,
                                               int iz, int ir, int gasGap) {
  // add diffusion onto the step dx
  long rest = 0, groups = 0, groupSize = 0;
  double sqrtdx = std::sqrt(dx);
  double r = m_rGrid[ir];

  GridNode *nd = &m_grid[iz][ir];

  // Diffuse in Groups of minimum m_dMinGroups groups a size groupSize:
  for (auto &size : m_vGroupSizes) {
    if (nElectron > m_dMinGroups * size) {
      groupSize = size;                            //< size of a single group
      groups = std::floor(nElectron / groupSize);  //< real # groups
      rest = nElectron - groups * groupSize;       //< rest
      break;
    }
  }

  // if electron are less than m_dMinGroups * GroupSize(-1) then we diffuse each
  // electron by itself:
  if (nElectron <= m_dMinGroups * m_vGroupSizes.back()) {
    groupSize = 1;
    groups = nElectron;
    rest = 0;
  }

  // calculate diffusion and add to transport step
  double sinTheta = 0.;
  double cosTheta = 1.;
  double MagEField = Mag(nd->eFieldZ + m_ezBkg[gasGap], nd->eFieldR);
  if (MagEField > 1.e-8) {
    cosTheta = (-(nd->eFieldZ + m_ezBkg[gasGap]) / MagEField);
    sinTheta = (-(nd->eFieldR) / MagEField);
  }

  for (int group = 0; group < groups; group++) {
    // in the last loop we add the rest to the groupSize.
    if (group == groups - 1) groupSize += rest;
    // diffuse each group as if it is a particle.
    // (U,V,W) Local coord system along E field.
    //<W is along E field, V is along e_phi and U perpendicular V and W
    const double dU = RndmGaussian(0, nd->dSigmaT * sqrtdx);
    const double dV = RndmGaussian(0, nd->dSigmaT * sqrtdx);
    const double dW = RndmGaussian(
        dx, nd->dSigmaL * sqrtdx);  //< along E-field i.e. mean = dx1
    // transform to avalanche coordinate system
    // (Z,R,Y) where R mimics an X axis and Y is perpendicular to R and Z
    const double dX = cosTheta * dU + sinTheta * dW;
    // dY = dV
    const double stepZ =
        cosTheta * dW -
        sinTheta * dU;  //< sign seems correct due to sign in cos- and sinTheta
    // calculate the change of radius
    const double stepR = std::sqrt((r + dX) * (r + dX) + dV * dV) -
                         r;  //< sign correct and stepR >= -r

    // distribute nodes and add electrons/ions to Holder
    DistributeCharges(
        groupSize, nPosIon * (double)groupSize / (double)nElectron,
        nNegIon * (double)groupSize / (double)nElectron, iz, ir, stepZ, stepR,
        gasGap);  // < fractional ion number is allowed otherwise loss of ions
  }
}

void AvalancheGridSpaceCharge::DistributeCharges(long nElectron, double nPosIon,
                                                 double nNegIon, int iz, int ir,
                                                 double stepZ, double stepR,
                                                 int gasGap) {
  // distributes the charges from a movement in Z and R direction
  // Compute the effective step size in r-direction
  stepR = std::abs(m_rGrid[ir] + stepR) - m_rGrid[ir];  
  double zRatio = (m_zGrid[iz] + stepZ - m_zGrid[0]) / m_zStepSize;
  double rRatio = std::abs(m_rGrid[ir] + stepR) /
                  m_rStepSize;  // can travel through r=0

  int izPost, irPost, signZstep, signRstep;
  if (stepZ < 0) {
    izPost = (int)ceil(zRatio);
    signZstep = -1;
  } else if (stepZ == 0) {
    izPost = iz;
    signZstep = 0;
  } else {
    izPost = (int)floor(zRatio);
    signZstep = +1;
  }

  if (stepR < 0) {
    irPost = (int)ceil(rRatio);
    if (irPost == 0) {
      // otherwise it travels to ir = -1
      signRstep = +1;
    } else {
      signRstep = -1;
    }
  } else if (stepR == 0) {
    // no movement in r direction
    irPost = ir;
    signRstep = 0;
  } else {
    irPost = (int)floor(rRatio);
    signRstep = +1;
  }

  // 4 point approximation:
  int izPost2 = izPost + signZstep;
  // always: irPost2 >= 0
  int irPost2 = irPost + signRstep;

  const double bz = std::abs(zRatio - (double)izPost);
  const double az = 1 - bz;
  const double br = std::abs(rRatio - (double)irPost);
  const double ar = 1 - br;

  if (az < 0 || az > 1) throw std::runtime_error("az not in range");
  if (bz < 0 || bz > 1) throw std::runtime_error("bz not in range");
  if (ar < 0 || ar > 1) throw std::runtime_error("ar not in range");
  if (br < 0 || br > 1) throw std::runtime_error("br not in range");

  // check if still in grid else place at boundaries (will be absorbed in next
  // step)?
  int izMin = m_zGasGapBoundaries[gasGap].front();
  int izMax = m_zGasGapBoundaries[gasGap].back();
  if (izPost < izMin) izPost = izMin;
  if (izPost2 < izMin) izPost2 = izMin;
  if (izPost > izMax) izPost = izMax;
  if (izPost2 > izMax) izPost2 = izMax;

  if (irPost > m_rSteps) irPost = m_rSteps;
  if (irPost2 > m_rSteps) irPost2 = m_rSteps;

  // add to the nodes the electrons travelled to (into nElectronHolder as they
  // will mix)
  if (nElectron > 200) {
    m_grid[izPost][irPost].nElectronHolder +=
        (long)std::round((double)nElectron * az * ar);
    m_grid[izPost][irPost2].nElectronHolder +=
        (long)std::round((double)nElectron * az * br);
    m_grid[izPost2][irPost].nElectronHolder +=
        (long)std::round((double)nElectron * bz * ar);
    m_grid[izPost2][irPost2].nElectronHolder +=
        (long)std::round((double)nElectron * bz * br);

    // add positive ions to the nodes (smeared values allowed)
    m_grid[izPost][irPost].nPosIonHolder += nPosIon * az * ar;
    m_grid[izPost][irPost2].nPosIonHolder += nPosIon * az * br;
    m_grid[izPost2][irPost].nPosIonHolder += nPosIon * bz * ar;
    m_grid[izPost2][irPost2].nPosIonHolder += nPosIon * bz * br;

    // add negative ions to the nodes (smeared values allowed)
    m_grid[izPost][irPost].nNegIonHolder += nNegIon * az * ar;
    m_grid[izPost][irPost2].nNegIonHolder += nNegIon * az * br;
    m_grid[izPost2][irPost].nNegIonHolder += nNegIon * bz * ar;
    m_grid[izPost2][irPost2].nNegIonHolder += nNegIon * bz * br;

  } else {
    // too large movement of few electrons -> only move it to 1 node (instead of
    // 4)
    izPost = (az >= bz) ? izPost : izPost2;
    irPost = (ar >= br) ? irPost : irPost2;

    m_grid[izPost][irPost].nElectronHolder += nElectron;
    m_grid[izPost][irPost].nPosIonHolder += nPosIon;
    m_grid[izPost][irPost].nNegIonHolder += nNegIon;
  }
}

void AvalancheGridSpaceCharge::GetLocalField(const int iz, const int ir,
                                             double &eFieldZ, double &eFieldR,
                                             const std::string &fieldOption,
                                             int gasGap) {
  // calculate space-charge (local field) at iz/ir
  eFieldZ = 0;
  eFieldR = 0;
  // HS: use enum instead of string.
  if (fieldOption == "coulomb") {
    if (!m_bImportElliptic) {
      throw std::runtime_error("::GetLocalField Elliptic values not imported.");
    }

    // loop over all cells with particles (except itself) and add fields
    for (int fz = 0; fz <= m_zSteps; fz++) {
      // continue if not in gas gap; only add field from charges in same gas gap
      int k = m_grid[fz][0].gasGapIndex;
      if (k == -1 || k != gasGap) continue;
      for (int fr = 0; fr <= m_rSteps; fr++) {
        // add electric field from charge at f at position i
        double N = -m_grid[fz][fr].nElectron + m_grid[fz][fr].nPosIon -
                   m_grid[fz][fr].nNegIon;
        if (std::abs(N) < 1.) continue;  //< N too small to consider
        AddFieldFromChargeAt(iz, ir, fz, fr, N, eFieldZ, eFieldR);
      }
    }
    // add prefactor (final field units V/cm)
    // HS: make the prefactor constexpr
    eFieldZ *= ElementaryCharge / (TwoPi * FourPiEpsilon0);
    eFieldR *= ElementaryCharge / (TwoPi * FourPiEpsilon0);
  } else if (fieldOption == "mirror") {
    // assume symmetric single layer rpc with equal permittivity resistive
    // layers.
    if (!m_bImportElliptic) {
      throw std::runtime_error("::GetLocalField Elliptic values not imported.");
    }

    if (m_vIndexGasGaps.size() > 1) {
      throw std::runtime_error(
          "::GetLocalField Mirror charge option implemented but not tested for "
          "MRPC.");
    }

    // get the rpc (ComponentParallelPlate)
    // HS: why?
    auto *rpc = m_ParallelPlate;

    // loop over all cells with particles and add fields
    int k;
    for (int fz = 0; fz <= m_zSteps; fz++) {
      // continue if not in gas gap; only add field from charges in same gas gap
      k = m_grid[fz][0].gasGapIndex;
      if (k == -1 || k != gasGap) continue;

      // HS: this can be done at initialization time...
      // get epsilon value from neighboring layer (assume both layers have same
      // eps)
      int IndexOfRightLayer = m_vIndexGasGaps[k] + 1;
      // int IndexOfLeftLayer = m_vIndexGasGaps[k] - 1;
      double eps = 1.;  //< neighbored resistive layer thickness from where?
      rpc->getPermittivityFromLayer(IndexOfRightLayer, eps);
      double alpha12 = (1. - eps) / (1. + eps);
      double beta12 = -4. * eps / ((eps + 1.) * (eps + 1) * alpha12);

      // Obtain bounds of current gas gap
      double zTop, zBottom;
      rpc->getZBoundFromLayer(m_vIndexGasGaps[k], zTop, zBottom);

      for (int fr = 0; fr <= m_rSteps; fr++) {
        // charge of interest at f, point of interest at i
        double N = -m_grid[fz][fr].nElectron + m_grid[fz][fr].nPosIon -
                   m_grid[fz][fr].nNegIon;
        if (std::abs(N) < 1.0) continue;  //< N too small to consider

        double zf = m_zGrid[fz];
        double rf = m_rGrid[fr];
        double zi = m_zGrid[iz];
        double ri = m_rGrid[ir];

        // direct charge interaction, delta_Q = 1 (except itself)
        AddFieldFromChargeAt(iz, ir, fz, fr, N, eFieldZ, eFieldR);

        // mirror charge interaction
        for (int i = 0; i < m_iFieldApprox; i++) {
          if (i == 0) {
            // 2a, alpha12 = delta_Q
            double zf0 = zf + 2. * (zTop - zf);
            AddFieldFromChargeAt(iz, ir, zf0, rf, N * alpha12, eFieldZ,
                                 eFieldR);

            // -2a', alpha12 = delta_Q
            zf0 = zf + 2. * (zBottom - zf);
            AddFieldFromChargeAt(iz, ir, zf0, rf, N * alpha12, eFieldZ,
                                 eFieldR);
          } else if (i == 1) {
            // TODO: higher order mirror charges
          } else {
            continue;
          }
        }
      }
    }
    // add prefactor (final field units V/cm)
    // HS: precompute (constexpr)
    eFieldZ *= ElementaryCharge / (TwoPi * FourPiEpsilon0);
    eFieldR *= ElementaryCharge / (TwoPi * FourPiEpsilon0);
  } else if (fieldOption == "relaxation") {
    // TODO: relaxation field method
  } else {
    // default
    eFieldZ = 0;
    eFieldR = 0;
  }
}

bool AvalancheGridSpaceCharge::AddFieldFromChargeAt(int iz, int ir, int fz,
                                                    int fr, double N,
                                                    double &eFieldZ,
                                                    double &eFieldR) {
  // charge of interest at f, point of interest at i
  if (fz == iz and fr == ir) return false;  //< field on itself is not included

  double zi = m_zGrid[iz];
  double ri = m_rGrid[ir];
  double zf = m_zGrid[fz];
  double intermediateEz = 0., intermediateEr = 0.;

  if (fr == 0) {
    // coulomb ball of radius dr / 2
    const double dist = std::sqrt((zi - zf) * (zi - zf) + ri * ri);

    intermediateEz = TwoPi / (dist * dist);
    // HS: do the division by dist in the expression above?
    intermediateEr = intermediateEz * ri / dist;
    intermediateEz *= (zi - zf) / dist;

  } else {  //< rf != 0
    // charged ring
    GetFreeChargedRing(iz, ir, fz, fr, intermediateEz, intermediateEr);
  }
  eFieldZ += intermediateEz * N;
  eFieldR += intermediateEr * N;
  return true;
}

bool AvalancheGridSpaceCharge::AddFieldFromChargeAt(int iz, int ir, double zf,
                                                    double rf, double N,
                                                    double &eFieldZ,
                                                    double &eFieldR) {
  // charge of interest at f, point of interest at i
  double zi = m_zGrid[iz];
  double ri = m_rGrid[ir];
  if (std::abs(zi - zf) / m_zStepSize < 1.e-3 &&
      std::abs(ri - rf) / m_rStepSize < 1.e-3) {
    return false;  //< field on itself is not included
  }
  double intermediateEz = 0, intermediateEr = 0;

  if (std::abs(rf) / m_rStepSize < 0.5) {
    // coulomb ball of radius dr / 2
    const double dist = std::sqrt((zi - zf) * (zi - zf) + ri * ri);
    intermediateEz = TwoPi / (dist * dist);
    intermediateEr = intermediateEz * ri / dist;
    intermediateEz *= (zi - zf) / dist;

  } else {  //< rf != 0
    // charged ring
    GetFreeChargedRing(zi, ri, zf, rf, intermediateEz, intermediateEr);
  }
  eFieldZ += intermediateEz * N;
  eFieldR += intermediateEr * N;
  return true;
}

void AvalancheGridSpaceCharge::GetFreeChargedRing(int iz, int ir, int fz,
                                                  int fr, double &eFieldZ,
                                                  double &eFieldR) {
  // calculate the electric field at point (zi, ri) form charged ring at (zf,
  // rf)

  // precondition
  if (iz == fz && ir == fr) {
    eFieldZ = 0;
    eFieldR = 0;
    return;
  }

  // transform to coordinates and get the field
  double ri = m_rGrid[ir];
  double rf = m_rGrid[fr];
  double zi = m_zGrid[iz];
  double zf = m_zGrid[fz];
  GetFreeChargedRing(zi, ri, zf, rf, eFieldZ, eFieldR);
}

void AvalancheGridSpaceCharge::GetFreeChargedRing(double zi, double ri,
                                                  double zf, double rf,
                                                  double &eFieldZ,
                                                  double &eFieldR) {
  // calculate the electric field at point (zi, ri) form charged ring at (zf,
  // rf)

  // precondition
  if (zi == zf && ri == rf) {
    eFieldZ = 0;
    eFieldR = 0;
    return;
  }

  double dz = zi - zf;  //< I double-checked that's the right sign

  // parameters (see Lippmann Diss.)
  const double a2 = (ri + rf) * (ri + rf) + dz * dz;
  const double b2 = (ri - rf) * (ri - rf) + dz * dz;
  const double b = std::sqrt(b2);
  const double c2 = ri * ri - rf * rf - dz * dz;
  // parameter for elliptic integrals
  const double x =
      -4 * ri * rf / b2;  //< x < 0, i.e. never near x = 1 (singularity)

  // calculation of elliptic integrals and fields (up to prefactor)
  double EllE, EllK;
  GetEllipticIntegrals(x, EllK, EllE);
  eFieldZ = EllE * 4. * dz / (a2 * b);
  eFieldR = c2 * EllE + a2 * EllK;
  // if ri = 0?
  if (ri < Small) {
    eFieldR = 0;
  } else {
    eFieldR *= 2 / (ri * a2 * b);
  }
}

void AvalancheGridSpaceCharge::GetGlobalCoordinates(double r, double z,
                                                    double phi, double &xg,
                                                    double &yg, double &zg,
                                                    int gasGap) {
  // wrt to where the center of electron number has been
  // negative r is allowed and for phi = 0 is just like the x-axis.
  // phi is wrt to local x
  // per definition local x is in global x direction.
  double xloc = r * std::cos(phi);
  // per definition local y is in global -z direction.
  double yloc = r * std::sin(phi);
  yg = z;
  xg = m_vCoNGasLayer[gasGap][0] + xloc;
  zg = m_vCoNGasLayer[gasGap][2] - yloc;
}

void AvalancheGridSpaceCharge::GetEllipticIntegrals(double x, double &K,
                                                    double &E) {
  // from x = 0 to 10 it is in steps of 1e-3. From 10 to 1e4 in steps of 1. Then
  // in steps of 1000 until 1e7.
  int arg;
  double stepSize;
  // HS: use 1. / stepSize
  if (-x < 1.e1) {
    stepSize = 1.e-3;
    arg = (int)(-x / stepSize);
  } else if (-x < 1.e4) {
    stepSize = 1.;
    arg = (int)(-x - 10) + 10000;
  } else if (-x < 1.e7) {
    stepSize = 1.e3;
    arg = (int)((-x - 1.e4) / stepSize) + 19990;
  } else {
    // not included in list.
    if (m_bDebug)
      std::cerr << m_className
                << "::GetEllipticIntegrals value not included in list.\n";
    K = m_vKElliptic.back();
    E = m_vEElliptic.back();
    return;
  }

  // linear interpolation:
  // HS: rewrite the linear interpolation,
  //     f * k[i] + (1. - f) * k[i + 1]
  K = m_vKElliptic.at(arg) +
      (-x - m_vXElliptic.at(arg)) *
          (m_vKElliptic.at(arg + 1) - m_vKElliptic.at(arg)) / (stepSize);
  E = m_vEElliptic.at(arg) +
      (-x - m_vXElliptic.at(arg)) *
          (m_vEElliptic.at(arg + 1) - m_vEElliptic.at(arg)) / (stepSize);
}

double AvalancheGridSpaceCharge::GetMeanDistance() {
  // Returns mean distance of electrons on the whole grid (doesn't work for
  // MRPCs)
  long nofElectrons = 0;
  double z = 0., meanDistance = 0.;
  for (int iz = 0; iz <= m_zSteps; iz++) {
    for (int ir = 0; ir <= m_rSteps; ir++) {
      const auto ne = m_grid[iz][ir].nElectron;
      if (ne < 0.5) continue;
      nofElectrons += ne;
      z += m_zGrid[iz] * ne;
    }
  }
  return z / (double)nofElectrons;
}

}  // namespace Garfield
