#if defined(__GPUCOMPILE__) || !defined(G_MEDIUM_MAGBOLTZ_11)

#if !defined(__GPUCOMPILE__) && !defined(G_MEDIUM_MAGBOLTZ_11)
#define G_MEDIUM_MAGBOLTZ_11
#endif

#ifndef __GPUCOMPILE__
#include <mutex>
#include <array>

#include "MagboltzInterface.hh"
#include "MediumGas.hh"
#endif

#ifndef __GPUCOMPILE__

class TPad;

namespace Garfield {

/// Interface to %Magboltz (version 11).
///  - http://magboltz.web.cern.ch/magboltz/

class MediumMagboltz : public MediumGas {
 public:
  /// Default constructor.
  MediumMagboltz();
  /// Constructor.
  MediumMagboltz(const std::string& gas1, const double f1 = 1.,
                 const std::string& gas2 = "", const double f2 = 0.,
                 const std::string& gas3 = "", const double f3 = 0.,
                 const std::string& gas4 = "", const double f4 = 0.,
                 const std::string& gas5 = "", const double f5 = 0.,
                 const std::string& gas6 = "", const double f6 = 0.);
  /// Destructor
  virtual ~MediumMagboltz() {}

  /// Set the highest electron energy to be included
  /// in the table of scattering rates.
  bool SetMaxElectronEnergy(const double e);
  /// Get the highest electron energy in the table of scattering rates.
  double GetMaxElectronEnergy() const { return m_eMax; }

  /// Set the highest photon energy to be included
  /// in the table of scattering rates.
  bool SetMaxPhotonEnergy(const double e);
  /// Get the highest photon energy in the table of scattering rates.
  double GetMaxPhotonEnergy() const { return m_eFinalGamma; }

  /// Switch on/off anisotropic scattering (enabled by default)
  void EnableAnisotropicScattering(const bool on = true) {
    m_useAnisotropic = on;
    m_isChanged = true;
  }

  /// Sample the secondary electron energy according to the Opal-Beaty model.
  void SetSplittingFunctionOpalBeaty();
  /// Sample the secondary electron energy according to the Green-Sawada model.
  void SetSplittingFunctionGreenSawada();
  /// Sample the secondary electron energy from a flat distribution.
  void SetSplittingFunctionFlat();

  /// Switch on (microscopic) de-excitation handling.
  void EnableDeexcitation();
  /// Switch off (microscopic) de-excitation handling.
  void DisableDeexcitation() { m_useDeexcitation = false; }
  /// Switch on discrete photoabsorption levels.
  void EnableRadiationTrapping();
  /// Switch off discrete photoabsorption levels.
  void DisableRadiationTrapping() { m_useRadTrap = false; }
  /// Set the number of emission line widths within which to apply 
  /// absorption by discrete lines.
  void SetLineWidth(const double n) { m_nAbsWidths = n; }

  bool EnablePenningTransfer() override;
  bool EnablePenningTransfer(const double r, const double lambda) override;
  bool EnablePenningTransfer(const double r, const double lambda,
                             std::string gasname) override;
  void DisablePenningTransfer() override;
  bool DisablePenningTransfer(std::string gasname) override;

  /// Write the gas cross-section table to a file during the initialisation.
  void EnableCrossSectionOutput(const bool on = true) { m_useCsOutput = on; }

  /// Empirical scaling of electron-multiplication probability.
  /// 1.0 keeps the official behaviour, values below 1 suppress a fraction
  /// of ionisation / Penning-created secondaries.
  void SetAvalancheProbabilityScale(const double scale) {
    m_avalancheProbabilityScale =
        scale < 0. ? 0. : (scale > 1. ? 1. : scale);
  }
  double GetAvalancheProbabilityScale() const {
    return m_avalancheProbabilityScale;
  }

  /// Multiply all excitation cross-sections by a uniform scaling factor.
  void SetExcitationScaling(const double r, std::string gasname);

  /// Initialise the table of scattering rates (called internally when a
  /// collision rate is requested and the gas mixture or other parameters
  /// have changed).
  bool Initialise(const bool verbose = false);

  void PrintGas() override;

  #endif

  #ifdef __GPUCOMPILE__
  /// Get the overall null-collision rate [ns-1].
  __device__ double GetElectronNullCollisionRate(const int band);
  #else
  double GetElectronNullCollisionRate(const int band) override;
  #endif

  #ifdef __GPUCOMPILE__
  __device__ double GetElectronCollisionRate__MediumMagboltz(const double e, const int band);
  #endif

  #ifndef __GPUCOMPILE__
  /// Get the (real) collision rate [ns-1] at a given electron energy e [eV].
  double GetElectronCollisionRate(const double e, const int band) override;
  /// Get the collision rate [ns-1] for a specific level.
  double GetElectronCollisionRate(const double e, const unsigned int level,
                                  const int band);
  /// Sample the collision type.
  bool ElectronCollision(const double e, int& type, int& level, double& e1,
                         double& dx, double& dy, double& dz,
                         std::vector<Secondary>& secondaries,
                         int& band) override;
  void ComputeDeexcitation(int iLevel, int& fLevel, 
                           std::vector<Secondary>& secondaries);

  double GetPhotonCollisionRate(const double e) override;
  bool PhotonCollision(const double e, int& type, int& level, double& e1,
                       double& ctheta, 
                       std::vector<Secondary>& secondaries) override;

  /// Switch on/off diagnostic collision counters.
  /// Disable this when a shared gas object is used only for transport.
  void EnableCollisionCounters(const bool on = true) {
    m_countCollisions = on;
  }
  /// Reset the collision counters.
  void ResetCollisionCounters();
  /// Get the total number of electron collisions.
  unsigned int GetNumberOfElectronCollisions() const;
  /// Get the number of collisions broken down by cross-section type.
  unsigned int GetNumberOfElectronCollisions(unsigned int& nElastic,
                                             unsigned int& nIonising,
                                             unsigned int& nAttachment,
                                             unsigned int& nInelastic,
                                             unsigned int& nExcitation,
                                             unsigned int& nSuperelastic) const;
  /// Get the number of cross-section terms.
  unsigned int GetNumberOfLevels();
  /// Get detailed information about a given cross-section term i.
  bool GetLevel(const unsigned int i, int& ngas, int& type, std::string& descr,
                double& e);
  /// Get the Penning transfer probability and distance of a specific level. 
  bool GetPenningTransfer(const unsigned int i, double& r, double& lambda);

  /// Get the number of collisions for a specific cross-section term.
  unsigned int GetNumberOfElectronCollisions(const unsigned int level) const;

  /// Get the number of Penning transfers that occured since the last reset.
  unsigned int GetNumberOfPenningTransfers() const { return m_nPenning; }

  /// Get the total number of photon collisions.
  unsigned int GetNumberOfPhotonCollisions() const;
  /// Get number of photon collisions by collision type.
  unsigned int GetNumberOfPhotonCollisions(unsigned int& nElastic,
                                           unsigned int& nIonising,
                                           unsigned int& nInelastic) const;

  /// Take the thermal motion of the gas at the selected temperature
  /// into account in the calculations done by Magboltz.
  /// By the default, this feature is off (static gas at 0 K).
  void EnableThermalMotion(const bool on = true) { m_useGasMotion = on; }
  /// Let Magboltz determine the upper energy limit (this is the default)
  /// or use the energy limit specified using SetMaxElectronEnergy).
  void EnableAutoEnergyLimit(const bool on = true) { m_autoEnergyLimit = on; }

  /** Run Magboltz for a given electric field, magnetic field and angle.
    * \param[in] e electric field
    * \param[in] b magnetic field
    * \param[in] btheta angle between electric and magnetic field
    * \param[in] ncoll number of collisions (in multiples of 10<sup>7</sup>)
                 to be simulated
    * \param[in] verbose verbosity flag
    * \param[out] vx,vy,vz drift velocity vector
    * \param[out] wv, wr flux and bulk drift velocity
    * \param[out] dl,dt diffusion cofficients
    * \param[out] alpha Townsend cofficient
    * \param[out] eta attachment cofficient
    * \param[out] riontof TOF ionisation rate
    * \param[out] ratttof TOF attachment rate
    * \param[out] lor Lorentz angle
    * \param[out] vxerr,vyerr,vzerr errors on drift velocity
    * \param[out] wverr, wrerr errors on TOF drifts
    * \param[out] dlerr,dterr errors on diffusion coefficients
    * \param[out] alphaerr,etaerr errors on Townsend/attachment coefficients
    * \param[out] riontoferr, ratttoferr errors on TOF rates
    * \param[out] lorerr error on Lorentz angle
    * \param[out] alphatof effective Townsend coefficient \f$(\alpha - \eta)\f$
    *             calculated using time-of-flight method
    * \param[out] difftens components of the diffusion tensor (zz, xx, yy, xz, yz, xy)
    */
  void RunMagboltz(const double e, const double b, const double btheta,
                   const int ncoll, bool verbose, double& vx, double& vy,
                   double& vz, double& wv, double& wr, double& dl, double& dt, double& alpha,
                   double& eta, double& riontof, double& ratttof, double& lor,
                   double& vxerr, double& vyerr, double& vzerr, double& wverr, double& wrerr,
                   double& dlerr, double& dterr, double& alphaerr, double& etaerr,
                   double& riontoferr, double& ratttoferr, double& lorerr,
                   double& alphatof, std::array<double, 6>& difftens);

  /// Generate a new gas table (can later be saved to file) by running
  /// Magboltz for all electric fields, magnetic fields, and
  /// angles in the currently set grid.
  void GenerateGasTable(const int numCollisions = 10,
                        const bool verbose = true);

  void PlotElectronCrossSections(const unsigned int i, TPad* pad);
  void PlotElectronCollisionRates(TPad* pad);
  void PlotElectronInverseMeanFreePath(TPad* pad);

  static int GetGasNumberMagboltz(const std::string& input);
  double CreateGPUTransferObject(MediumGPU *&med_gpu) override;
 private:
  #endif
  static constexpr int nEnergyStepsLog = 1000;
  static constexpr int nEnergyStepsGamma = 5000;
  static constexpr int nCsTypes = 7;

  #ifndef __GPUCOMPILE__
  static constexpr int nCsTypesGamma = 4;

  static const int DxcTypeRad;
  static const int DxcTypeCollIon;
  static const int DxcTypeCollNonIon;

  /// Mutex.
  std::mutex m_mutex;

  /// Simulate thermal motion of the gas or not (when running Magboltz).
  bool m_useGasMotion = false;
  /// Automatic calculation of the energy limit by Magboltz or not.
  bool m_autoEnergyLimit = true;
  #endif

/// Max. electron energy in the collision rate tables.
  double m_eMax;
  /// Energy spacing in the linear part of the collision rate tables.
  double m_eStep;
  /// Inverse energy spacing.
  double m_eStepInv;
  double m_eHigh, m_eHighLog;
  double m_lnStep;

/// Flag enabling/disabling output of cross-section table to file
  bool m_useCsOutput = false;
  double m_avalancheProbabilityScale = 1.;
  /// Number of different cross-section types in the current gas mixture
  unsigned int m_nTerms = 0;
  #ifndef __GPUCOMPILE__
  /// Mass
  std::array<double, m_nMaxGases> m_mgas;
  /// Recoil energy parameter
  std::array<double, m_nMaxGases> m_rgas;
  std::array<double, m_nMaxGases> m_s2;
   /// Opal-Beaty-Peterson splitting parameter [eV]
  std::array<double, Magboltz::nMaxLevels> m_wOpalBeaty;
  /// Green-Sawada splitting parameters [eV]
  /// (&Gamma;s, &Gamma;b, Ts, Ta, Tb).
  std::array<std::array<double, 5>, m_nMaxGases> m_parGreenSawada;
  std::array<bool, m_nMaxGases> m_hasGreenSawada;
  #endif
  /// Sample secondary electron energies using Opal-Beaty parameterisation
  bool m_useOpalBeaty = true;
  /// Sample secondary electron energies using Green-Sawada parameterisation
  bool m_useGreenSawada = false;

#ifdef __GPUCOMPILE__
  int m_csType[Magboltz::nMaxLevels];
  double m_energyLoss[Magboltz::nMaxLevels];

  double** m_scatPar{nullptr};
  int* m_numscatParIdx{nullptr};
  int m_numscatPar{0};
  
  double** m_scatCut{nullptr};
  int* m_numscatCutIdx{nullptr};
  int m_numscatCut{0};

  double** m_scatParLog{nullptr};
  int* m_numscatParLogIdx{nullptr};
  int m_numscatParLog{0};

  double** m_scatCutLog{nullptr};
  int* m_numscatCutLogIdx{nullptr};
  int m_numscatCutLog{0};
  int  m_scatModel[Magboltz::nMaxLevels];

  double* m_cfTot{nullptr};
  int m_numcfTot{0};
  double* m_cfTotLog{nullptr};
  int m_numcfTotLog{0};

  double **m_cf{nullptr};
  int* m_numcfIdx{nullptr};
  int m_numcf{0};

  double** m_cfLog{nullptr};
  int* m_numcfLogIdx{nullptr};
  int m_numcfLog{0};
#else
  /// Energy loss
  std::array<double, Magboltz::nMaxLevels> m_energyLoss;
  /// Cross-section type
  std::array<int, Magboltz::nMaxLevels> m_csType;

  /// Fluorescence yield
  std::array<double, Magboltz::nMaxLevels> m_yFluorescence;
  /// Number of Auger electrons produced in a collision
  std::array<unsigned int, Magboltz::nMaxLevels> m_nAuger1;
  std::array<unsigned int, Magboltz::nMaxLevels> m_nAuger2;
  /// Energy imparted to Auger electrons
  std::array<double, Magboltz::nMaxLevels> m_eAuger1;
  std::array<double, Magboltz::nMaxLevels> m_eAuger2;
  std::array<unsigned int, Magboltz::nMaxLevels> m_nFluorescence;
  std::array<double, Magboltz::nMaxLevels> m_eFluorescence;

  // Parameters for calculation of scattering angles
  std::vector<std::vector<double> > m_scatPar;
  std::vector<std::vector<double> > m_scatCut;
  std::vector<std::vector<double> > m_scatParLog;
  std::vector<std::vector<double> > m_scatCutLog;
  std::array<int, Magboltz::nMaxLevels> m_scatModel;

  /// Level description
  std::vector<std::string> m_description;

  // Total collision frequency
  std::vector<double> m_cfTot;
  std::vector<double> m_cfTotLog;
  // Collision frequencies
  std::vector<std::vector<double> > m_cf;
  std::vector<std::vector<double> > m_cfLog;

#endif
 bool m_useAnisotropic = true;
  /// Null-collision frequency
  double m_cfNull = 0.;
#ifdef __GPUCOMPILE__
  double m_wOpalBeaty[Magboltz::nMaxLevels];
  double m_yFluorescence[Magboltz::nMaxLevels];

  unsigned int m_nAuger1[Magboltz::nMaxLevels];
  unsigned int m_nAuger2[Magboltz::nMaxLevels];
  /// Energy imparted to Auger electrons
  double m_eAuger1[Magboltz::nMaxLevels];
  double m_eAuger2[Magboltz::nMaxLevels];

  unsigned int  m_nFluorescence[Magboltz::nMaxLevels];
  double m_eFluorescence[Magboltz::nMaxLevels];
  double m_rgas[m_nMaxGases];
  double m_s2[m_nMaxGases];
#endif

#ifndef __GPUCOMPILE__
  /// Collision counters
  /// 0: elastic
  /// 1: ionisation
  /// 2: attachment
  /// 3: inelastic
  /// 4: excitation
  /// 5: super-elastic
  std::array<unsigned int, nCsTypes> m_nCollisions;
  /// Number of collisions for each cross-section term
  std::vector<unsigned int> m_nCollisionsDetailed;
  /// Flag enabling diagnostic collision counters during transport.
  bool m_countCollisions = true;

  // Penning transfer
  /// Penning transfer probability (by level)
  std::array<double, Magboltz::nMaxLevels> m_rPenning;
  /// Mean distance of Penning ionisation (by level)
  std::array<double, Magboltz::nMaxLevels> m_lambdaPenning;
  /// Number of Penning ionisations
  unsigned int m_nPenning = 0;

  // Deexcitation
  /// Flag enabling/disabling detailed simulation of de-excitation process
  bool m_useDeexcitation = false;
  /// Flag enabling/disable radiation trapping
  /// (absorption of photons discrete excitation lines)
  bool m_useRadTrap = true;

  struct Deexcitation {
    // Gas component
    int gas;
    // Associated cross-section term
    int level;
    // Level description
    std::string label;
    // Energy
    double energy;
    // Branching ratios
    std::vector<double> p;
    // Final levels
    std::vector<int> final;
    // Type of transition
    std::vector<int> type;
    // Oscillator strength
    double osc;
    // Total decay rate
    double rate;
    // Doppler broadening
    double sDoppler;
    // Pressure broadening
    double gPressure;
    // Effective width
    double width;
    // Integrated absorption collision rate
    double cf;
  };
  std::vector<Deexcitation> m_deexcitations;
  // Mapping between deexcitations and cross-section terms.
  std::array<int, Magboltz::nMaxLevels> m_iDeexcitation;
  // Number of emission widths within which to conosider discrete line
  // absorption.
  // TODO: default value?
  double m_nAbsWidths = 1000.;

  /// Ionisation potentials of each component
  std::array<double, m_nMaxGases> m_ionPot;
  /// Minimum ionisation potential
  double m_minIonPot = -1.;

  // Scaling factor for excitation cross-sections
  std::array<double, m_nMaxGases> m_scaleExc;

  // Energy spacing of photon collision rates table
  double m_eFinalGamma, m_eStepGamma;
  // Number of photon collision cross-section terms
  unsigned int m_nPhotonTerms = 0;
  // Total photon collision frequencies
  std::vector<double> m_cfTotGamma;
  // Photon collision frequencies
  std::vector<std::vector<double> > m_cfGamma;
  std::vector<int> csTypeGamma;
  // Photon collision counters
  // 0: elastic
  // 1: ionisation
  // 2: inelastic
  // 3: excitation
  std::array<unsigned int, nCsTypesGamma> m_nPhotonCollisions;

  bool Update(const bool verbose = false);
  bool Mixer(const bool verbose = false);
  void SetupGreenSawada();

  void GetExcitationIonisationLevels();

  void ComputeDeexcitationTable(const bool verbose);
  void AddPenningDeexcitation(Deexcitation& dxc, const double rate,
                              const double pPenning) {
    dxc.p.push_back(rate * (1. - pPenning));
    dxc.type.push_back(DxcTypeCollNonIon);
    dxc.final.push_back(-1);
    if (pPenning > 0.) { 
      dxc.p.push_back(rate * pPenning);
      dxc.type.push_back(DxcTypeCollIon);
      dxc.final.push_back(-1);
    }
  }
  double RateConstantWK(const double energy, const double osc,
                        const double pacs, const int igas1,
                        const int igas2) const;
  double RateConstantHardSphere(const double r1, const double r2,
                                const int igas1, const int igas2) const;
  double CalcDiscreteLineCf(const Deexcitation& dxc,
                            double e, double cfOth) const;
  void ComputeDeexcitationInternal(int iLevel, int& fLevel,
      std::vector<Secondary>& secondaries);
  bool ComputePhotonCollisionTable(const bool verbose);
};
}
#endif
#endif
