//------------------------------------------------------------------//
// Usage:
//   ./thgem [N] [chunk_size] [num_threads] [total_events]
//   ./thgem 1 1 32 10000
//------------------------------------------------------------------//

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <thread>
#include <vector>

#ifdef __linux__
  #include <pthread.h>
  #include <sched.h>
#endif

#include <TFile.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"

using namespace Garfield;

// Set to 0 if the target machine does not allow CPU affinity.
#define USE_THREAD_PINNING 1

// Runtime defaults. The command line can override all three values.
static int gChunkSize = 1;
static int gNumThreads = 32;
static int gTotalEvents = 10000;

// Geometry in cm.
static const double pitch = 0.1;
static const double metal = 0.0012;
static const double thick = 0.08;
static const double induce = 0.2;
static const double driftz = 0.5;

static const double Periodicity = 20.;
static const double xmin = -pitch / 4. * Periodicity;
static const double xmax = pitch / 4. * Periodicity;
static const double ymin = -pitch / 4. * Periodicity * std::sqrt(3.);
static const double ymax = pitch / 4. * Periodicity * std::sqrt(3.);
static const double zmin = 0.;
static const double zmax = thick + induce + driftz + metal * 2.;

struct EventMeta {
  uint32_t eOff = 0;
  uint32_t iOff = 0;
  uint32_t nep = 0;
  uint32_t nip = 0;
  int evt = -1;
  int ne = 0;
  int ni = 0;
};

struct ElectronEndpoint {
  double x1 = 0.;
  double y1 = 0.;
  double z1 = 0.;
  double t1 = 0.;
  double e1 = 0.;
  double x2 = 0.;
  double y2 = 0.;
  double z2 = 0.;
  double t2 = 0.;
  double e2 = 0.;
};

struct IonEndpoint {
  double x1 = 0.;
  double y1 = 0.;
  double z1 = 0.;
  double t1 = 0.;
  double x2 = 0.;
  double y2 = 0.;
  double z2 = 0.;
  double t2 = 0.;
};

struct ThreadResult {
  std::vector<EventMeta> events;
  std::vector<ElectronEndpoint> epts;
  std::vector<IonEndpoint> ipts;
  int eventsProcessed = 0;
};

static uint64_t SplitMix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

static uint64_t MakeRunSeed() {
  std::random_device rd;
  uint64_t seed = (uint64_t(rd()) << 32) ^ uint64_t(rd());
  seed ^= uint64_t(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
  seed ^= 0x544847454d5f5448ULL;
  return SplitMix64(seed);
}

static uint64_t ThreadSeed(const uint64_t runSeed, const int threadId,
                           const int cpuToPin) {
  uint64_t seed = runSeed;
  seed ^= uint64_t(threadId + 1) * 0x9e3779b97f4a7c15ULL;
  seed ^= uint64_t(cpuToPin + 1) * 0xbf58476d1ce4e5b9ULL;
  return SplitMix64(seed);
}

#ifdef __linux__
static void SetAffinitySingleCpu(const int cpu) {
#if USE_THREAD_PINNING
  if (cpu < 0) return;
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(cpu, &cpuset);
  (void)pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#else
  (void)cpu;
#endif
}

template <class F>
static void RunPinned(const int cpu, F&& fn) {
#if USE_THREAD_PINNING
  cpu_set_t oldSet;
  CPU_ZERO(&oldSet);
  (void)pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), &oldSet);
  SetAffinitySingleCpu(cpu);
  fn();
  (void)pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &oldSet);
#else
  (void)cpu;
  fn();
#endif
}
#else
static void SetAffinitySingleCpu(const int) {}
template <class F>
static void RunPinned(const int, F&& fn) {
  fn();
}
#endif

static std::vector<int> BuildCpuOrder() {
  unsigned int n = std::thread::hardware_concurrency();
  if (n == 0) n = 1;

  std::vector<int> cpus;
  cpus.reserve(n);
  for (unsigned int c = 0; c < n; c += 2) cpus.push_back(static_cast<int>(c));
  for (unsigned int c = 1; c < n; c += 2) cpus.push_back(static_cast<int>(c));
  return cpus;
}

static void BuildFieldAndGasPinned(const TString& path, const int pinCpu,
                                   std::unique_ptr<ComponentComsol>& fm,
                                   std::unique_ptr<MediumMagboltz>& gas) {
  RunPinned(pinCpu, [&]() {
    fm.reset(new ComponentComsol());
    const TString mfile = path + "/THGEM.mphtxt";
    const TString dfile = path + "/dielectrics.dat";
    const TString ffile = path + "/THGEM.txt";

    fm->Initialise(mfile.Data(), dfile.Data(), ffile.Data());
    fm->EnableMirrorPeriodicityX();
    fm->EnableMirrorPeriodicityY();
    fm->PrintRange();
    fm->EnableConvergenceWarnings(false);

    gas.reset(new MediumMagboltz());
    gas->SetComposition("ne", 95., "ch4", 5.);
    gas->SetTemperature(293.15);
    gas->SetPressure(760.);
    gas->SetMaxElectronEnergy(200.);
    gas->SetMaxPhotonEnergy(200.);
    gas->Initialise(true);

    const double rPenning = 0.40;
    const double lambdaPenning = 0.;
    gas->EnablePenningTransfer(rPenning, lambdaPenning, "ne");

    const char* garfieldInstall = std::getenv("GARFIELD_INSTALL");
    if (garfieldInstall) {
      const std::string mobilityFile =
          std::string(garfieldInstall) +
          "/share/Garfield/Data/IonMobility_Ne+_Ne.txt";
      gas->LoadIonMobility(mobilityFile);
    }

    fm->SetGas(gas.get());
    fm->PrintMaterials();
  });
}

static TString RootFileName(const int n) {
  return TString::Format("./result/result%d.root", n);
}

static void PrintProgressLine(const int done, const int total) {
  const int safeTotal = std::max(1, total);
  const int clampedDone = std::max(0, std::min(done, safeTotal));
  const int width = 32;
  const double frac = double(clampedDone) / double(safeTotal);
  const int filled = std::max(0, std::min(width, int(frac * width + 0.5)));

  const std::ios::fmtflags oldFlags = std::cout.flags();
  const std::streamsize oldPrecision = std::cout.precision();

  std::string bar(width, '-');
  for (int i = 0; i < filled; ++i) bar[i] = '=';

  std::cout << "\rProgress     : [" << bar << "] " << clampedDone << "/"
            << total << " (" << std::fixed << std::setprecision(1)
            << frac * 100.0 << "%)" << std::flush;

  std::cout.flags(oldFlags);
  std::cout.precision(oldPrecision);
}

struct WorkerConfig {
  int threadId = 0;
  int cpuToPin = 0;
  uint64_t rngSeed = 0;
  ComponentComsol* field = nullptr;
  std::atomic<int>* nextEvent = nullptr;
  std::atomic<int>* completedEvents = nullptr;
  ThreadResult* out = nullptr;
};

static void WorkerFunc(WorkerConfig cfg) {
  SetAffinitySingleCpu(cfg.cpuToPin);
  Garfield::Random::SetThreadSeed(cfg.rngSeed);

  std::unique_ptr<Sensor> sensor(new Sensor());
  sensor->AddComponent(cfg.field);
  sensor->SetArea(xmin, ymin, zmin, xmax, ymax, zmax);

  std::unique_ptr<AvalancheMicroscopic> aval(new AvalancheMicroscopic());
  aval->SetSensor(sensor.get());
  aval->SetShowProgress(false);

  std::unique_ptr<AvalancheMC> driftSim(new AvalancheMC());
  driftSim->SetSensor(sensor.get());
  driftSim->SetDistanceSteps(2.e-4);

  const int approxPerThread =
      (gTotalEvents + std::max(1, gNumThreads) - 1) / std::max(1, gNumThreads);
  cfg.out->events.reserve(approxPerThread);
  cfg.out->epts.reserve(static_cast<size_t>(approxPerThread) * 200);
  cfg.out->ipts.reserve(static_cast<size_t>(approxPerThread) * 200);

  while (true) {
    const int startIdx =
        cfg.nextEvent->fetch_add(gChunkSize, std::memory_order_relaxed);
    if (startIdx >= gTotalEvents) break;
    const int endIdx = std::min(startIdx + gChunkSize, gTotalEvents);

    for (int evt = startIdx; evt < endIdx; ++evt) {
      const double x0 = -pitch / 4. + pitch / 2. * Garfield::RndmUniform();
      const double y0 =
          (-pitch / 4. + pitch / 2. * Garfield::RndmUniform()) * std::sqrt(3.);
      const double z0 = zmax - 0.2 * driftz;

      aval->AvalancheElectron(x0, y0, z0, 0., 25.85, 0., 0., 0.);

      int ne = 0;
      int ni = 0;
      aval->GetAvalancheSize(ne, ni);

      const int nElectronEndpoints =
          static_cast<int>(aval->GetNumberOfElectronEndpoints());
      const int nIonEndpoints = nElectronEndpoints;

      EventMeta meta;
      meta.evt = evt;
      meta.ne = ne;
      meta.ni = ni;
      meta.nep = static_cast<uint32_t>(nElectronEndpoints);
      meta.nip = static_cast<uint32_t>(nIonEndpoints);
      meta.eOff = static_cast<uint32_t>(cfg.out->epts.size());
      meta.iOff = static_cast<uint32_t>(cfg.out->ipts.size());

      cfg.out->epts.resize(cfg.out->epts.size() +
                           static_cast<size_t>(nElectronEndpoints));
      cfg.out->ipts.resize(cfg.out->ipts.size() +
                           static_cast<size_t>(nIonEndpoints));

      for (int j = 0; j < nElectronEndpoints; ++j) {
        double xe1 = 0.;
        double ye1 = 0.;
        double ze1 = 0.;
        double te1 = 0.;
        double en1 = 0.;
        double xe2 = 0.;
        double ye2 = 0.;
        double ze2 = 0.;
        double te2 = 0.;
        double en2 = 0.;
        int status = 0;

        aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, en1, xe2, ye2, ze2,
                                  te2, en2, status);

        ElectronEndpoint ep;
        ep.x1 = xe1;
        ep.y1 = ye1;
        ep.z1 = ze1;
        ep.t1 = te1;
        ep.e1 = en1;
        ep.x2 = xe2;
        ep.y2 = ye2;
        ep.z2 = ze2;
        ep.t2 = te2;
        ep.e2 = en2;
        cfg.out->epts[meta.eOff + static_cast<uint32_t>(j)] = ep;
      }

      for (int j = 0; j < nIonEndpoints; ++j) {
        const ElectronEndpoint& src =
            cfg.out->epts[meta.eOff + static_cast<uint32_t>(j)];

        IonEndpoint ip;
        ip.x1 = src.x1;
        ip.y1 = src.y1;
        ip.z1 = src.z1;
        ip.t1 = src.t1;

        driftSim->DriftIon(src.x1, src.y1, src.z1, src.t1);

        double xi1 = 0.;
        double yi1 = 0.;
        double zi1 = 0.;
        double ti1 = 0.;
        double xi2 = 0.;
        double yi2 = 0.;
        double zi2 = 0.;
        double ti2 = 0.;
        int status = 0;
        driftSim->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2,
                                 status);

        ip.x2 = xi2;
        ip.y2 = yi2;
        ip.z2 = zi2;
        ip.t2 = ti2;
        cfg.out->ipts[meta.iOff + static_cast<uint32_t>(j)] = ip;
      }

      cfg.out->events.push_back(meta);
      ++cfg.out->eventsProcessed;
      cfg.completedEvents->fetch_add(1, std::memory_order_relaxed);
    }
  }
}

static int ParsePositiveInt(const char* value, const int fallback) {
  try {
    return std::max(1, std::stoi(value));
  } catch (...) {
    return fallback;
  }
}

int main(int argc, char* argv[]) {
  int runIndex = 0;
  if (argc >= 2) runIndex = ParsePositiveInt(argv[1], 1);
  if (argc >= 3) gChunkSize = ParsePositiveInt(argv[2], gChunkSize);
  if (argc >= 4) gNumThreads = ParsePositiveInt(argv[3], gNumThreads);
  if (argc >= 5) gTotalEvents = ParsePositiveInt(argv[4], gTotalEvents);

  TString rootFile = "./result/result.root";
  if (runIndex > 0) rootFile = RootFileName(runIndex);

  const std::vector<int> cpuOrder = BuildCpuOrder();
  int numThreads = std::min(gNumThreads, gTotalEvents);
  numThreads = std::max(1, numThreads);

  std::vector<int> threadCpu(numThreads, 0);
  for (int i = 0; i < numThreads; ++i) {
    threadCpu[i] = cpuOrder[i % cpuOrder.size()];
  }

  const TString fieldPath = "./simdata";
  const auto setupStart = std::chrono::high_resolution_clock::now();
  std::unique_ptr<ComponentComsol> field;
  std::unique_ptr<MediumMagboltz> gas;
  BuildFieldAndGasPinned(fieldPath, threadCpu[0], field, gas);
  const auto setupEnd = std::chrono::high_resolution_clock::now();
  const double setupSeconds =
      std::chrono::duration<double>(setupEnd - setupStart).count();

  const uint64_t runSeed = MakeRunSeed();
  std::vector<ThreadResult> results(numThreads);
  std::atomic<int> nextEvent(0);
  std::atomic<int> completedEvents(0);
  std::atomic<bool> progressDone(false);

  std::vector<std::thread> workers;
  workers.reserve(numThreads);

  const auto simStart = std::chrono::high_resolution_clock::now();
  PrintProgressLine(0, gTotalEvents);
  std::thread progressThread([&]() {
    while (!progressDone.load(std::memory_order_relaxed)) {
      PrintProgressLine(completedEvents.load(std::memory_order_relaxed),
                        gTotalEvents);
      std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
    PrintProgressLine(completedEvents.load(std::memory_order_relaxed),
                      gTotalEvents);
    std::cout << "\n";
  });

  for (int i = 0; i < numThreads; ++i) {
    WorkerConfig cfg;
    cfg.threadId = i;
    cfg.cpuToPin = threadCpu[i];
    cfg.rngSeed = ThreadSeed(runSeed, i, cfg.cpuToPin);
    cfg.field = field.get();
    cfg.nextEvent = &nextEvent;
    cfg.completedEvents = &completedEvents;
    cfg.out = &results[i];
    workers.push_back(std::thread(WorkerFunc, cfg));
  }

  for (size_t i = 0; i < workers.size(); ++i) {
    workers[i].join();
  }
  progressDone.store(true, std::memory_order_relaxed);
  progressThread.join();

  const auto simEnd = std::chrono::high_resolution_clock::now();
  const double simSeconds =
      std::chrono::duration<double>(simEnd - simStart).count();

  std::cout << "==== THGEM charging-up simulation ====\n";
  std::cout << "Output file  : " << rootFile.Data() << "\n";
  std::cout << "Total events : " << gTotalEvents << "\n";
  std::cout << "Threads used : " << numThreads << "\n";
  std::cout << "Chunk size   : " << gChunkSize << "\n";
  std::cout << "Setup time   : " << setupSeconds << " s\n";
  std::cout << "Sim time     : " << simSeconds << " s\n";
  std::cout << "Avg rate     : "
            << (simSeconds > 0. ? gTotalEvents / simSeconds : 0.) << " evt/s\n";

  gSystem->mkdir("result", true);
  TFile outFile(rootFile.Data(), "RECREATE");
  if (outFile.IsZombie()) {
    std::cerr << "Error: failed to create ROOT output file " << rootFile.Data()
              << "\n";
    return 1;
  }

  TTree tree("Tree", "THGEM endpoints for charging-up analysis");
  tree.SetAutoSave(0);

  std::vector<double> e1hit;
  std::vector<double> e2hit;
  std::vector<double> i1hit;
  std::vector<double> i2hit;
  int evtId = -1;
  int ne = 0;
  int ni = 0;

  tree.Branch("e1hit", &e1hit);
  tree.Branch("e2hit", &e2hit);
  tree.Branch("i1hit", &i1hit);
  tree.Branch("i2hit", &i2hit);
  tree.Branch("evt", &evtId, "evt/I");
  tree.Branch("ne", &ne, "ne/I");
  tree.Branch("ni", &ni, "ni/I");

  for (int ti = 0; ti < numThreads; ++ti) {
    const ThreadResult& tr = results[ti];
    for (size_t evIndex = 0; evIndex < tr.events.size(); ++evIndex) {
      const EventMeta& ev = tr.events[evIndex];
      evtId = ev.evt;
      ne = ev.ne;
      ni = ev.ni;

      const size_t nep = static_cast<size_t>(ev.nep);
      const size_t nip = static_cast<size_t>(ev.nip);

      e1hit.resize(nep * 5);
      e2hit.resize(nep * 5);
      i1hit.resize(nip * 4);
      i2hit.resize(nip * 4);

      for (size_t j = 0; j < nep; ++j) {
        const ElectronEndpoint& ep =
            tr.epts[ev.eOff + static_cast<uint32_t>(j)];
        const size_t base = j * 5;
        e1hit[base + 0] = ep.x1;
        e1hit[base + 1] = ep.y1;
        e1hit[base + 2] = ep.z1;
        e1hit[base + 3] = ep.t1;
        e1hit[base + 4] = ep.e1;
        e2hit[base + 0] = ep.x2;
        e2hit[base + 1] = ep.y2;
        e2hit[base + 2] = ep.z2;
        e2hit[base + 3] = ep.t2;
        e2hit[base + 4] = ep.e2;
      }

      for (size_t j = 0; j < nip; ++j) {
        const IonEndpoint& ip = tr.ipts[ev.iOff + static_cast<uint32_t>(j)];
        const size_t base = j * 4;
        i1hit[base + 0] = ip.x1;
        i1hit[base + 1] = ip.y1;
        i1hit[base + 2] = ip.z1;
        i1hit[base + 3] = ip.t1;
        i2hit[base + 0] = ip.x2;
        i2hit[base + 1] = ip.y2;
        i2hit[base + 2] = ip.z2;
        i2hit[base + 3] = ip.t2;
      }

      tree.Fill();
    }
  }

  outFile.cd();
  tree.Write();
  outFile.Close();

  return 0;
}
