//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file ionBeamSim.cc
/// \brief Main program of the H+ ion beam / Rutherford scattering simulation
///
/// Usage:
///   ionBeamSim                 - interactive/animated run in a viewer window,
///                                 kDefaultAnimationEvents particles.
///   ionBeamSim -vis [N]         - same as above, N particles.
///   ionBeamSim -headless [N]    - headless production run, N particles
///                                 (default kDefaultProductionEvents),
///                                 multithreaded, no visualization: the fast
///                                 path for real statistics.
///   ionBeamSim <N>              - headless production run, N particles.
///   ionBeamSim -i [N]           - interactive Geant4 console (the "Idle>"
///                                 prompt window), N particles (default
///                                 kDefaultProductionEvents), multithreaded,
///                                 no visualization/trajectories: same speed
///                                 as the headless run but leaves you at a
///                                 live command prompt afterward.
///   ionBeamSim <file.mac>       - legacy: execute an arbitrary macro file.
///   ionBeamSim -measure-bore [out.csv]
///                               - no simulation: scans gdml/my_model.gdml's
///                                 real interior along the beam axis and
///                                 writes the clear bore radius per Z slice
///                                 to out.csv (default bore_profile.csv), so
///                                 kChamberRadius/BeamlineLayout.hh can be
///                                 set from the real model instead of a
///                                 guess. Prints a banded summary too.

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "IonBeamSim.hh"

#include "G4ExceptionHandler.hh"
#include "G4RunManagerFactory.hh"
#include "G4StateManager.hh"
#include "G4SteppingVerbose.hh"
#include "G4Threading.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include <algorithm>
#include <cstdlib>
#include <string>

using namespace B1;

namespace {

// The CAD geometry's per-track navigator/transportation warnings (e.g. a
// track getting "stuck" at a boundary the navigator is confused about) are
// G4Exceptions of severity JustWarning, printed once per occurrence with no
// built-in rate limit -- at production event counts (hundreds of thousands
// to a million) that printing is itself a meaningful chunk of the runtime.
// This handler swallows only JustWarning; everything more serious (fatal
// errors, aborted events/runs) still goes through Geant4's normal handling
// via the base class. Installed only for the tracking phase (see below), so
// the one-time geometry-construction diagnostics -- including the CAD
// overlap-check report -- still print normally.
class QuietWarningExceptionHandler : public G4ExceptionHandler
{
  public:
  G4bool Notify(const char* originOfException, const char* exceptionCode,
                G4ExceptionSeverity severity, const char* description) override
  {
    if (severity == JustWarning) return false;
    return G4ExceptionHandler::Notify(originOfException, exceptionCode, severity, description);
  }
};

}  // namespace

namespace {

// Statistics target for a headless production run: enough events for good
// large-angle Rutherford-scattering statistics in the foil.
constexpr G4int kDefaultProductionEvents = 1'000'000;

// Default event count for an animated/visualized run. EventAction only ever
// stores/draws trajectories for the first kMaxVisualizedEvents events of a
// run regardless of how many are requested (see EventAction.hh), so there is
// no visual benefit to asking for more than that.
constexpr G4int kDefaultAnimationEvents = 200;

// Parses a positive integer from a command-line argument. Returns false
// (leaving outCount untouched) if the argument isn't one, so callers can
// fall back to treating it as a macro filename instead.
G4bool ParseEventCount(const char* arg, G4int& outCount)
{
  char* end = nullptr;
  long value = std::strtol(arg, &end, 10);
  if (end == arg || *end != '\0' || value <= 0) return false;
  outCount = static_cast<G4int>(value);
  return true;
}

// Common setup shared by every run mode: detector, physics, and user
// actions. Visualization and threading are configured separately by each
// mode since they differ (headless/MT vs interactive/single-thread).
void RegisterUserInitializations(G4RunManager* runManager, G4int physicsVerbose)
{
  runManager->SetUserInitialization(new DetectorConstruction());

  auto physicsList = new IonPhysicsList;
  physicsList->SetVerboseLevel(physicsVerbose);
  runManager->SetUserInitialization(physicsList);

  runManager->SetUserInitialization(new ActionInitialization());
}

// Headless, multithreaded, no visualization: this is the fast path meant for
// real statistics. G4MTRunManager parallelizes events across all available
// cores; skipping visualization entirely means no trajectory bookkeeping or
// viewer redraw ever happens, and verbosity is kept low so terminal I/O
// doesn't become the bottleneck at high event counts.
void RunProduction(G4int nEvents)
{
  G4SteppingVerbose::UseBestUnit(4);

  auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
  runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());

  RegisterUserInitializations(runManager, /*physicsVerbose=*/0);

  auto UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/run/verbose 1");
  UImanager->ApplyCommand("/event/verbose 0");
  UImanager->ApplyCommand("/tracking/verbose 0");
  const G4int progressInterval = std::max(1, nEvents / 20);
  UImanager->ApplyCommand("/run/printProgress " + std::to_string(progressInterval));

  G4cout << G4endl << "Running " << nEvents << " particles headless on "
         << G4Threading::G4GetNumberOfCores() << " threads..." << G4endl;

  runManager->Initialize();
  // Installed after Initialize() (geometry/physics construction, including
  // the CAD overlap-check report) so only the tracking phase's per-track
  // warning spam is silenced -- see QuietWarningExceptionHandler above.
  G4StateManager::GetStateManager()->SetExceptionHandler(new QuietWarningExceptionHandler());
  runManager->BeamOn(nEvents);

  delete runManager;
}

// Interactive Geant4 console, no visualization: gives you the G4UIExecutive
// "Idle>" prompt window so you can inspect state or issue further commands
// after the run, but skips G4VisExecutive entirely. With no vis manager
// attached, EventAction's trajectory-storage/animation-pacing logic never
// triggers (it no-ops whenever G4VVisManager::GetConcreteInstance() is
// null), so this runs at full headless speed -- just like RunProduction --
// while still leaving a live console open. Multithreaded like RunProduction,
// since there is no viewer/animation to serialize around.
void RunInteractiveNoVis(G4int nEvents, int argc, char** argv)
{
  G4SteppingVerbose::UseBestUnit(4);

  auto ui = new G4UIExecutive(argc, argv);

  auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
  runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());

  RegisterUserInitializations(runManager, /*physicsVerbose=*/0);

  auto UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/run/verbose 1");
  UImanager->ApplyCommand("/event/verbose 0");
  UImanager->ApplyCommand("/tracking/verbose 0");
  const G4int progressInterval = std::max(1, nEvents / 20);
  UImanager->ApplyCommand("/run/printProgress " + std::to_string(progressInterval));

  G4cout << G4endl << "Running " << nEvents << " particles in the interactive console (no "
         << "visualization) on " << G4Threading::G4GetNumberOfCores() << " threads..." << G4endl;

  runManager->Initialize();
  // See QuietWarningExceptionHandler above: installed after Initialize() so
  // only the tracking phase's per-track warning spam is silenced.
  G4StateManager::GetStateManager()->SetExceptionHandler(new QuietWarningExceptionHandler());
  runManager->BeamOn(nEvents);

  ui->SessionStart();

  delete ui;
  delete runManager;
}

// Interactive, single-threaded, visualized: same physics/geometry as
// RunProduction(), but paced for the eye (see EventAction::EndOfEventAction)
// and capped to the first kMaxVisualizedEvents events being drawn. Forced to
// one thread because the animation pacing sleep and the shared viewer are
// not safe/meaningful to parallelize across worker threads.
void RunAnimated(G4int nEvents, int argc, char** argv)
{
  G4SteppingVerbose::UseBestUnit(4);

  auto ui = new G4UIExecutive(argc, argv);

  auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
  runManager->SetNumberOfThreads(1);

  RegisterUserInitializations(runManager, /*physicsVerbose=*/1);

  auto visManager = new G4VisExecutive(argc, argv);
  visManager->Initialize();

  auto UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/control/execute init_vis.mac");
  UImanager->ApplyCommand("/run/beamOn " + std::to_string(nEvents));

  ui->SessionStart();

  delete ui;
  delete visManager;
  delete runManager;
}

// Legacy path: execute an arbitrary macro file in batch mode, exactly as
// the original main() did for any non-flag, non-numeric argument. Kept for
// run1.mac/run2.mac-style workflows.
void RunMacroFile(const G4String& fileName, int argc, char** argv)
{
  auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  RegisterUserInitializations(runManager, /*physicsVerbose=*/1);

  auto visManager = new G4VisExecutive(argc, argv);
  visManager->Initialize();

  auto UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/control/execute " + fileName);

  delete visManager;
  delete runManager;
}

}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  if (argc >= 2) {
    const G4String flag = argv[1];
    if (flag == "-measure-bore" || flag == "--measure-bore") {
      G4String outPath = (argc >= 3) ? argv[2] : "bore_profile.csv";
      MeasureCadBoreProfile("gdml/my_model.gdml", outPath);
      return 0;
    }

    if (flag == "-vis" || flag == "--vis" || flag == "--animate") {
      G4int nEvents = kDefaultAnimationEvents;
      if (argc >= 3) ParseEventCount(argv[2], nEvents);
      RunAnimated(nEvents, argc, argv);
      return 0;
    }

    if (flag == "-i" || flag == "--interactive" || flag == "--console") {
      G4int nEvents = kDefaultProductionEvents;
      if (argc >= 3) ParseEventCount(argv[2], nEvents);
      RunInteractiveNoVis(nEvents, argc, argv);
      return 0;
    }

    if (flag == "-headless" || flag == "--headless" || flag == "--batch") {
      G4int nEvents = kDefaultProductionEvents;
      if (argc >= 3) ParseEventCount(argv[2], nEvents);
      RunProduction(nEvents);
      return 0;
    }

    G4int nEvents = 0;
    if (ParseEventCount(argv[1], nEvents)) {
      RunProduction(nEvents);
      return 0;
    }

    RunMacroFile(flag, argc, argv);
    return 0;
  }

  // Default with no arguments: open a viewer window and animate a run, so
  // just running the executable shows something instead of silently
  // crunching a million events headless.
  RunAnimated(kDefaultAnimationEvents, argc, argv);
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
