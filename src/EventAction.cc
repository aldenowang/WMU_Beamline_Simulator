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
/// \file EventAction.cc
/// \brief Implementation of the B1::EventAction class

#include "EventAction.hh"

#include "RunAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4TrackingManager.hh"
#include "G4VVisManager.hh"

#include <algorithm>
#include <chrono>
#include <thread>

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction) : fRunAction(runAction) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{
  fEdep = 0.;
  fVisualizeThisEvent = event->GetEventID() < kMaxVisualizedEvents;

  // Suppress trajectory storage past the first kMaxVisualizedEvents events
  // so a multi-million-event run doesn't pay to build and hold trajectories
  // for events the viewer will never show.
  if (G4VVisManager::GetConcreteInstance() != nullptr) {
    auto* trackingManager = G4EventManager::GetEventManager()->GetTrackingManager();
    if (fOnStoreTrajectoryValue < 0) {
      // Capture whatever type vis.mac configured (e.g. 2 = smooth) the
      // first time we see it, before we ever touch it ourselves.
      fOnStoreTrajectoryValue = trackingManager->GetStoreTrajectory();
    }
    trackingManager->SetStoreTrajectory(fVisualizeThisEvent ? fOnStoreTrajectoryValue : 0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);

  // If a viewer is attached and this event is one of the ones being
  // animated, pace the trajectory buildup so the animated subset takes
  // roughly kTargetAnimationSeconds, instead of flashing on screen almost
  // instantly. This assumes events are processed serially (init_vis.mac
  // forces /run/numberOfThreads 1): with multiple worker threads sleeping
  // in parallel, the wall-clock total shrinks to roughly
  // kTargetAnimationSeconds / numberOfThreads. Events beyond
  // kMaxVisualizedEvents run at full speed with no delay.
  if (fVisualizeThisEvent && G4VVisManager::GetConcreteInstance() != nullptr) {
    const G4int nEvents =
      G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
    const G4int nAnimated = std::min(nEvents, kMaxVisualizedEvents);
    if (nAnimated > 0) {
      const auto delay = std::chrono::duration<G4double>(kTargetAnimationSeconds / nAnimated);
      std::this_thread::sleep_for(delay);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
