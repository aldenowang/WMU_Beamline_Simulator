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
/// \file SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "ScatteringRecorder.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) {
    const auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  }

  // get volume of the current step
  G4LogicalVolume* volume =
    step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  // Rutherford scattering-angle scoring, restricted to primaries
  // (secondaries were never on "the beam path" to begin with) and to
  // non-animated runs: RunProduction/RunInteractiveNoVis (ionBeamSim.cc)
  // never attach a vis manager, so GetConcreteInstance()==nullptr is exactly
  // the same "is this run animated" test EventAction already uses. The
  // angle is always measured against that particle's own initial (vertex)
  // direction, not a fixed beam axis, since the gun already has a small
  // angular spread of its own.
  if (step->GetTrack()->GetParentID() == 0 && G4VVisManager::GetConcreteInstance() == nullptr) {
    const G4Track* track = step->GetTrack();
    const G4ThreeVector& vertexDir = track->GetVertexMomentumDirection();

    // The instant a primary leaves the gold foil -- the true
    // scattering event, over the full angular range including the rare
    // large-angle single-Coulomb-scattering tail Rutherford's law predicts.
    // Uses the POST-step point: the scattering happened during this step,
    // inside the foil, so the direction has to be read after it, not before.
    if (volume == fScoringVolume && step->IsLastStepInVolume()) {
      G4double angleDeg = step->GetPostStepPoint()->GetMomentumDirection().angle(vertexDir) / deg;
      ScatteringRecorder::Record(ScatteringRecorder::Channel::FoilExit, angleDeg);
    }
  }

  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
