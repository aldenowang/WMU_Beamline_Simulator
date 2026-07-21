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
/// \file EventAction.hh
/// \brief Definition of the B1::EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;

namespace B1
{

class RunAction;

/// Event action class

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    void AddEdep(G4double edep) { fEdep += edep; }

  private:
    RunAction* fRunAction = nullptr;
    G4double fEdep = 0.;
    G4bool fVisualizeThisEvent = true;

    // Trajectory type (e.g. 2 = G4SmoothTrajectory) that vis.mac configured
    // via "/vis/scene/add/trajectories", captured from the first event so we
    // can restore it for visualized events after suppressing it for others.
    // -1 means "not yet captured".
    G4int fOnStoreTrajectoryValue = -1;

    // Target wall-clock length, in seconds, of the on-screen trajectory
    // buildup for a run when a viewer is attached (split evenly across
    // however many events are actually animated), so viewers can watch
    // particles travel instead of seeing the whole run flash on screen at
    // once.
    static constexpr G4double kTargetAnimationSeconds = 5.0;

    // Only the first this-many events of a run get a trajectory stored and
    // drawn. This lets a run fire millions of events for real statistics
    // while the viewer only ever has to build/animate a couple hundred
    // trajectories.
    static constexpr G4int kMaxVisualizedEvents = 200;
};

}  // namespace B1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
