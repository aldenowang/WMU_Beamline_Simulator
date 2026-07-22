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
/// \file IonBeamSim.cc
/// \brief Implementation of the B1::IonPhysicsList and B1::IonPrimaryGeneratorAction classes

#include "IonBeamSim.hh"

#include "BeamlineLayout.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

namespace B1
{

namespace
{
// Van de Graaff H+ beam parameters
constexpr G4double kBeamEnergyPerNucleon = 12.0 * MeV;  
constexpr G4double kBeamEnergySpread = 0.0005 * kBeamEnergyPerNucleon;  // 0.05% 
constexpr G4double kBeamSpotSigma = 0.5 * mm;  // beam radius
constexpr G4double kBeamAngularSigma = 0.5 * mrad;  
}  // namespace


IonPhysicsList::IonPhysicsList()
{
    //physics parameters for our beamline rutheford experiment
  RegisterPhysics(new G4EmStandardPhysics());  // standard EM/ionization physics (fast default)
  RegisterPhysics(new G4DecayPhysics());  // standard particle decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());  // nuclear decay (harmless, covers edge cases)
  RegisterPhysics(new G4IonPhysics());  // ion-specific stopping-power models
  RegisterPhysics(new G4HadronElasticPhysics());  // elastic nuclear scattering
  RegisterPhysics(new G4StoppingPhysics());  // stopped-particle physics

  // Use Geant4's own default range cut (1 mm) rather than an artificially
  // fine one -- the fine cut was forcing a huge number of low-energy
  // secondaries to be tracked individually, which is what was making runs
  // so slow.
}


IonPrimaryGeneratorAction::IonPrimaryGeneratorAction()
{
  fParticleGun = new G4ParticleGun(1);

  //define the particle type in the beamline
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("proton");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.)); //moving in 1 direction initially (z)
  fParticleGun->SetParticleEnergy(kBeamEnergyPerNucleon);
}


IonPrimaryGeneratorAction::~IonPrimaryGeneratorAction()
{
  delete fParticleGun;
}


void IonPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // DetectorConstruction centers the "Chamber" volume on the real GDML
  // model's bounding box (not the world origin), and places the foil inside
  // it at local z = Layout::kFoilLocalZ. Look the chamber up once (it won't
  // move between events) to get its actual world position.
  if (!fChamber) {
    fChamber = G4PhysicalVolumeStore::GetInstance()->GetVolume("Chamber");
    if (!fChamber) {
      G4Exception("IonPrimaryGeneratorAction::GeneratePrimaries()", "MyCode0007", FatalException,
                  "Chamber volume not found; cannot place beam relative to the foil.");
    }
  }
  G4ThreeVector foilPos = fChamber->GetObjectTranslation() + G4ThreeVector(0., 0., Layout::kFoilLocalZ);

  //How far the beam starts relative to the foil
  G4double z0 = foilPos.z() - Layout::kDistanceStartToFoil;

  // Transverse beam profile: independent Gaussian in x and y about the
  // foil's transverse position, 1 sigma = kBeamSpotSigma
  G4double x0 = foilPos.x() + G4RandGauss::shoot(0., kBeamSpotSigma);
  G4double y0 = foilPos.y() + G4RandGauss::shoot(0., kBeamSpotSigma);

  // Angular divergence: small-angle Gaussian tilt
  //geant4 methods, uses probabilitty and steps to determine how the particle deflects
  G4double thetaX = G4RandGauss::shoot(0., kBeamAngularSigma);
  G4double thetaY = G4RandGauss::shoot(0., kBeamAngularSigma);
  G4ThreeVector direction(thetaX, thetaY, 1.); 

  // Energy spread: Gaussian about the nominal beam energy
  G4double energy = G4RandGauss::shoot(kBeamEnergyPerNucleon, kBeamEnergySpread);

  //update particle information after striking foil
  fParticleGun->SetParticleMomentumDirection(direction.unit());
  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  fParticleGun->GeneratePrimaryVertex(event);
}


}  // namespace B1
