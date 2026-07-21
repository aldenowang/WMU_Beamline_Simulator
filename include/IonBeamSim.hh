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
/// \file IonBeamSim.hh
/// \brief Definition of the B1::IonPhysicsList and B1::IonPrimaryGeneratorAction classes

#ifndef B1IonPhysicsList_h
#define B1IonPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;
class G4VPhysicalVolume;

namespace B1
{

/// Physics list for H+ ion beam / Rutherford scattering simulation.
///
/// Combines high-precision EM physics (single Coulomb scattering, needed
/// for the large-angle Rutherford tail), decay physics, ion stopping-power
/// models, hadronic elastic scattering, and stopped-particle physics.

class IonPhysicsList : public G4VModularPhysicsList
{
  public:
    IonPhysicsList();
    ~IonPhysicsList() override = default;
};

/// Primary generator for the H+ (proton) Van de Graaff beam.
///
/// Fires protons along +z from 100 cm upstream of the foil (positioned per
/// BeamlineLayout::kFoilLocalZ within the "Chamber" volume, which
/// DetectorConstruction centers on the real GDML model), with a Gaussian
/// beam spot (1 sigma = 0.5 mm), Gaussian angular divergence (1 sigma =
/// 0.5 mrad), and a Gaussian energy spread (6.0 MeV/u, sigma = 0.05%)
/// matching the Van de Graaff beam spec.

class IonPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    IonPrimaryGeneratorAction();
    ~IonPrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

  private:
    G4ParticleGun* fParticleGun = nullptr;
    G4VPhysicalVolume* fChamber = nullptr;  // cached lookup, see GeneratePrimaries()
};

}  // namespace B1

#endif
