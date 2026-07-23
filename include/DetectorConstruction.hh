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
/// \file DetectorConstruction.hh
/// \brief Definition of the B1::DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;

namespace B1
{

// Standalone diagnostic (no G4RunManager/physics needed): scans the CAD
// mesh's real interior along the beam (Z) axis and reports the clear bore
// radius at each Z slice, so kChamberRadius/BeamlineLayout.hh's constants
// can be set to match the real modeled chamber instead of an assumed value.
// Writes a per-slice CSV to csvOutPath and a banded summary to G4cout. Radii
// are found by marching outward from the axis at numAngles evenly-spaced
// angles and taking the closest hit via the CAD solid's own Inside() test --
// the same test G4Navigator uses -- so the measurement matches what Geant4
// will actually do at run time.
void MeasureCadBoreProfile(const G4String& gdmlFile, const G4String& csvOutPath,
                            G4double zStep = 5. * mm, G4int numAngles = 16,
                            G4double radialStep = 2. * mm, G4double maxRadius = 700. * mm);

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction() = default;
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

    // The CAD-imported chamber shell (gdml/my_model.gdml), "CADModel" --
    // vacuum, rendering-only geometry (see Construct()). Exposed so
    // SteppingAction can stop tracks on contact for the visual run without
    // giving this volume a real material, which would perturb the physics.
    G4LogicalVolume* GetCadVolume() const { return fCadVolume; }

    // The CAD shell's own G4VSolid (placed with identity rotation/translation
    // directly under World, so World-frame coordinates need no transform to
    // use here). Exposed so SteppingAction can test a step's actual segment
    // against the solid directly, as a fallback for thin/tapered mesh
    // features the navigator's per-step volume bookkeeping can tunnel
    // through undetected -- see the comment at its use site.
    G4VSolid* GetCadSolid() const { return fCadSolid; }

  protected:
    G4LogicalVolume* fScoringVolume = nullptr;
    G4LogicalVolume* fCadVolume = nullptr;
    G4VSolid* fCadSolid = nullptr;
};

}  // namespace B1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
