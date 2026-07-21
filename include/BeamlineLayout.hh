//
/// \file BeamlineLayout.hh
/// \brief Shared beam-path spec constants used by both DetectorConstruction
///        (to size/place the vacuum chamber, foil and Faraday cup) and
///        IonPrimaryGeneratorAction (to place the proton gun). Keeping a
///        single source of truth means the analytic physics geometry and
///        the simulated beam can't drift out of sync with each other.

#ifndef B1BeamlineLayout_h
#define B1BeamlineLayout_h 1

#include "G4SystemOfUnits.hh"
#include "G4Types.hh"

namespace B1::Layout
{

// Van de Graaff beamline hardware spec, measured from the gold foil target:
constexpr G4double kDistanceStartToFoil = 100. * cm;  // gun to foil (upstream)
constexpr G4double kFoilToCupDistance = 120. * cm;  // foil to Faraday cup front face
constexpr G4double kCupHalfDepth = 15. * cm / 2.;  // Faraday cup depth (15.0 cm)

// The full beam path (gun -> foil -> back of cup) that the vacuum chamber
// must contain, split into the upstream/downstream legs measured from the
// foil:
constexpr G4double kBeamPathUpstream = kDistanceStartToFoil;
constexpr G4double kBeamPathDownstream = kFoilToCupDistance + 2. * kCupHalfDepth;

// The chamber (an analytic G4Tubs) is sized to exactly span that path and
// centered on the real CAD/GDML model's bounding box (see
// DetectorConstruction::Construct()). Because the upstream and downstream
// legs are not equal, the foil (and cup) sit off-center within the chamber
// so both ends land exactly on the chamber's local +/- half-length:
constexpr G4double kChamberAxialHalfLength = (kBeamPathUpstream + kBeamPathDownstream) / 2.;
constexpr G4double kFoilLocalZ = (kBeamPathUpstream - kBeamPathDownstream) / 2.;
constexpr G4double kCupLocalZ = kFoilLocalZ + kFoilToCupDistance + kCupHalfDepth;

}  // namespace B1::Layout

#endif
