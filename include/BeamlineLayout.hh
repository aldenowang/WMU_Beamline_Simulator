//
/// \file BeamlineLayout.hh
/// \brief Shared beam-path spec constants used by both DetectorConstruction
///        (to size/place the vacuum chamber and foil) and
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

// The foil sits at the exact center of the chamber, 1 m from each wall, so
// the chamber's axial half-length is just the gun-to-foil distance and the
// foil has no offset within it.
constexpr G4double kChamberAxialHalfLength = kDistanceStartToFoil;
constexpr G4double kFoilLocalZ = 0.;

}  // namespace B1::Layout

#endif
