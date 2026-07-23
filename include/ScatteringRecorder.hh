//
/// \file ScatteringRecorder.hh
/// \brief Thread-safe collector for per-particle Rutherford scattering
///        angles, recorded in SteppingAction and drained to CSV by RunAction.

#ifndef B1ScatteringRecorder_h
#define B1ScatteringRecorder_h 1

#include "G4String.hh"
#include "G4Types.hh"

namespace B1::ScatteringRecorder
{

// Recording point along the beamline:
//  - FoilExit: angle the instant a primary leaves the silicon foil -- the true
//    scattering event, over the full angular range (including the rare
//    large-angle single-Coulomb-scattering tail Rutherford's law predicts).
enum class Channel { FoilExit };

// Records one particle's scattering angle (degrees, relative to that
// particle's own initial direction) on the given channel. Safe to call from
// any worker thread.
void Record(Channel channel, G4double angleDeg);

// Writes every angle recorded so far on the given channel as a
// single-column CSV (header "angle_deg", one row per particle) to filePath,
// creating parent directories as needed, then clears that channel's
// recorded data. Call once per run, from the master thread only, after all
// worker threads have finished.
void WriteCsvAndClear(Channel channel, const G4String& filePath);

}  // namespace B1::ScatteringRecorder

#endif
