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

// Two independent recording points along the beamline:
//  - FoilExit: angle the instant a primary leaves the gold foil -- the true
//    scattering event, over the full angular range (including the rare
//    large-angle single-Coulomb-scattering tail Rutherford's law predicts).
//  - FaradayCup: angle the instant a primary reaches the Faraday cup -- the
//    last thing in the beamline that actually stops it. Since the cup is a
//    narrow tube on the beam axis, this channel only ever sees the small
//    forward cone (dominated by multiple scattering, not the Rutherford
//    tail), but it is still the literal "hits the thing that stops it"
//    measurement.
enum class Channel { FoilExit, FaradayCup };

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
