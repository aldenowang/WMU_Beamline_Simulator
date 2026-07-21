//
/// \file ScatteringRecorder.cc
/// \brief Implementation of the B1::ScatteringRecorder collector.

#include "ScatteringRecorder.hh"

#include <array>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <vector>

namespace B1::ScatteringRecorder
{

namespace
{
constexpr std::size_t kNumChannels = 1;

std::array<std::mutex, kNumChannels> gMutex;
std::array<std::vector<G4double>, kNumChannels> gAnglesDeg;
}  // namespace

void Record(Channel channel, G4double angleDeg)
{
  const auto idx = static_cast<std::size_t>(channel);
  std::lock_guard<std::mutex> lock(gMutex[idx]);
  gAnglesDeg[idx].push_back(angleDeg);
}

void WriteCsvAndClear(Channel channel, const G4String& filePath)
{
  const auto idx = static_cast<std::size_t>(channel);
  std::vector<G4double> angles;
  {
    std::lock_guard<std::mutex> lock(gMutex[idx]);
    angles.swap(gAnglesDeg[idx]);
  }

  std::filesystem::path path(filePath.c_str());
  if (path.has_parent_path()) {
    std::filesystem::create_directories(path.parent_path());
  }

  std::ofstream csv(path);
  csv << "angle_deg\n";
  for (G4double angleDeg : angles) {
    csv << angleDeg << "\n";
  }
}

}  // namespace B1::ScatteringRecorder
