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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "BeamlineLayout.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4TessellatedSolid.hh"
#include "G4ThreeVector.hh"
#include "G4TriangularFacet.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <queue>
#include <regex>
#include <unordered_map>
#include <vector>

namespace B1
{

namespace
{
// Vacuum scattering chamber. G4Tubs is built along its local Z axis, which we
// keep aligned with the beam (matching the primary generator and foil
// placement below). Its axial half-length comes from BeamlineLayout.hh,
// sized so the foil sits at the exact chamber center, 1 m from each wall;
// the chamber is then centered on the actual CAD/GDML model's bounding box
// (computed at load time below) so the physics geometry fits inside the
// real model rather than an independent guess. "Height" becomes the
// transverse radius.
constexpr G4double kChamberRadius = 100. * cm / 2.;  // 100.0 cm height spec

// Gold foil target: single sheet, perpendicular to the beam, at the
// chamber's own center (BeamlineLayout::kFoilLocalZ == 0)
constexpr G4double kFoilHalfXY = 1.905 * cm / 2.;  // 1.905 cm x 1.905 cm
constexpr G4double kFoilHalfZ = 100. * um / 2.;  // 1.0 um thickness

// World must contain both the analytic chamber assembly above and the CAD
// overlay loaded from gdml/my_model.gdml below. As authored, that mesh's
// long axis (its chamber pipe) runs along Y, not Z, so it is rotated by
// kCadRotationAngle below to bring the pipe onto the beam axis. After that
// rotation the mesh's vertex bounding box is roughly x:[-1284,1216]
// y:[-1064,1286] z:[-350,2194] mm (it is not centered on the origin), so the
// world is sized to that, with margin, rather than just around the chamber.
constexpr G4double kWorldHalfX = 1500. * mm;
constexpr G4double kWorldHalfY = 2500. * mm;
constexpr G4double kWorldHalfZ = 2500. * mm;

// The CAD model is authored with its chamber pipe along Y; rotating 90 deg
// about X swaps Y and Z (with a sign) so the pipe lines up with the beam,
// which travels along +Z (see IonPrimaryGeneratorAction), with the pipe/gun
// end landing on the beam-start side and the rest of the assembly (wall)
// downstream. If the model comes out mirrored front-to-back once rendered
// (beam entering through the wall end instead of the pipe), flip the sign
// to -90. * deg.
constexpr G4double kCadRotationAngle = 90. * deg;

// Axis-aligned bounding box of the CAD mesh's vertices, gathered while
// parsing the GDML below.
struct GdmlBounds
{
  G4double minX = 0., maxX = 0., minY = 0., maxY = 0., minZ = 0., maxZ = 0.;
  G4bool valid = false;

  G4ThreeVector Center() const
  {
    return G4ThreeVector((minX + maxX) / 2., (minY + maxY) / 2., (minZ + maxZ) / 2.);
  }
};

// Encodes an undirected edge between welded-vertex ids (a,b) into a single
// key, independent of which order (a,b) is asked in.
std::uint64_t EdgeKey(G4int a, G4int b)
{
  if (a > b) std::swap(a, b);
  return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint32_t>(b);
}

// Repairs facet winding so the whole mesh has a single consistent notion of
// "outward". CAD software merging several parts (pipe, flanges, ports) into
// one exported mesh commonly leaves a subset of facets with inverted
// winding; G4TessellatedSolid::SetSolidClosed() already flags this ("some
// facets have wrong orientation") but does not fix it, and an inconsistently
// wound mesh is exactly what corrupts inside/outside classification for the
// navigator -- confirmed here by G4PVPlacement::CheckOverlaps() reporting
// the CAD shell overlapping the analytic Chamber by tens of centimetres, and
// by G4Navigator::ComputeStep() logging repeated "Stuck Track ... Likely
// geometry overlap" warnings for tracks crossing into CADModel.
//
// This uses the standard mesh-repair approach (not a per-facet/convexity
// heuristic, which was tried first and made no measurable difference on
// this non-convex, elongated geometry): flood-fill the facet-adjacency graph
// so every pair of facets sharing an edge is wound in opposite directions
// along that edge (the definition of consistent orientation for a manifold
// mesh), then, per connected component, pick the one global in/out flip
// that makes its net (divergence-theorem) signed volume positive.
//
// Vertex adjacency is matched by *welded* (position-rounded) id rather than
// by raw GDML vertex name/index, because parts that were merged without
// welding can have several distinctly-named vertices sitting at the same
// point -- matching by raw index would see those as unconnected and miss
// the seam entirely.
void RepairFacetOrientation(const std::vector<G4ThreeVector>& vertexPos,
                             std::vector<std::array<G4int, 3>>& faceVerts, std::size_t& flippedOut)
{
  const std::size_t numFaces = faceVerts.size();

  // Weld vertices that sit at (numerically) the same point so adjacency
  // across a badly-merged seam is still found.
  std::map<std::array<long long, 3>, G4int> weldMap;
  std::vector<G4int> weldOf(vertexPos.size());
  for (std::size_t i = 0; i < vertexPos.size(); ++i) {
    const G4ThreeVector& p = vertexPos[i];
    std::array<long long, 3> key{std::llround(p.x() * 1000.), std::llround(p.y() * 1000.),
                                  std::llround(p.z() * 1000.)};  // 1 micron buckets
    auto [it, inserted] = weldMap.try_emplace(key, static_cast<G4int>(i));
    weldOf[i] = it->second;
  }

  struct EdgeOcc
  {
    G4int face;
    G4int from, to;  // this face's original directed order along the edge
  };
  std::unordered_map<std::uint64_t, std::vector<EdgeOcc>> edgeMap;
  edgeMap.reserve(numFaces * 3);
  for (std::size_t f = 0; f < numFaces; ++f) {
    G4int a = weldOf[faceVerts[f][0]], b = weldOf[faceVerts[f][1]], c = weldOf[faceVerts[f][2]];
    edgeMap[EdgeKey(a, b)].push_back({static_cast<G4int>(f), a, b});
    edgeMap[EdgeKey(b, c)].push_back({static_cast<G4int>(f), b, c});
    edgeMap[EdgeKey(c, a)].push_back({static_cast<G4int>(f), c, a});
  }

  std::vector<G4bool> flipped(numFaces, false);
  std::vector<G4bool> visited(numFaces, false);

  for (std::size_t seed = 0; seed < numFaces; ++seed) {
    if (visited[seed]) continue;

    std::vector<std::size_t> component;
    std::queue<std::size_t> toVisit;
    visited[seed] = true;
    toVisit.push(seed);

    while (!toVisit.empty()) {
      std::size_t f = toVisit.front();
      toVisit.pop();
      component.push_back(f);

      G4int v[3] = {weldOf[faceVerts[f][0]], weldOf[faceVerts[f][1]], weldOf[faceVerts[f][2]]};
      if (flipped[f]) std::swap(v[1], v[2]);  // this face's *effective* winding

      for (int e = 0; e < 3; ++e) {
        G4int x = v[e], y = v[(e + 1) % 3];
        for (const EdgeOcc& occ : edgeMap[EdgeKey(x, y)]) {
          if (static_cast<std::size_t>(occ.face) == f || visited[occ.face]) continue;
          // Consistent orientation means the shared edge runs opposite ways
          // across its two facets. If this neighbour's *original* direction
          // matches ours instead, it must be flipped to agree with us.
          flipped[occ.face] = (occ.from == x && occ.to == y);
          visited[occ.face] = true;
          toVisit.push(occ.face);
        }
      }
    }

    // Decide the one remaining global bit for this component -- in vs. out
    // -- from the sign of its net signed volume (divergence theorem), which
    // stays reliable even for a concave/elongated shape since it's an
    // aggregate over every facet in the component rather than a local,
    // per-facet or per-centroid guess.
    G4double signedVolume6 = 0.;
    for (std::size_t f : component) {
      G4int a = faceVerts[f][0], b = faceVerts[f][1], c = faceVerts[f][2];
      if (flipped[f]) std::swap(b, c);
      signedVolume6 += vertexPos[a].dot(vertexPos[b].cross(vertexPos[c]));
    }
    if (signedVolume6 < 0.) {
      for (std::size_t f : component) flipped[f] = !flipped[f];
    }
  }

  for (std::size_t f = 0; f < numFaces; ++f) {
    if (flipped[f]) {
      std::swap(faceVerts[f][1], faceVerts[f][2]);
      ++flippedOut;
    }
  }
}

// Parses the subset of GDML used by gdml/my_model.gdml (a single tessellated
// mesh exported from CAD: <position> vertices + <triangular> facets, all in
// mm) directly, without a full GDML/XML parser. This Geant4 build has no
// GDML support (no Xerces-C), so G4GDMLParser is unavailable here.
//
// Also fills `bounds` with the axis-aligned bounding box of the mesh's
// vertices, so the caller can size and place the analytic physics geometry
// to fit inside the real (imported) model instead of an independent guess.
G4VSolid* BuildTessellatedSolidFromGdml(const G4String& filename, GdmlBounds& bounds)
{
  std::ifstream in(filename);
  if (!in.is_open()) {
    G4Exception("DetectorConstruction::BuildTessellatedSolidFromGdml", "MyCode0003",
                FatalException, ("Could not open GDML file: " + filename).c_str());
  }

  static const std::regex posRe(
    R"re(<position\s+name="([^"]+)"\s+x="([^"]+)"\s+y="([^"]+)"\s+z="([^"]+)")re");
  static const std::regex triRe(
    R"re(<triangular\s+vertex1="([^"]+)"\s+vertex2="([^"]+)"\s+vertex3="([^"]+)")re");

  std::unordered_map<std::string, G4int> vertexIndex;
  vertexIndex.reserve(30000);
  std::vector<G4ThreeVector> vertexPos;
  vertexPos.reserve(30000);
  std::vector<std::array<G4int, 3>> faceVerts;
  faceVerts.reserve(60000);

  std::smatch match;
  std::string line;

  const G4double cosCadRot = std::cos(kCadRotationAngle);
  const G4double sinCadRot = std::sin(kCadRotationAngle);

  while (std::getline(in, line)) {
    if (std::regex_search(line, match, posRe)) {
      G4double x = std::stod(match[2].str()) * mm;
      G4double yRaw = std::stod(match[3].str()) * mm;
      G4double zRaw = std::stod(match[4].str()) * mm;
      // Rotate about X by kCadRotationAngle: y' = y*cos - z*sin, z' = y*sin + z*cos
      G4double y = cosCadRot * yRaw - sinCadRot * zRaw;
      G4double z = sinCadRot * yRaw + cosCadRot * zRaw;
      if (!bounds.valid) {
        bounds.minX = bounds.maxX = x;
        bounds.minY = bounds.maxY = y;
        bounds.minZ = bounds.maxZ = z;
        bounds.valid = true;
      }
      else {
        bounds.minX = std::min(bounds.minX, x);
        bounds.maxX = std::max(bounds.maxX, x);
        bounds.minY = std::min(bounds.minY, y);
        bounds.maxY = std::max(bounds.maxY, y);
        bounds.minZ = std::min(bounds.minZ, z);
        bounds.maxZ = std::max(bounds.maxZ, z);
      }
      vertexIndex.emplace(match[1].str(), static_cast<G4int>(vertexPos.size()));
      vertexPos.emplace_back(x, y, z);
      continue;
    }
    if (std::regex_search(line, match, triRe)) {
      auto v1 = vertexIndex.find(match[1].str());
      auto v2 = vertexIndex.find(match[2].str());
      auto v3 = vertexIndex.find(match[3].str());
      if (v1 == vertexIndex.end() || v2 == vertexIndex.end() || v3 == vertexIndex.end()) continue;
      faceVerts.push_back({v1->second, v2->second, v3->second});
    }
  }

  if (!bounds.valid) {
    G4Exception("DetectorConstruction::BuildTessellatedSolidFromGdml", "MyCode0004",
                FatalException, ("GDML file has no vertices: " + filename).c_str());
  }

  std::size_t flippedFacetCount = 0;
  RepairFacetOrientation(vertexPos, faceVerts, flippedFacetCount);

  auto solid = new G4TessellatedSolid("CADModel_solid");
  for (const auto& f : faceVerts) {
    solid->AddFacet(
      new G4TriangularFacet(vertexPos[f[0]], vertexPos[f[1]], vertexPos[f[2]], ABSOLUTE));
  }
  solid->SetSolidClosed(true);

  G4cout << "DetectorConstruction: loaded CAD model \"" << filename << "\" (" << faceVerts.size()
         << " facets, " << flippedFacetCount << " re-oriented)" << G4endl;

  return solid;
}
}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MeasureCadBoreProfile(const G4String& gdmlFile, const G4String& csvOutPath, G4double zStep,
                            G4int numAngles, G4double radialStep, G4double maxRadius)
{
  GdmlBounds bounds;
  auto solidCad = BuildTessellatedSolidFromGdml(gdmlFile, bounds);

  std::ofstream csv(csvOutPath);
  if (!csv.is_open()) {
    G4Exception("MeasureCadBoreProfile", "MyCode0008", FatalException,
                ("Could not open output CSV for writing: " + csvOutPath).c_str());
  }
  csv << "z_mm,clear_radius_mm";
  for (G4int a = 0; a < numAngles; ++a) csv << ",r_theta" << a << "_mm";
  csv << "\n";

  G4cout << "MeasureCadBoreProfile: scanning z=[" << bounds.minZ / mm << "," << bounds.maxZ / mm
         << "] mm, " << numAngles << " angles, radial step " << radialStep / mm
         << " mm, cap " << maxRadius / mm << " mm" << G4endl;

  // Direction unit vectors for each sampled angle, precomputed once.
  std::vector<G4ThreeVector> dirs(numAngles);
  for (G4int a = 0; a < numAngles; ++a) {
    G4double theta = 2. * CLHEP::pi * a / numAngles;
    dirs[a] = G4ThreeVector(std::cos(theta), std::sin(theta), 0.);
  }

  // Banded summary: report contiguous z-runs where the clear radius is below
  // vs. at/above the currently-configured kChamberRadius, so the mismatch
  // (or a real nozzle/chamber-body transition) is easy to see without
  // reading the whole CSV.
  G4bool bandOpen = false;
  G4bool bandBelowThreshold = false;
  G4double bandStartZ = 0.;
  G4double bandMinRadius = 0., bandMaxRadius = 0.;

  auto flushBand = [&](G4double endZ) {
    if (!bandOpen) return;
    G4cout << "  z=[" << bandStartZ / mm << ", " << endZ / mm << "] mm : clear radius "
           << bandMinRadius / mm << "-" << bandMaxRadius / mm << " mm "
           << (bandBelowThreshold ? "(BELOW kChamberRadius)" : "(fits kChamberRadius)") << G4endl;
  };

  for (G4double z = bounds.minZ; z <= bounds.maxZ; z += zStep) {
    G4double clearRadius = maxRadius;
    std::vector<G4double> hitRadii(numAngles, maxRadius);

    for (G4int a = 0; a < numAngles; ++a) {
      G4double hit = maxRadius;
      for (G4double r = 0.; r <= maxRadius; r += radialStep) {
        G4ThreeVector p(0., 0., z);
        p += r * dirs[a];
        if (solidCad->Inside(p) != kOutside) {
          hit = r;
          break;
        }
      }
      hitRadii[a] = hit;
      clearRadius = std::min(clearRadius, hit);
    }

    csv << z / mm << "," << clearRadius / mm;
    for (G4int a = 0; a < numAngles; ++a) csv << "," << hitRadii[a] / mm;
    csv << "\n";

    G4bool belowThreshold = clearRadius < kChamberRadius;
    if (!bandOpen) {
      bandOpen = true;
      bandBelowThreshold = belowThreshold;
      bandStartZ = z;
      bandMinRadius = bandMaxRadius = clearRadius;
    }
    else if (belowThreshold != bandBelowThreshold) {
      flushBand(z);
      bandBelowThreshold = belowThreshold;
      bandStartZ = z;
      bandMinRadius = bandMaxRadius = clearRadius;
    }
    else {
      bandMinRadius = std::min(bandMinRadius, clearRadius);
      bandMaxRadius = std::max(bandMaxRadius, clearRadius);
    }
  }
  flushBand(bounds.maxZ);

  G4cout << "MeasureCadBoreProfile: done, wrote \"" << csvOutPath << "\" (reference kChamberRadius="
         << kChamberRadius / mm << " mm)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");  // high vacuum chamber

  G4bool checkOverlaps = true;

  //
  // World
  //
  auto solidWorld = new G4Box("World", kWorldHalfX, kWorldHalfY, kWorldHalfZ);
  auto logicWorld = new G4LogicalVolume(solidWorld, vacuum, "World");

  auto physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false,
                                     0, checkOverlaps);

  //
  // CAD model (parsed first: the analytic physics geometry below is sized
  // and placed to fit inside its bounding box, so it stays consistent with
  // the real imported model instead of an independent guess)
  //
  GdmlBounds gdmlBounds;
  auto solidCad = BuildTessellatedSolidFromGdml("gdml/my_model.gdml", gdmlBounds);

  // Transverse (world X,Y) placement: NOT the mesh's bounding-box centroid.
  // The CAD assembly's overall bounding box is skewed off the true pipe bore
  // by asymmetric flanges/ports elsewhere on the chamber body (measured
  // centroid offset ~(-34, 111) mm), but the pipe itself -- the narrow entry
  // bore the beam actually needs to travel through -- is modeled dead-center
  // on the mesh's own local (x=0, z=0) axis (verified against the raw GDML
  // vertices: the pipe's cross-section slices centroid to (0,0) to within
  // rounding). Since only X,Y_raw<->Z_raw are mixed by kCadRotationAngle
  // (rotation about X), local x=0 maps straight to world x=0, and local
  // z=0 maps to world y=0. So the analytic chamber -- and everything placed
  // relative to it (foil, cup, and the beam gun via
  // IonPrimaryGeneratorAction) -- must sit at world (x=0, y=0) transversely;
  // only the axial (Z) position comes from the mesh's bounding box.
  G4ThreeVector chamberPos(0., 0., gdmlBounds.Center().z());

  // Verify the beam-path chamber (sized from BeamlineLayout.hh) actually
  // fits inside the real model before placing anything, rather than
  // silently clipping through its walls.
  G4double transverseMargin =
    std::min({gdmlBounds.maxX - chamberPos.x(), chamberPos.x() - gdmlBounds.minX,
              gdmlBounds.maxY - chamberPos.y(), chamberPos.y() - gdmlBounds.minY});
  if (kChamberRadius > transverseMargin) {
    G4Exception("DetectorConstruction::Construct", "MyCode0005", FatalException,
                "Chamber radius does not fit within the GDML model's transverse extent");
  }
  G4double axialMargin = (gdmlBounds.maxZ - gdmlBounds.minZ) / 2.;
  if (Layout::kChamberAxialHalfLength > axialMargin) {
    G4Exception("DetectorConstruction::Construct", "MyCode0006", FatalException,
                "Beam path (gun to foil) does not fit within the GDML model's axial extent");
  }

  //
  // Chamber (vacuum), centered on the CAD model's bounding box
  //
  auto solidChamber =
    new G4Tubs("Chamber", 0., kChamberRadius, Layout::kChamberAxialHalfLength, 0., 360. * deg);
  auto logicChamber = new G4LogicalVolume(solidChamber, vacuum, "Chamber");

  new G4PVPlacement(nullptr, chamberPos, logicChamber, "Chamber", logicWorld, false, 0,
                    checkOverlaps);

  //
  // Gold foil target (single sheet; swap "G4_Au" for "G4_C" for a carbon foil)
  //
  G4Material* foilMat = nist->FindOrBuildMaterial("G4_Au");

  auto solidFoil = new G4Box("Foil", kFoilHalfXY, kFoilHalfXY, kFoilHalfZ);
  auto logicFoil = new G4LogicalVolume(solidFoil, foilMat, "Foil");

  // Give the foil an actual gold colour and force it solid; at 1 um thick
  // it's otherwise so thin it can disappear entirely under default vis
  // attributes, especially once seen through the semi-transparent CAD shell.
  auto foilVis = new G4VisAttributes(G4Colour(1.0, 0.84, 0.0));
  foilVis->SetForceSolid(true);
  logicFoil->SetVisAttributes(foilVis);

  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., Layout::kFoilLocalZ), logicFoil, "Foil",
                    logicChamber, false, 0, checkOverlaps);

  //
  // CAD model overlay (rendering only; kept as "vacuum" material rather than
  // the GDML's own G4_Al so it cannot perturb the beam physics above even if
  // its mesh geometrically clips the chamber/foil/cup). Placed at the
  // world origin since its vertices are already in the model's own
  // (world) coordinate system.
  //
  auto logicCad = new G4LogicalVolume(solidCad, vacuum, "CADModel");

  // Alpha kept moderate so the shell reads as solid in the UI while the
  // gold foil target and beam inside remain visible through it.
  auto cadVis = new G4VisAttributes(G4Colour(0.6, 0.6, 0.6, 0.35));
  cadVis->SetForceSolid(true);
  logicCad->SetVisAttributes(cadVis);

  // checkOverlaps=false here: the default checkOverlaps=true path samples
  // only 1000 points, too sparse for a ~58k-facet mesh with thin features
  // like the beam-pipe bore to reliably probe. Check explicitly below with a
  // much higher sample count instead.
  auto physCad = new G4PVPlacement(nullptr, G4ThreeVector(), logicCad, "CADModel", logicWorld,
                                   false, 0, false);
  physCad->CheckOverlaps(100000, 0., true);

  // Score energy deposited in the gold foil
  fScoringVolume = logicFoil;

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
