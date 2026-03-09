#pragma once
#include "nbody.hpp"
#include <algorithm>
#include <memory>
#include <vector>

namespace nbody {

// Represents a 3D bounding box
struct BoundingBox {
  double min_x, min_y, min_z;
  double max_x, max_y, max_z;

  bool contains(double x, double y, double z) const {
    return x >= min_x && x <= max_x && y >= min_y && y <= max_y && z >= min_z &&
           z <= max_z;
  }

  // Size of the box (maximum dimension assuming roughly cubic shape)
  double size() const {
    return std::max({max_x - min_x, max_y - min_y, max_z - min_z});
  }
};

// Octant indices (Bottom/Top, South/North, West/East)
// BSW=0, BSE=1, BNW=2, BNE=3, TSW=4, TSE=5, TNW=6, TNE=7
enum Octant {
  // Bottom (Z min)
  BSW = 0,
  BSE = 1,
  BNW = 2,
  BNE = 3,
  // Top (Z max)
  TSW = 4,
  TSE = 5,
  TNW = 6,
  TNE = 7
};

class OctNode {
public:
  BoundingBox bbox;

  // Center of mass for this node
  double com_x = 0.0;
  double com_y = 0.0;
  double com_z = 0.0;
  double total_mass = 0.0;

  int particle_idx = -1; // -1 if this is an internal node

  std::unique_ptr<OctNode> children[8];

  explicit OctNode(const BoundingBox &box) : bbox(box) {}

  bool isLeaf() const { return children[0] == nullptr; }
};

class BarnesHutTree {
public:
  explicit BarnesHutTree(const Particles &particles, double theta = 0.5);

  void build();
  void computeForces(std::vector<double> &fx, std::vector<double> &fy,
                     std::vector<double> &fz) const;

private:
  const Particles &particles_;
  std::unique_ptr<OctNode> root_;
  double theta_; // Barnes-Hut opening angle criterion

  BoundingBox computeRootBoundingBox() const;
  void insert(OctNode *node, int p_idx, const BoundingBox &box);
  void computeCenterOfMass(OctNode *node);

  void calculateForceOnParticle(int p_idx, OctNode *node, double &fx,
                                double &fy, double &fz) const;
};

} // namespace nbody
