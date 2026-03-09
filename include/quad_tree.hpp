#pragma once
#include "nbody.hpp"
#include <memory>
#include <vector>

namespace nbody {

// Represents a 2D bounding box
struct BoundingBox {
  double min_x, min_y;
  double max_x, max_y;

  bool contains(double x, double y) const {
    return x >= min_x && x <= max_x && y >= min_y && y <= max_y;
  }

  // Size of the box (width assuming square or max dimension)
  double size() const { return std::max(max_x - min_x, max_y - min_y); }
};

// Quadrant indices
enum Quadrant { NW = 0, NE = 1, SW = 2, SE = 3 };

class QuadNode {
public:
  BoundingBox bbox;

  // Center of mass for this node
  double com_x = 0.0;
  double com_y = 0.0;
  double total_mass = 0.0;

  int particle_idx = -1; // -1 if this is an internal node

  std::unique_ptr<QuadNode> children[4];

  explicit QuadNode(const BoundingBox &box) : bbox(box) {}

  bool isLeaf() const { return children[0] == nullptr; }
};

class BarnesHutTree {
public:
  explicit BarnesHutTree(const Particles &particles, double theta = 0.5);

  void build();
  void computeForces(std::vector<double> &fx, std::vector<double> &fy)
      const; // z is ignored for 2D Barnes-Hut

private:
  const Particles &particles_;
  std::unique_ptr<QuadNode> root_;
  double theta_; // Barnes-Hut opening angle criterion

  BoundingBox computeRootBoundingBox() const;
  void insert(QuadNode *node, int p_idx, const BoundingBox &box);
  void computeCenterOfMass(QuadNode *node);

  void calculateForceOnParticle(int p_idx, QuadNode *node, double &fx,
                                double &fy) const;
};

} // namespace nbody
