#include "oct_tree.hpp"
#include <cmath>
#include <iostream>

namespace nbody {

BarnesHutTree::BarnesHutTree(const Particles &particles, double theta)
    : particles_(particles), theta_(theta) {
  build();
}

BoundingBox BarnesHutTree::computeRootBoundingBox() const {
  if (particles_.size == 0)
    return {0, 0, 0, 0, 0, 0};

  double min_x = particles_.x[0], max_x = particles_.x[0];
  double min_y = particles_.y[0], max_y = particles_.y[0];
  double min_z = particles_.z[0], max_z = particles_.z[0];

  for (std::size_t i = 1; i < particles_.size; ++i) {
    if (particles_.x[i] < min_x)
      min_x = particles_.x[i];
    if (particles_.x[i] > max_x)
      max_x = particles_.x[i];

    if (particles_.y[i] < min_y)
      min_y = particles_.y[i];
    if (particles_.y[i] > max_y)
      max_y = particles_.y[i];

    if (particles_.z[i] < min_z)
      min_z = particles_.z[i];
    if (particles_.z[i] > max_z)
      max_z = particles_.z[i];
  }

  // Make it a cube to simplify octree subdivision
  double width = max_x - min_x;
  double height = max_y - min_y;
  double depth = max_z - min_z;
  double max_dim = std::max({width, height, depth});

  return {min_x,           min_y,           min_z,
          min_x + max_dim, min_y + max_dim, min_z + max_dim};
}

void BarnesHutTree::build() {
  if (particles_.size == 0)
    return;

  BoundingBox root_box = computeRootBoundingBox();
  root_ = std::make_unique<OctNode>(root_box);

  for (std::size_t i = 0; i < particles_.size; ++i) {
    insert(root_.get(), i, root_box);
  }

  computeCenterOfMass(root_.get());
}

void BarnesHutTree::insert(OctNode *node, int p_idx, const BoundingBox &box) {
  if (node->isLeaf()) {
    if (node->particle_idx == -1) {
      // Empty leaf, just put the particle here
      node->particle_idx = p_idx;
    } else {
      // Collision: Leaf already has a particle. We must subdivide.
      int existing_p = node->particle_idx;
      node->particle_idx = -1; // It becomes an internal node

      double mid_x = (box.min_x + box.max_x) / 2.0;
      double mid_y = (box.min_y + box.max_y) / 2.0;
      double mid_z = (box.min_z + box.max_z) / 2.0;

      // Bottom (Z min to mid_z)
      node->children[BSW] = std::make_unique<OctNode>(
          BoundingBox{box.min_x, box.min_y, box.min_z, mid_x, mid_y, mid_z});
      node->children[BSE] = std::make_unique<OctNode>(
          BoundingBox{mid_x, box.min_y, box.min_z, box.max_x, mid_y, mid_z});
      node->children[BNW] = std::make_unique<OctNode>(
          BoundingBox{box.min_x, mid_y, box.min_z, mid_x, box.max_y, mid_z});
      node->children[BNE] = std::make_unique<OctNode>(
          BoundingBox{mid_x, mid_y, box.min_z, box.max_x, box.max_y, mid_z});

      // Top (mid_z to Z max)
      node->children[TSW] = std::make_unique<OctNode>(
          BoundingBox{box.min_x, box.min_y, mid_z, mid_x, mid_y, box.max_z});
      node->children[TSE] = std::make_unique<OctNode>(
          BoundingBox{mid_x, box.min_y, mid_z, box.max_x, mid_y, box.max_z});
      node->children[TNW] = std::make_unique<OctNode>(
          BoundingBox{box.min_x, mid_y, mid_z, mid_x, box.max_y, box.max_z});
      node->children[TNE] = std::make_unique<OctNode>(
          BoundingBox{mid_x, mid_y, mid_z, box.max_x, box.max_y, box.max_z});

      // Re-insert the old particle
      insert(node, existing_p, box);
      // Insert the new particle
      insert(node, p_idx, box);
    }
  } else {
    // Internal node: find which octant the particle belongs to
    double mid_x = (box.min_x + box.max_x) / 2.0;
    double mid_y = (box.min_y + box.max_y) / 2.0;
    double mid_z = (box.min_z + box.max_z) / 2.0;

    bool is_east = particles_.x[p_idx] > mid_x;  // E=1, W=0
    bool is_north = particles_.y[p_idx] > mid_y; // N=2, S=0
    bool is_top = particles_.z[p_idx] > mid_z;   // T=4, B=0

    int oct_idx = (is_east ? 1 : 0) + (is_north ? 2 : 0) + (is_top ? 4 : 0);

    insert(node->children[oct_idx].get(), p_idx, node->children[oct_idx]->bbox);
  }
}

void BarnesHutTree::computeCenterOfMass(OctNode *node) {
  if (!node)
    return;

  if (node->isLeaf()) {
    if (node->particle_idx != -1) {
      node->com_x = particles_.x[node->particle_idx];
      node->com_y = particles_.y[node->particle_idx];
      node->com_z = particles_.z[node->particle_idx];
      node->total_mass = particles_.mass[node->particle_idx];
    }
    return;
  }

  node->com_x = 0;
  node->com_y = 0;
  node->com_z = 0;
  node->total_mass = 0;

  for (int i = 0; i < 8; ++i) {
    if (node->children[i]) {
      computeCenterOfMass(node->children[i].get());
      double m = node->children[i]->total_mass;
      if (m > 0) {
        node->com_x += node->children[i]->com_x * m;
        node->com_y += node->children[i]->com_y * m;
        node->com_z += node->children[i]->com_z * m;
        node->total_mass += m;
      }
    }
  }

  if (node->total_mass > 0) {
    node->com_x /= node->total_mass;
    node->com_y /= node->total_mass;
    node->com_z /= node->total_mass;
  }
}

void BarnesHutTree::calculateForceOnParticle(int p_idx, OctNode *node,
                                             double &fx, double &fy,
                                             double &fz) const {
  if (!node || node->total_mass == 0)
    return;

  // If it's a leaf containing the current particle itself, skip
  if (node->isLeaf() && node->particle_idx == p_idx)
    return;

  double dx = node->com_x - particles_.x[p_idx];
  double dy = node->com_y - particles_.y[p_idx];
  double dz = node->com_z - particles_.z[p_idx];
  double dist_sqr = dx * dx + dy * dy + dz * dz + SOFTENING;
  double dist = std::sqrt(dist_sqr);

  // Opening angle criterion (MAC: Multipole Acceptance Criterion)
  if (node->isLeaf() || (node->bbox.size() / dist < theta_)) {
    // Approximate as a single cluster
    double inv_dist = 1.0 / dist;
    double inv_dist3 = inv_dist * inv_dist * inv_dist;

    double f = G * particles_.mass[p_idx] * node->total_mass * inv_dist3;
    fx += f * dx;
    fy += f * dy;
    fz += f * dz;
  } else {
    // Node is too close, recurse into children
    for (int i = 0; i < 8; ++i) {
      calculateForceOnParticle(p_idx, node->children[i].get(), fx, fy, fz);
    }
  }
}

void BarnesHutTree::computeForces(std::vector<double> &fx,
                                  std::vector<double> &fy,
                                  std::vector<double> &fz) const {
  if (!root_)
    return;

// Can be parallelized safely because we are reading the immutable tree
#pragma omp parallel for schedule(dynamic, 64)
  for (std::size_t i = 0; i < particles_.size; ++i) {
    double current_fx = 0.0;
    double current_fy = 0.0;
    double current_fz = 0.0;

    calculateForceOnParticle(static_cast<int>(i), root_.get(), current_fx,
                             current_fy, current_fz);

    fx[i] = current_fx;
    fy[i] = current_fy;
    fz[i] = current_fz;
  }
}

} // namespace nbody
