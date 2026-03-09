#include "nbody.hpp"
#include <algorithm>
#include <cmath>
#include <random>

#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif

#include "oct_tree.hpp"

namespace nbody {

// Removed G and SOFTENING as they are now in the header

Simulator::Simulator(std::size_t num_particles, double target_mass) {
  particles_.size = num_particles;

  particles_.x.reserve(num_particles);
  particles_.y.reserve(num_particles);
  particles_.z.reserve(num_particles);
  particles_.vx.reserve(num_particles);
  particles_.vy.reserve(num_particles);
  particles_.vz.reserve(num_particles);
  particles_.mass.reserve(num_particles);

  std::mt19937_64 rng(42); // Fixed seed for reproducibility
  std::uniform_real_distribution<double> pos_dist(-1.0, 1.0);
  std::uniform_real_distribution<double> vel_dist(-0.1, 0.1);

  for (std::size_t i = 0; i < num_particles; ++i) {
    particles_.x.push_back(pos_dist(rng));
    particles_.y.push_back(pos_dist(rng));
    particles_.z.push_back(pos_dist(rng)); // Randomize Z now
    particles_.vx.push_back(vel_dist(rng));
    particles_.vy.push_back(vel_dist(rng));
    particles_.vz.push_back(vel_dist(rng));
    particles_.mass.push_back(target_mass);
  }
}

void Simulator::computeForces(std::vector<double> &fx, std::vector<double> &fy,
                              std::vector<double> &fz) {
  std::fill(fx.begin(), fx.end(), 0.0);
  std::fill(fy.begin(), fy.end(), 0.0);
  std::fill(fz.begin(), fz.end(), 0.0);

  // Build the 3D Barnes-Hut Tree for the current step
  // O(N log N) construction
  BarnesHutTree tree(particles_, 0.5); // theta = 0.5 (standard MAC)

  // O(N log N) force calculation across all 3 spatial dimensions
  tree.computeForces(fx, fy, fz);
}

void Simulator::step(double dt) {
  std::size_t n = particles_.size;
  std::vector<double> fx(n, 0.0), fy(n, 0.0), fz(n, 0.0);

  computeForces(fx, fy, fz);

// Semi-implicit Euler integration
#pragma omp parallel for
  for (std::size_t i = 0; i < n; ++i) {
    particles_.vx[i] += (fx[i] / particles_.mass[i]) * dt;
    particles_.vy[i] += (fy[i] / particles_.mass[i]) * dt;
    particles_.vz[i] += (fz[i] / particles_.mass[i]) * dt;

    particles_.x[i] += particles_.vx[i] * dt;
    particles_.y[i] += particles_.vy[i] * dt;
    particles_.z[i] += particles_.vz[i] * dt;
  }
}

const Particles &Simulator::getParticles() const { return particles_; }

} // namespace nbody
