#pragma once
#include <cstddef>
#include <vector>

namespace nbody {

inline constexpr double G = 1.0;          // Normalized gravitational constant
inline constexpr double SOFTENING = 1e-9; // To prevent division by zero

struct Particles {
  std::size_t size = 0;
  std::vector<double> x, y, z;
  std::vector<double> vx, vy, vz;
  std::vector<double> mass;
};

class Simulator {
public:
  explicit Simulator(std::size_t num_particles, double target_mass = 1.0);

  // Disable copy for safety
  Simulator(const Simulator &) = delete;
  Simulator &operator=(const Simulator &) = delete;

  void step(double dt);

  // For testing and reading state
  const Particles &getParticles() const;

private:
  Particles particles_;
  void computeForces(std::vector<double> &fx, std::vector<double> &fy,
                     std::vector<double> &fz);
};

} // namespace nbody
