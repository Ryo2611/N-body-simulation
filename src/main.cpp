#include "nbody.hpp"
#include <chrono>
#include <iostream>

int main() {
  const std::size_t NUM_PARTICLES = 100000;
  const int NUM_STEPS = 10;
  const double DT = 0.01;

  std::cout << "Starting C++ N-Body Simulation" << std::endl;
  std::cout << "Particles: " << NUM_PARTICLES << std::endl;
  std::cout << "Steps: " << NUM_STEPS << std::endl;

  nbody::Simulator sim(NUM_PARTICLES);

  auto start_time = std::chrono::high_resolution_clock::now();

  for (int step = 0; step < NUM_STEPS; ++step) {
    sim.step(DT);
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> duration = end_time - start_time;
  double ms = duration.count();

  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Simulation completed in " << ms << " ms." << std::endl;
  std::cout << "Average step time: " << ms / NUM_STEPS << " ms." << std::endl;

  return 0;
}
