#include "nbody.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  bool visualize = false;
  if (argc > 1 && std::string(argv[1]) == "--visualize") {
    visualize = true;
  }

  const std::size_t NUM_PARTICLES = visualize ? 1000 : 100000;
  const int NUM_STEPS = visualize ? 200 : 10;
  const double DT = 0.01;

  std::cout << "Starting C++ N-Body Simulation" << std::endl;
  std::cout << "Mode: " << (visualize ? "Visualization" : "Benchmark")
            << std::endl;
  std::cout << "Particles: " << NUM_PARTICLES << std::endl;
  std::cout << "Steps: " << NUM_STEPS << std::endl;

  nbody::Simulator sim(NUM_PARTICLES);
  std::ofstream out;

  if (visualize) {
    out.open("output.csv");
    out << "step,id,x,y,z\n";
  }

  auto start_time = std::chrono::high_resolution_clock::now();

  for (int step = 0; step < NUM_STEPS; ++step) {
    if (visualize) {
      const auto &p = sim.getParticles();
      for (std::size_t i = 0; i < p.size; ++i) {
        out << step << "," << i << "," << p.x[i] << "," << p.y[i] << ","
            << p.z[i] << "\n";
      }
    }

    sim.step(DT);
  }

  if (visualize) {
    out.close();
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> duration = end_time - start_time;
  double ms = duration.count();

  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Simulation completed in " << ms << " ms." << std::endl;
  std::cout << "Average step time: " << ms / NUM_STEPS << " ms." << std::endl;

  return 0;
}
