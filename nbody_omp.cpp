#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <omp.h> // Header for OpenMP
// Constants
const double G = 1.0;
const double DT = 0.01;
const int STEPS = 10;        // reduced for quick test
const int N_PARTICLES = 100; // reduced for quick test

struct Point3D
{
    double x, y, z;
};

class Body
{
public:
    Point3D pos;
    Point3D vel;
    Point3D acc;
    double mass;

    // Use default constructor for array initialization
    Body() : pos({0, 0, 0}), vel({0, 0, 0}), acc({0, 0, 0}), mass(1.0) {}
};

// Function to compute forces
void computeAccelerations(std::vector<Body> &bodies)
{
    size_t N = bodies.size();

// The "schedule(static)" clause divides the loop into equal chunks for each core.
// We parallelize the OUTER loop (i). Each thread gets a private 'i'.
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < N; ++i)
    {

        // Use local variables to accumulate force to avoid accessing cache constantly
        double ax = 0.0;
        double ay = 0.0;
        double az = 0.0;

        // Inner loop: We must iterate all j to avoid Race Conditions on writing to j.
        for (size_t j = 0; j < N; ++j)
        {
            if (i == j)
                continue; // Skip self

            double dx = bodies[j].pos.x - bodies[i].pos.x;
            double dy = bodies[j].pos.y - bodies[i].pos.y;
            double dz = bodies[j].pos.z - bodies[i].pos.z;

            // Softening parameter prevents division by zero if particles overlap
            double distSq = dx * dx + dy * dy + dz * dz + 1e-9;
            double dist = std::sqrt(distSq);

            // Force magnitude: F = G * m_i * m_j / r^2
            // Acceleration: a = F / m_i = G * m_j / r^2
            // Vector component: a_x = a * (dx / r) = G * m_j * dx / r^3
            double f_over_r = (G * bodies[j].mass) / (distSq * dist);

            ax += f_over_r * dx;
            ay += f_over_r * dy;
            az += f_over_r * dz;
        }

        // Write back to main memory once per particle
        bodies[i].acc.x = ax;
        bodies[i].acc.y = ay;
        bodies[i].acc.z = az;
    }
}

void updatePositions(std::vector<Body> &bodies)
{
// These loops are "embarrassingly parallel" - no dependencies
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < bodies.size(); ++i)
    {
        bodies[i].vel.x += 0.5 * bodies[i].acc.x * DT;
        bodies[i].vel.y += 0.5 * bodies[i].acc.y * DT;
        bodies[i].vel.z += 0.5 * bodies[i].acc.z * DT;

        bodies[i].pos.x += bodies[i].vel.x * DT;
        bodies[i].pos.y += bodies[i].vel.y * DT;
        bodies[i].pos.z += bodies[i].vel.z * DT;
    }
}

void updateVelocities(std::vector<Body> &bodies)
{
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < bodies.size(); ++i)
    {
        bodies[i].vel.x += 0.5 * bodies[i].acc.x * DT;
        bodies[i].vel.y += 0.5 * bodies[i].acc.y * DT;
        bodies[i].vel.z += 0.5 * bodies[i].acc.z * DT;
    }
}

int main()
{
    // Generate random particles for a stress test
    std::vector<Body> bodies(N_PARTICLES);
    for (auto &b : bodies)
    {
        b.pos = {(double)(rand() % 100), (double)(rand() % 100), (double)(rand() % 100)};
        b.mass = 1.0;
    }

    std::cout << "Running Simulation with " << N_PARTICLES << " particles..." << std::endl;
    std::cout << "Max Threads available: " << omp_get_max_threads() << std::endl;

    double start_time = omp_get_wtime(); // Start Timer

    computeAccelerations(bodies);

    for (int t = 0; t < STEPS; ++t)
    {
        updatePositions(bodies);
        computeAccelerations(bodies);
        updateVelocities(bodies);
    }

    double end_time = omp_get_wtime(); // End Timer

    std::cout << "Simulation complete." << std::endl;
    std::cout << "Time elapsed: " << (end_time - start_time) << " seconds." << std::endl;

    // Performance Metric: Pairs per second
    double total_interactions = (double)N_PARTICLES * N_PARTICLES * STEPS;
    std::cout << "Performance: " << (total_interactions / (end_time - start_time)) / 1e6 << " Million pairs/sec" << std::endl;

    return 0;
}