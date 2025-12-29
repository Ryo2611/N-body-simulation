#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

// Constants
const double G = 1.0; // Normalized Gravitational Constant for simplicity
const double DT = 0.01; // Time step
const int STEPS = 1000; // Total simulation steps

struct Point3D {
    double x, y, z;
};

class Body {
public:
    Point3D pos;
    Point3D vel;
    Point3D acc;
    double mass;

    Body(double x, double y, double z, double vx, double vy, double vz, double m) {
        pos = {x, y, z};
        vel = {vx, vy, vz};
        acc = {0.0, 0.0, 0.0};
        mass = m;
    }
};

// Function to compute forces and update accelerations
// Complexity: O(N^2)
void computeAccelerations(std::vector<Body>& bodies) {
    size_t N = bodies.size();

    // Reset accelerations
    for (auto& b : bodies) {
        b.acc = {0.0, 0.0, 0.0};
    }

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            double dx = bodies[j].pos.x - bodies[i].pos.x;
            double dy = bodies[j].pos.y - bodies[i].pos.y;
            double dz = bodies[j].pos.z - bodies[i].pos.z;

            double distSq = dx*dx + dy*dy + dz*dz + 1e-9; // softening to avoid division by zero
            double dist = std::sqrt(distSq);
            double f = (G * bodies[i].mass * bodies[j].mass) / (distSq * dist);

            // Newton's 3rd Law: Force on i is equal and opposite to force on j
            double fx = f * dx;
            double fy = f * dy;
            double fz = f * dz;

            bodies[i].acc.x += fx / bodies[i].mass;
            bodies[i].acc.y += fy / bodies[i].mass;
            bodies[i].acc.z += fz / bodies[i].mass;

            bodies[j].acc.x -= fx / bodies[j].mass;
            bodies[j].acc.y -= fy / bodies[j].mass;
            bodies[j].acc.z -= fz / bodies[j].mass;
        }
    }
}

// Velocity Verlet Integration
// 1. v(t + 0.5dt) = v(t) + 0.5 * a(t) * dt
// 2. x(t + dt) = x(t) + v(t + 0.5dt) * dt
// 3. Calculate new a(t + dt)
// 4. v(t + dt) = v(t + 0.5dt) + 0.5 * a(t + dt) * dt
void updatePositions(std::vector<Body>& bodies) {
    for (auto& b : bodies) {
        b.vel.x += 0.5 * b.acc.x * DT;
        b.vel.y += 0.5 * b.acc.y * DT;
        b.vel.z += 0.5 * b.acc.z * DT;

        b.pos.x += b.vel.x * DT;
        b.pos.y += b.vel.y * DT;
        b.pos.z += b.vel.z * DT;
    }
}

void updateVelocities(std::vector<Body>& bodies) {
    for (auto& b : bodies) {
        b.vel.x += 0.5 * b.acc.x * DT;
        b.vel.y += 0.5 * b.acc.y * DT;
        b.vel.z += 0.5 * b.acc.z * DT;
    }
}

int main() {
    // Initialize standard N-body setup (e.g., Figure-8 or random)
    std::vector<Body> bodies;
    
    // Example: 3-Body Problem (Random initialization for demo)
    bodies.emplace_back(1.0, 0.0, 0.0,  0.0, 0.5, 0.0,  1.0); // Body 1
    bodies.emplace_back(-1.0, 0.0, 0.0, 0.0, -0.5, 0.0, 1.0); // Body 2
    bodies.emplace_back(0.0, 1.0, 0.0, -0.5, 0.0, 0.0,  1.0); // Body 3

    std::ofstream file("nbody_output.csv");
    file << "step,b1_x,b1_y,b2_x,b2_y,b3_x,b3_y\n";

    // Main Loop
    computeAccelerations(bodies); // Initial acceleration calculation

    for (int t = 0; t < STEPS; ++t) {
        updatePositions(bodies);
        computeAccelerations(bodies); // Calculate a(t+dt)
        updateVelocities(bodies);

        // Output data for visualization
        file << t << ",";
        for (size_t i = 0; i < bodies.size(); ++i) {
             file << bodies[i].pos.x << "," << bodies[i].pos.y;
             if(i < bodies.size() - 1) file << ",";
        }
        file << "\n";
    }

    file.close();
    std::cout << "Simulation complete. Data written to nbody_output.csv" << std::endl;
    return 0;
}