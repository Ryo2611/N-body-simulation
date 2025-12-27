#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>


const double G = 1.0;
const double DT = 0.01;
const int STEPS = 1000;
const int N_PARTICLES = 100;

struct Point3D {
    double x, y, z;
};

class Particle {
public:
    Point3D pos;
    Point3D vel;
    Point3D acc;
    double mass;
    int id;

    Particle(double m, int i): mass(m), id(i) {
        pos = {0.0, 0.0, 0.0};
        vel = {0.0, 0.0, 0.0};
        acc = {0.0, 0.0, 0.0};
    }
};

class NBodySystem {
private:
    std::vector<Particle> particles;

public:
    NBodySystem(int n) {
        particles.reserve(n);
        initialize_particles(n);
    }

    void initialize_particles(int n) {
        std::mt19937 gen(42);
        std::uniform_real_distribution<double> dis(-10.0, 10.0);
        std::uniform_real_distribution<double> mass_dis(0.1, 5.0);

        for (int i = 0; i < n; ++i) {
            Particle p(mass_dis(gen), i);
            p.pos = {dis(gen) * 10.0, dis(gen) * 10.0, dis(gen) * 10.0};
            particles.push_back(p);
        }
    }

    void update_positions() {
        for (auto& p : particles) {
            p.pos.x += p.vel.x * DT + 0.5 * p.acc.x * DT * DT;
            p.pos.y += p.vel.y + DT + 0.5 * p.acc.y * DT * DT;
            p.pos.z += p.vel.z * DT + 0.5 * p.acc.z * DT * DT;

            p.vel.x += 0.5 * p.acc.x * DT;
            p.vel.y += 0.5 * p.acc.y * DT;
            p.vel.x += 0.5 * p.acc.z * DT;
        }
    }

    void compute_forces() {
        for (auto& p : particles) {
            p.acc = {0.0, 0.0, 0.0};
        }

        for (int i = 0; i < particles.size(); ++i) {
            for (int j = i + 1; j < particles.size(); ++j) {
                double dx = particles[j].pos.x - particles[i].pos.x;
                double dy = particles[j].pos.y - particles[i].pos.y;
                double dz = particles[j].pos.z - particles[i].pos.z;
                double dist_sq = dx * dx + dy * dy + dz * dz + 1e-9;
                double dist = std::sqrt(dist_sq);
                double f_mag = (G * particles[i].mass * particles[j].mass) / dist_sq;

                double fx = f_mag * (dx / dist);
                double fy = f_mag * (dy / dist);
                double fz = f_mag * (dz / dist);

                particles[i].acc.x += fx / particles[i].mass;
                particles[i].acc.y += fy / particles[i].mass;
                particles[i].acc.z += fz / particles[i].mass;

                particles[j].acc.x -= fx / particles[j].mass;
                particles[j].acc.y -= fy / particles[j].mass;
                particles[j].acc.z -= fz / particles[j].mass;
            }
        }
    }

    void update_velocities() {
        for (auto& p : particles){
            p.vel.x += p.acc.x * DT * 0.5;
            p.vel.y += p.acc.y * DT * 0.5;
            p.vel.z += p.acc.z * DT * 0.5;
        }
    }

    void run() {
        std::ofstream file("output.csv");
        file << "t,id,x,y,z\n";

        compute_forces();

        for (int t = 0; t < STEPS; ++t) {
            update_positions();
            compute_forces();
            update_velocities();

            if (t % 10 == 0) {
                file << t * DT << "," << particles[0].id << ","
                    << particles[0].pos.x << "," << particles[0].pos.y << "," << particles[0].pos.z << "\n";
            }
        }
        file.close();
        std::cout << "Simulation complete. Data saved to output.csv" << std::endl;

    }
};


int main()
{
    NBodySystem system(N_PARTICLES);
    system.run();

    return 0;
}