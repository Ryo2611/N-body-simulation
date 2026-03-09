# High-Performance 3D N-Body Simulation (C++17)

![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)
![CMake](https://img.shields.io/badge/CMake-3.14%2B-brightgreen.svg)
![OpenMP](https://img.shields.io/badge/OpenMP-Enabled-orange.svg)
![ARM NEON](https://img.shields.io/badge/SIMD-ARM%20NEON-red.svg)

A highly optimized 3D N-body gravity simulator designed as a comprehensive software engineering and performance algorithm portfolio piece for quantitative finance and high-performance computing (HPC) research.

This project iteratively optimizes a naive $O(N^2)$ algorithm down to an $O(N \log N)$ spatial partitioning complexity. By leveraging cache locality (SoA), hardware vectorization (SIMD), multi-core CPU parallelism (OpenMP), and a 3D Barnes-Hut Octree, the engine scales effortlessly to compute the gravitational interactions of **100,000 particles in 3D-space in ~255 milliseconds per step** on an M-series Apple Silicon chip.

---

## 🚀 Optimization Journey & Benchmarks

This simulation demonstrates addressing massive computational complexity through a low-level systems engineering approach. All benchmarks were recorded on an Apple Silicon (`arm64`) architecture with `-O3 -mcpu=apple-m1` compiler flags enabled.

### 1. Naive Baseline Implementation ($O(N^2)$)
- **Data Structure**: Array of Structures (AoS) - `std::vector<Particle>`
- **Algorithm**: Brute-force all-pairs gravity calculation.
- **Performance**: Computing $N=2,000$ particles took **~51 ms** per step. Scaling to 10,000+ particles was computationally impractical.

### 2. Cache Locality Improvement (Structure of Arrays - SoA)
- **Data Structure**: Structure of Arrays (SoA) - Refactored `Particles` to hold contiguous vectors (`std::vector<double>`) for `x`, `y`, `z`, `vx`, `vy`, and `vz`.
- **Impact**: Drastically reduced L1/L2 cache misses. When computing forces along the X-axis, the CPU can load solely position vectors into cache lines without unpacking unnecessary interleaved data (like velocities or Y/Z components).

### 3. SIMD Vectorization & Multi-threading (ARM NEON + OpenMP)
- **SIMD (Single Instruction, Multiple Data)**: Integrated **ARM NEON Intrinsics** (`<arm_neon.h>`) to process two double-precision floats (`float64x2_t`) simultaneously in a single 128-bit clock cycle.
- **Multi-threading**: Utilized `#pragma omp parallel for` via Homebrew's `libomp` to distribute the outermost gravity evaluation loops across all physical CPU cores dynamically.
- **Performance**: The synergy of SoA, SIMD, and OpenMP reduced the step time for $N=2,000$ to **~1.7 ms** **(a ~30x geometric speedup from baseline)**.

### 4. Algorithmic Scaling (3D Barnes-Hut Octree)
- **Algorithm**: To break through hardware limitations, we implemented an **Octree**, recursively dividing the 3D space into 8 octants. Distant particle clusters are approximated as a single massive center of mass using the Multipole Acceptance Criterion ($\theta < 0.5$).
- **Complexity**: Slashed algorithmic cost from $O(N^2)$ down to **$O(N \log N)$**.
- **Performance**: The engine can now simulate **100,000 particles (in 3D)** at **~255 ms** per step, bringing calculations that would take hours on the naive engine down to near real-time timescales.

---

## 🛠 Project Architecture

- **Language**: Standard C++17
- **Compiler Config**: Strict clean builds targeting `-Wall -Wextra -Werror -pedantic`.
- **Build System**: CMake (3.14+)
- **Unit Testing**: GoogleTest (configured via `FetchContent`). Ensures physical integrity (e.g., verifying Momentum Conservation within the acceptable tolerance bounds of the algorithmic MAC approximation).

```text
├── CMakeLists.txt
├── include/
│   ├── nbody.hpp      # Simulator engine definitions
│   └── oct_tree.hpp   # 3D Barnes-Hut octree definitions
├── src/
│   ├── main.cpp       # Entry point & benchmark suite
│   ├── nbody.cpp      # Integrator & core simulation body
│   └── oct_tree.cpp   # Octree construction & force approximation algorithms
└── tests/
    └── test_nbody.cpp # GoogleTest Suite
```

---

## 💻 Quick Start Guide

### Prerequisites (macOS / Apple Silicon)
Ensure OpenMP and CMake are installed via Homebrew.
```bash
brew install cmake libomp
```

### Build Instructions
Generate the highly optimized production build output.
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
```

### Run the Benchmark Simulation
Execute the core logic directly via the terminal.
```bash
./nbody_sim
```

### Run the Test Suite
Validate the physical integrity of the N-Body numerical integrator.
```bash
./nbody_tests
```
