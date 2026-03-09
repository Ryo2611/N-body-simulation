# High-Performance N-Body Simulation (C++17)

A highly optimized N-body gravity simulator, designed as a comprehensive portfolio piece demonstrating low-level performance engineering and algorithmic complexity management. 

This project iteratively optimizes a naive $O(N^2)$ algorithm down to $O(N \log N)$ complexity, leveraging cache locality, hardware vectorization, and multi-core parallelism. The engine scales effortlessly to compute the gravitational interactions of $100,000$ particles in under 100 milliseconds per step on an M-series Apple Silicon chip.

## 🚀 Optimization Journey & Benchmarks

The simulation underwent four major optimization phases. All benchmarks were performed on Apple Silicon (`arm64`) using `AppleClang` with `-O3 -mcpu=apple-m1` flags.

### 1. The Naive Baseline ($O(N^2)$)
- **Data Structure**: Array of Structures (AoS) - `std::vector<Particle>` (where `Particle` contains `pos`, `vel`, `mass`).
- **Algorithm**: Direct all-pairs summation.
- **Performance**: ~51 ms per step for $N=2,000$. Computing macroscopic scales ($N > 10,000$) is computationally implausible.

### 2. Cache-Locality & Memory Layout (Structure of Arrays)
- **Data Structure**: Structure of Arrays (SoA) - `Particles` struct containing separate contiguous `std::vector<double>` arrays for `x, y, z`, `vx, vy, vz`, and `mass`.
- **Why**: Drastically reduces L1/L2 cache misses. When computing forces along the X-axis, the CPU fetches contiguous X-coordinates into the cache line without loading unused Y/Z or velocity data.

### 3. SIMD Vectorization & Multi-threading (ARM NEON + OpenMP)
- **SIMD**: Integrated **ARM NEON Intrinsics** (`<arm_neon.h>`) to process two `double` precision floats simultaneously in 128-bit registers (`float64x2_t`). Since true `rsqrt` for `f64` isn't natively supported on NEON, a manual fallback using multiplication/sqrt was applied to the vector. 
- **Multi-threading**: Utilized `#pragma omp parallel for` (via Homebrew's `libomp`) to distribute the outer force-calculation loop across all independent physical cores.
- **Performance**: The combination of SoA, SIMD, and OpenMP brought the execution time down to **~1.7 ms per step** for $N=2,000$ (A **~30x speedup** over the baseline).

### 4. Algorithmic Optimization (Barnes-Hut QuadTree)
- **Algorithm**: The mechanical limitation of $O(N^2)$ cannot be bypassed solely with hardware. By partitioning the 2D space recursively into a QuadTree, distant groups of particles are approximated as a single massive center of mass.
- **Complexity**: $O(N \log N)$.
- **Performance**: The simulation can now handle **$100,000$ particles at ~95 ms per step**. Running 100,000 particles on the naive engine would take hours per frame.

## 🛠 Project Structure & Tools
- **Language**: C++17 (Strictly enforced `-Wall -Werror -pedantic`)
- **Build System**: CMake (3.14+)
- **Testing**: GoogleTest (integrated via `FetchContent`) to guarantee physical invariants like **Momentum Conservation**.

## 💻 Building and Running

### Prerequisites (macOS)
```bash
brew install cmake libomp
```

### Build Instructions
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
```

### Run the Simulation
```bash
./nbody_sim
```

### Run the Test Suite
```bash
./nbody_tests
```
