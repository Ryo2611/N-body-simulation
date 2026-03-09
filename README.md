# High-Performance 3D N-Body Simulation (C++17)

![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)
![CMake](https://img.shields.io/badge/CMake-3.14%2B-brightgreen.svg)
![OpenMP](https://img.shields.io/badge/OpenMP-Enabled-orange.svg)
![ARM NEON](https://img.shields.io/badge/SIMD-ARM%20NEON-red.svg)

A highly optimized 3D N-body gravity simulator designed as a comprehensive software engineering and performance algorithm portfolio piece for quantitative finance and high-performance computing (HPC) research.

This project iteratively optimizes a naive $O(N^2)$ algorithm down to an $O(N \log N)$ spatial partitioning complexity. By leveraging cache locality (SoA), hardware vectorization (SIMD), multi-core CPU parallelism (OpenMP), and a 3D Barnes-Hut Octree, the engine scales effortlessly to compute the gravitational interactions of **100,000 particles in 3D-space in ~255 milliseconds per step** on an M-series Apple Silicon chip.

---

## 🚀 最適化プロセスとベンチマーク結果 (Optimization Journey)

当シミュレーションは、計算コストが膨大となるN体問題をソフトウェアエンジニアリングおよび低レイヤー最適化の観点から高速化したものです。すべてのベンチマークはApple Silicon (`arm64`) 環境（`-O3 -mcpu=apple-m1` フラグ有効化）で計測されています。

### 1. ナイーブなベースライン実装 ($O(N^2)$)
- **データ構造**: Array of Structures (AoS) - `std::vector<Particle>`
- **アルゴリズム**: すべての粒子間の重力を直接計算する総当たり方式
- **パフォーマンス**: $N=2,000$ 粒子の場合、1ステップあたり **~51 ms**。10,000粒子以上の大規模シミュレーションは実用不可能な計算量。

### 2. キャッシュ局所性の向上 (Structure of Arrays - SoA)
- **データ構造**: 要素ごとの配列 (SoA) - `Particles` 構造体内に `x, y, z` 座標や `vx, vy, vz` 速度を個別の等間隔な配列(`std::vector<double>`)として保存
- **目的と効果**: L1/L2キャッシュミスを劇的に削減。例えば、X座標の力を計算する際、CPUは使用しないY/Zや速度のデータをメモリからフェッチすることなく、必要なベクトルのみをキャッシュラインにロードできます。

### 3. SIMD ベクトル化 & マルチスレッド化 (ARM NEON + OpenMP)
- **SIMD (単一命令複数データ処理)**: **ARM NEON Intrinsics** (`<arm_neon.h>`) を導入し、128-bitレジスタ（`float64x2_t`）で2つの倍精度浮動小数点数（`double`）を1クロックで同時処理。
- **マルチスレッド**: `#pragma omp parallel for` を使用して、重力計算の最外ループをCPUの物理コア全体（Homebrewの `libomp` 使用）にスレッド分散。
- **パフォーマンス**: SoA・SIMD・OpenMPの相乗効果により、$N=2,000$ の処理速度が **~1.7 ms** に短縮 **(ベースライン比 約30倍の高速化)**。

### 4. アルゴリズム最適化 (3D Barnes-Hut Octree)
- **アルゴリズム**: ハードウェア限界を突破するため、3D空間を再帰的に8分割する **Octree（八分木）** を実装。遠方にある粒子群の重心（Center of Mass）を一つの巨大な質点として近似計算（Multipole Acceptance Criterion: $\theta < 0.5$）。
- **計算量**: 空間分割木により $O(N^2)$ から **$O(N \log N)$** に削減。
- **パフォーマンス**: **100,000粒子 (3D空間)** を1ステップ **~255 ms** で処理可能に。$O(N^2)$ エンジンでは数時間かかる計算をリアルタイムに近い速度で処理できる次元へとスケーリングを証明しました。

---

## 🛠 プロジェクト構成 (Architecture & Tools)
- **言語**: C++17
- **コンパイル設定**: 厳格な警告基準によるクリーンビルド (`-Wall -Wextra -Werror -pedantic`)
- **ビルドシステム**: CMake (3.14+)
- **ユニットテスト**: GoogleTest (`FetchContent` にて自動構成)。「運動量保存則」などの物理的整合性が保たれているかを自動的にテストします（Barnes-Hut のマクロレベルの近似誤差許容も考慮済）。

```text
├── CMakeLists.txt
├── include/
│   ├── nbody.hpp      # Simulator engine definitions
│   └── oct_tree.hpp   # 3D Barnes-Hut octree definitions
├── src/
│   ├── main.cpp       # Entry point
│   ├── nbody.cpp      # Integrator & Simulation body
│   └── oct_tree.cpp   # Octree construction & forces calculation
└── tests/
    └── test_nbody.cpp # GoogleTest Suite
```

---

## 💻 ビルドと実行方法 (Quick Start)

### 前提要件 (macOS / Apple Silicon 向け)
OpenMPライブラリをインストールします。
```bash
brew install cmake libomp
```

### ビルド手順
CMakeを使って最適化ビルド(`Release`)を生成します。
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
```

### シミュレーションの実行
```bash
./nbody_sim
```

### ユニットテストの実行
シミュレーションの整合性チェックを行います。
```bash
./nbody_tests
```
