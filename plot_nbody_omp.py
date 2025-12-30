#!/usr/bin/env python3
"""
plot_nbody_omp.py

`nbody_omp` の実行ログ（標準出力をファイルに保存したもの）を解析し、
粒子数（particles）に対する性能（Million pairs/sec）をプロットします。

使い方:
  python3 plot_nbody_omp.py               # デフォルトで nbody_omp.log を解析
  python3 plot_nbody_omp.py logs/*.log    # 複数ログを指定
  python3 plot_nbody_omp.py --out perf.png

ログ中の以下の行を解析します:
  Running Simulation with <N> particles...
  Time elapsed: <sec> seconds.
  Performance: <value> Million pairs/sec

依存: matplotlib, numpy
"""
import re
import argparse
import glob
import os
from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np


Entry = namedtuple('Entry', ['file', 'particles', 'time', 'perf', 'threads'])


def parse_log(path):
    particles = None
    time = None
    perf = None
    threads = None
    with open(path, 'r') as f:
        for line in f:
            if particles is None:
                m = re.search(r'Running Simulation with\s*(\d+)\s*particles', line)
                if m:
                    particles = int(m.group(1))
            if threads is None:
                m = re.search(r'Max Threads available:\s*(\d+)', line)
                if m:
                    threads = int(m.group(1))
            m = re.search(r'Time elapsed:\s*([0-9.]+)\s*seconds', line)
            if m:
                time = float(m.group(1))
            m = re.search(r'Performance:\s*([0-9.]+)\s*Million pairs/sec', line)
            if m:
                perf = float(m.group(1))
    return Entry(file=os.path.basename(path), particles=particles, time=time, perf=perf, threads=threads)


def collect_entries(paths):
    entries = []
    for p in paths:
        for f in glob.glob(p):
            try:
                e = parse_log(f)
                entries.append(e)
            except Exception as ex:
                print(f'Warning: failed parse {f}: {ex}')
    # Filter out entries without performance
    entries = [e for e in entries if e.perf is not None]
    return entries


def plot(entries, out='nbody_omp_performance.png', show=True):
    if not entries:
        print('No performance data to plot.')
        return

    particles = np.array([e.particles for e in entries], dtype=float)
    perf = np.array([e.perf for e in entries], dtype=float)
    labels = [e.file for e in entries]

    # Sort by particles
    idx = np.argsort(particles)
    particles = particles[idx]
    perf = perf[idx]
    labels = [labels[i] for i in idx]

    plt.figure(figsize=(7,5))
    plt.plot(particles, perf, '-o')
    for x,y,l in zip(particles, perf, labels):
        plt.annotate(l, (x,y), textcoords='offset points', xytext=(5,3), fontsize=8)
    plt.xlabel('Particles')
    plt.ylabel('Performance (Million pairs/sec)')
    plt.title('nbody_omp performance')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    print(f'Saved plot to {out}')
    if show:
        plt.show()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('paths', nargs='*', default=['nbody_omp.log'], help='Log file paths or glob patterns')
    p.add_argument('--out', default='nbody_omp_performance.png', help='Output image path')
    p.add_argument('--no-show', dest='show', action='store_false', help="Don't show interactive window")
    args = p.parse_args()

    entries = collect_entries(args.paths)
    if not entries:
        print('No valid log entries found. Check paths or log content.')
        return
    plot(entries, out=args.out, show=args.show)


if __name__ == '__main__':
    main()
