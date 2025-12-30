#!/usr/bin/env python3
"""
plot_nbody.py

読み込む CSV の想定形式:
  step,b1_x,b1_y[,b1_z],b2_x,b2_y[,b2_z],...

実行例:
  python3 plot_nbody.py            # nbody_output.csv を読みプロット
  python3 plot_nbody.py --csv output.csv --out traj.png

依存: matplotlib, numpy
"""
import argparse
import csv
import math
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import os


def read_csv(path):
    with open(path, newline='') as f:
        reader = csv.reader(f)
        header = next(reader)
        cols = header[1:]
        # Determine coordinate suffixes per body (x,y[,z])
        bodies = defaultdict(lambda: {})
        for idx, name in enumerate(cols):
            # name example: b1_x
            if '_' in name:
                bname, comp = name.rsplit('_', 1)
            else:
                # fallback
                bname = name
                comp = str(idx)
            bodies[bname][comp] = idx

        # Initialize arrays
        steps = []
        data = {b: defaultdict(list) for b in bodies}

        for row in reader:
            steps.append(int(row[0]))
            for bname, compmap in bodies.items():
                for comp, col_idx in compmap.items():
                    val = float(row[1 + col_idx])
                    data[bname][comp].append(val)

    return steps, data


def plot_trajectories(data, outpath=None, show=True):
    # Detect 2D or 3D
    sample = next(iter(data.values()))
    has_z = 'z' in sample

    plt.figure(figsize=(8, 6))
    ax = plt.gca()

    colors = plt.cm.tab10
    for i, (bname, comps) in enumerate(sorted(data.items())):
        x = np.array(comps.get('x') or comps.get('X'))
        y = np.array(comps.get('y') or comps.get('Y'))
        if x.size == 0 or y.size == 0:
            continue
        ax.plot(x, y, '-', color=colors(i % 10), label=bname)
        ax.plot(x[0], y[0], 'o', color=colors(i % 10), markersize=6)  # start
        ax.plot(x[-1], y[-1], 's', color=colors(i % 10), markersize=6)  # end

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Particle trajectories')
    ax.legend()
    ax.set_aspect('equal', adjustable='datalim')
    plt.tight_layout()

    if outpath:
        plt.savefig(outpath, dpi=200)
        print(f'Saved plot to {outpath}')
    if show:
        plt.show()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--csv', default='nbody_output.csv', help='CSV file to read')
    p.add_argument('--out', default='nbody_trajectories.png', help='Output image path')
    p.add_argument('--no-show', dest='show', action='store_false', help="Don't show interactive window")
    args = p.parse_args()

    if not os.path.exists(args.csv):
        print(f'CSV file not found: {args.csv}')
        return

    steps, data = read_csv(args.csv)
    if len(steps) == 0:
        print('No data in CSV')
        return

    plot_trajectories(data, outpath=args.out, show=args.show)


if __name__ == '__main__':
    main()
