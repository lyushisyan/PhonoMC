#!/usr/bin/env python3
"""Create an explicitly synthetic Callaway mode bank, never a material prediction."""
import argparse
from pathlib import Path

import h5py
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "output/callaway_demo/material")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if (args.out / "kappa-fbz.hdf5").exists() or (args.out / "POSCAR").exists():
        parser.error("Output already contains material files; choose a new --out directory.")
    (args.out / "POSCAR").write_text("SYNTHETIC Callaway demonstration, not a real material\n1\n5 0 0\n0 5 0\n0 0 5\nX\n1\nDirect\n0 0 0\n")
    normal = np.array([[[0.03, 0.05, 0.09], [0.03, 0.05, 0.09], [0.03, 0.05, 0.09]]])
    umklapp = np.array([[[0.002, 0.004, 0.01], [0.002, 0.004, 0.01], [0.002, 0.004, 0.01]]])
    with h5py.File(args.out / "kappa-fbz.hdf5", "x") as f:
        f.attrs["synthetic_demo"] = True
        f.attrs["description"] = "Artificial positive-frequency mode bank; no first-principles interpretation."
        for key, value in {
            "mesh": [3, 1, 1], "weight": [1, 1, 1], "temperature": [300.],
            "qpoint": [[0., 0., 0.], [1/3, 0., 0.], [-1/3, 0., 0.]],
            "frequency": [[1., 2., 3.]] * 3,
            "group_velocity": [[[0., 0., 0.]] * 3, [[10., 0., 0.], [8., 0., 0.], [6., 0., 0.]],
                               [[-10., 0., 0.], [-8., 0., 0.], [-6., 0., 0.]]],
            "gamma": normal + umklapp, "gamma_N": normal, "gamma_U": umklapp,
        }.items():
            f.create_dataset(key, data=value)
    print(f"Synthetic demonstration material written to {args.out.resolve()}")


if __name__ == "__main__":
    main()
