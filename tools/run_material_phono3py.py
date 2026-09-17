#!/usr/bin/env python3
"""Compute intrinsic three-phonon RTA N/U data from a supplied material pack.

Run in an empty job directory using the phono3py Python environment. Original
force constants are read unchanged; no force-constant refitting is performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

import h5py
import numpy as np
import phono3py
import phonopy
from phonopy.interface.vasp import write_vasp
from phonopy.phonon.grid import get_ir_grid_points


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--job', type=Path, required=True)
    p.add_argument('--mesh', type=int, nargs=3, required=True)
    p.add_argument('--temperatures', type=float, nargs='+', required=True)
    p.add_argument('--resume', action='store_true', help='Validate and reuse completed per-q checkpoints')
    p.add_argument('--exact-zero-mask', action='store_true', help='Skip only exactly zero FC3 atom-triplet tensors')
    args = p.parse_args()
    source, job = args.source.resolve(), args.job.resolve()
    job.mkdir(parents=True, exist_ok=True)
    if (job / 'manifest.json').exists() and not args.resume:
        raise RuntimeError('Job already initialized; use a new directory to preserve results.')
    inp = job / 'input'
    inp.mkdir(exist_ok=True)
    hashes = {}
    for name in ('POSCAR', 'phono3py.yaml', 'phono3py_disp.yaml', 'fc2.hdf5', 'fc3.hdf5'):
        if not args.resume:
            shutil.copy2(source / name, inp / name)
        hashes[name] = digest(inp / name)
    manifest = {
        'source': str(source), 'input_sha256': hashes,
        'mesh': args.mesh, 'temperatures_K': args.temperatures,
        'phono3py': phono3py.__version__, 'phonopy': phonopy.__version__,
        'method': 'RTA', 'scattering': 'intrinsic three-phonon N and U only',
        'integration': 'tetrahedron', 'backend': 'C',
        'symmetry_tolerance_A': 1e-5, 'cutoff_frequency_THz': 1e-4,
        'make_r0_average': True, 'symmetrize_fc3q': False,
        'force_constants_modified': False,
        'gamma_convention': 'HWHM THz; inverse lifetime in ps^-1 = 4*pi*gamma',
        'environment': {key: os.environ.get(key) for key in
                        ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')},
        'status': 'running', 'started_unix': time.time(),
    }
    manifest_path = job / 'manifest.json'
    if args.resume:
        previous = json.loads(manifest_path.read_text())
        for key in ('source', 'input_sha256', 'mesh', 'temperatures_K', 'phono3py', 'phonopy'):
            if previous[key] != manifest[key]:
                raise RuntimeError(f'Resume input mismatch: {key}')
        manifest['started_unix'] = previous['started_unix']
        manifest['resumed_unix'] = time.time()
    manifest['pid'] = os.getpid()
    manifest['exact_zero_fc3_mask'] = args.exact_zero_mask
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    os.chdir(job)
    try:
        ph = phono3py.load(
            inp / 'phono3py.yaml', fc2_filename=inp / 'fc2.hdf5',
            fc3_filename=inp / 'fc3.hdf5', produce_fc=False,
            is_nac=False, make_r0_average=True, lang='C', log_level=1,
        )
        if args.exact_zero_mask:
            ph.fc3_nonzero_indices = np.array(np.any(ph.fc3 != 0, axis=(-3, -2, -1)), dtype='byte', order='C')
            manifest['fc3_nonzero_atom_triplet_fraction'] = float(ph.fc3_nonzero_indices.mean())
        write_vasp('POSCAR.primitive', ph.primitive)
        ph.mesh_numbers = args.mesh
        ph.init_phph_interaction()
        ph.run_phonon_solver()
        grid = ph.grid
        ir_grg, ir_weight, ir_map = get_ir_grid_points(grid)
        freq, eigvec, _ = ph.get_phonon_data()
        bz = grid.grg2bzg
        q = grid.addresses[bz] @ grid.QDinv.T
        with h5py.File('harmonic-grid.hdf5', 'w') as f:
            for key, value in {
                'qpoint': q, 'frequency': freq[bz], 'mesh': args.mesh,
                'grg2bzg': bz, 'ir_grg': ir_grg, 'ir_weight': ir_weight,
                'ir_map': ir_map, 'primitive_lattice': ph.primitive.cell,
            }.items():
                f.create_dataset(key, data=value)
        manifest['primitive_volume_A3'] = ph.primitive.volume
        manifest['num_ir_grid_points'] = len(ir_grg)
        manifest['num_full_grid_points'] = len(bz)
        manifest['minimum_frequency_THz'] = float(freq[bz].min())
        manifest['negative_frequency_modes_below_minus_0_01_THz'] = int((freq[bz] < -0.01).sum())
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
        print('COMPUTE', json.dumps({k: manifest[k] for k in
              ('mesh', 'temperatures_K', 'num_ir_grid_points', 'minimum_frequency_THz')}), flush=True)
        common = dict(temperatures=args.temperatures, is_LBTE=False, is_N_U=True,
                      is_isotope=False, boundary_mfp=None, log_level=1)
        if args.resume:
            remaining = []
            mesh_str = ''.join(map(str, args.mesh))
            for gp in grid.grg2bzg[ir_grg]:
                checkpoint = job / f'kappa-m{mesh_str}-g{gp}.hdf5'
                if not checkpoint.exists():
                    remaining.append(int(gp))
                    continue
                with h5py.File(checkpoint) as f:
                    np.testing.assert_array_equal(f['temperature'][:], args.temperatures)
                    np.testing.assert_array_equal(f['mesh'][:], args.mesh)
                    np.testing.assert_allclose(f['frequency'][:], freq[gp], atol=5e-6, rtol=1e-7)
                    for key in ('gamma', 'gamma_N', 'gamma_U'):
                        assert f[key].shape == (len(args.temperatures), len(ph.primitive) * 3)
                        assert np.isfinite(f[key][:]).all() and (f[key][:] >= 0).all()
                    np.testing.assert_allclose(f['gamma_N'][:] + f['gamma_U'][:], f['gamma'][:], atol=1e-14, rtol=1e-7)
            print(f'RESUME: {len(ir_grg) - len(remaining)} complete; {len(remaining)} remaining.', flush=True)
            if remaining:
                ph.run_thermal_conductivity(grid_points=remaining, write_gamma=True, **common)
            ph.run_thermal_conductivity(read_gamma=True, write_kappa=True, **common)
        else:
            ph.run_thermal_conductivity(write_gamma=True, write_kappa=True, **common)
        manifest['status'] = 'calculated'
        manifest['finished_unix'] = time.time()
        manifest['elapsed_seconds'] = manifest['finished_unix'] - manifest['started_unix']
    except BaseException as exc:
        manifest['status'] = 'failed'
        manifest['error'] = repr(exc)
        raise
    finally:
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
