#!/usr/bin/env python3
"""Validate a run_material_phono3py job and export a complete PhonoMC mesh.

Uses the exact GRG -> irreducible mapping saved by phono3py, with one point per
mesh address. Velocities are recomputed from FC2 at every full-grid q point.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil

import h5py
import numpy as np

EXTRA_LINEWIDTHS = ('gamma_isotope', 'gamma_impurity', 'gamma_defect')


def expand_extra_linewidths(data, src):
    """Preserve independent static or temperature-dependent disorder tables."""
    result = {}
    for name in EXTRA_LINEWIDTHS:
        if name not in data:
            continue
        values = np.asarray(data[name])
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f'{name} must contain finite nonnegative THz linewidths')
        if values.shape == data['frequency'].shape:
            result[name] = values[src, :]
        elif values.shape == data['gamma'].shape:
            result[name] = values[:, src, :]
        else:
            raise ValueError(f'{name} shape mismatch')
    return result


def write_rate_summary(destination: Path):
    """Heat-capacity-weighted inverse lifetimes, not inverse mean lifetimes."""
    with h5py.File(destination / 'kappa-fbz.hdf5') as f:
        frequency = f['frequency'][:]
        capacity = f['heat_capacity'][:]
        capacity[:, frequency < 1e-4] = 0
        total_capacity = capacity.sum(axis=(1, 2))
        normal = (capacity * f['rate_N'][:]).sum(axis=(1, 2)) / total_capacity
        umklapp = (capacity * f['rate_U'][:]).sum(axis=(1, 2)) / total_capacity
        temperatures = f['temperature'][:]
    with (destination / 'scattering_summary.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['temperature_K', 'Cv_weighted_N_rate_ps^-1', 'Cv_weighted_U_rate_ps^-1', 'N_over_U_weighted_rate'])
        writer.writerows(zip(temperatures, normal, umklapp, normal / umklapp, strict=True))


def export(job: Path, destination: Path):
    import phonopy
    job, destination = job.resolve(), destination.resolve()
    manifest = json.loads((job / 'manifest.json').read_text())
    if manifest['status'] != 'calculated':
        raise RuntimeError('The phono3py calculation has not completed.')
    if destination.exists():
        raise RuntimeError(f'Refusing to overwrite {destination}')
    with h5py.File(job / 'harmonic-grid.hdf5') as f:
        harmonic = {key: f[key][()] for key in f}
    mesh = harmonic['mesh']
    raw = job / ('kappa-m' + ''.join(map(str, mesh)) + '.hdf5')
    with h5py.File(raw) as f:
        data = {key: f[key][()] for key in f}
    q = harmonic['qpoint']
    nq = int(np.prod(mesh))
    np.testing.assert_array_equal(mesh, manifest['mesh'])
    np.testing.assert_array_equal(data['temperature'], manifest['temperatures_K'])
    np.testing.assert_array_equal(data['grid_point'], harmonic['grg2bzg'][harmonic['ir_grg']])
    np.testing.assert_array_equal(data['weight'], harmonic['ir_weight'])
    assert data['weight'].sum() == nq == len(q)
    assert np.unique(np.rint(q * mesh).astype(int) % mesh, axis=0).shape[0] == nq
    np.testing.assert_allclose(q * mesh, np.rint(q * mesh), atol=1e-9)
    lookup = {int(g): i for i, g in enumerate(harmonic['ir_grg'])}
    src = np.array([lookup[int(g)] for g in harmonic['ir_map']], dtype=np.int64)
    np.testing.assert_array_equal(np.bincount(src), data['weight'])
    freq = data['frequency'][src]
    freq_error = float(np.max(np.abs(freq - harmonic['frequency'])))
    np.testing.assert_allclose(freq, harmonic['frequency'], atol=5e-6, rtol=1e-7)
    for name in ('gamma', 'gamma_N', 'gamma_U', 'kappa', 'heat_capacity'):
        assert np.isfinite(data[name]).all(), name
    for name in ('gamma', 'gamma_N', 'gamma_U'):
        assert (data[name] >= 0).all(), name
    np.testing.assert_allclose(data['gamma_N'] + data['gamma_U'], data['gamma'], atol=1e-14, rtol=1e-7)
    sum_error = float(np.max(np.abs(data['gamma_N'] + data['gamma_U'] - data['gamma'])))
    extras = expand_extra_linewidths(data, src)

    ph = phonopy.load(job / 'input/phono3py.yaml',
                      force_constants_filename=job / 'input/fc2.hdf5',
                      produce_fc=False, is_nac=False, lang='C')
    np.testing.assert_allclose(ph.primitive.cell, harmonic['primitive_lattice'], atol=1e-12)
    ph.run_qpoints(q, with_group_velocities=True)
    phonons = ph.get_qpoints_dict()
    np.testing.assert_allclose(phonons['frequencies'], harmonic['frequency'], atol=5e-6, rtol=1e-7)
    gv = phonons['group_velocities']
    assert np.isfinite(gv).all()
    # phono3py's own kappa includes isotope scattering when provided. Custom
    # impurity/defect tables need a separate total-rate conductivity below.
    gamma = data['gamma'][:, src, :] + extras.get('gamma_isotope', 0)
    cv = data['heat_capacity'][:, src, :]
    inverse_2gamma = np.divide(1, 2 * gamma, out=np.zeros_like(gamma), where=gamma > 0)
    inverse_2gamma[:, freq < 1e-4] = 0
    pairs = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))
    vv = np.stack([gv[..., a] * gv[..., b] for a, b in pairs], axis=-1)
    recovered = np.einsum('tqb,qbc->tc', cv * inverse_2gamma, vv) * data['kappa_unit_conversion'] / nq
    kappa_error = float(np.max(np.abs(recovered - data['kappa'])) / np.max(np.abs(data['kappa'])))
    # This checks the complete expansion, units, and independently computed velocities.
    np.testing.assert_allclose(recovered, data['kappa'], rtol=1e-4, atol=1e-6)
    total_gamma = data['gamma'][:, src, :].copy()
    for values in extras.values():
        total_gamma += values
    if not np.isfinite(total_gamma).all():
        raise ValueError('Total linewidth overflow including resistive scattering')
    total_inverse = np.divide(1, 2 * total_gamma, out=np.zeros_like(total_gamma), where=total_gamma > 0)
    total_inverse[:, freq < 1e-4] = 0
    total_kappa = (np.einsum('tqb,qbc->tc', cv * total_inverse, vv)
                   * data['kappa_unit_conversion'] / nq) if extras else data['kappa']
    if (freq[np.linalg.norm(q, axis=1) > 1e-10] <= 0).any():
        raise RuntimeError('Non-Gamma imaginary/zero modes require review before transport export.')
    # PhonoMC selects positive frequencies, whereas phono3py selects >= cutoff.
    # Make the mode sets identical, retaining all original values separately.
    transport_frequency = freq.copy()
    transport_frequency[freq < 1e-4] = 0
    validation = {
        'full_grid_points': nq, 'irreducible_grid_points': len(data['weight']),
        'temperatures_K': data['temperature'].tolist(),
        'NU_sum_max_absolute_error_THz': sum_error,
        'symmetry_frequency_max_absolute_error_THz': freq_error,
        'full_grid_RTA_kappa_max_relative_error': kappa_error,
        'negative_frequency_modes_below_minus_0_01_THz': int((freq < -0.01).sum()),
        'minimum_frequency_THz': float(freq.min()),
        'transport_modes_excluded_by_phono3py_cutoff': int((freq < 1e-4).sum()),
        'minimum_non_Gamma_frequency_THz': float(freq[np.linalg.norm(q, axis=1) > 1e-10].min()),
        'Gamma_frequencies_THz': freq[np.linalg.norm(q, axis=1) < 1e-10].tolist(),
        'spacegroup_at_symprec_1e-5': ph.primitive_symmetry.dataset.international,
        'raw_kappa_file': str(raw),
        'extra_resistive_linewidths': list(extras),
    }
    destination.mkdir(parents=True)
    shutil.copy2(job / 'POSCAR.primitive', destination / 'POSCAR')
    with h5py.File(destination / 'kappa-fbz.hdf5', 'w') as f:
        f.attrs['source_job'] = str(job)
        f.attrs['method'] = 'RTA; three-phonon plus ' + ', '.join(extras) if extras else 'RTA; intrinsic three-phonon only'
        f.attrs['gamma_units'] = 'THz HWHM; rate_ps^-1 = 4*pi*gamma'
        f.attrs['group_velocity_units'] = 'THz Angstrom (= 100 m/s)'
        f.attrs['qpoint_convention'] = 'fractional primitive reciprocal coordinates; one GRG representative in first BZ'
        f.attrs['transport_frequency_cutoff_THz'] = 1e-4
        values = {'qpoint': q, 'mesh': mesh, 'weight': np.ones(nq, dtype=np.int64),
                  'frequency': transport_frequency, 'frequency_raw_phono3py': freq, 'group_velocity': gv,
                  'temperature': data['temperature'], 'heat_capacity': cv,
                  'kappa': total_kappa, 'kappa_phono3py': data['kappa'], '_fbz_src_q_index': src}
        values.update(extras)
        for name in ('gamma', 'gamma_N', 'gamma_U'):
            values[name] = data[name][:, src, :]
        for name in ('N', 'U'):
            values['rate_' + name] = 4 * np.pi * values['gamma_' + name]
        resistive_gamma = values['gamma_U'].copy()
        for extra in extras.values():
            resistive_gamma += extra
        values['rate_R'] = 4 * np.pi * resistive_gamma
        for key, value in values.items():
            f.create_dataset(key, data=value, compression='gzip')
        f['rate_N'].attrs['units'] = 'ps^-1'
        f['rate_U'].attrs['units'] = 'ps^-1'
        f['rate_R'].attrs['units'] = 'ps^-1'
    write_rate_summary(destination)
    sheet = mesh[2] == 1
    thickness_m = abs(np.linalg.det(ph.primitive.cell)) / np.linalg.norm(np.cross(ph.primitive.cell[0], ph.primitive.cell[1])) * 1e-10
    with (destination / 'bulk_kappa_RTA.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['temperature_K'] + [f'kappa_{s}_W_mK' for s in ('xx', 'yy', 'zz', 'yz', 'xz', 'xy')]
                        + (['sheet_xx_W_K', 'sheet_yy_W_K', 'kappa_xx_at_3.35A_W_mK', 'kappa_yy_at_3.35A_W_mK'] if sheet else []))
        for temp, kappa in zip(data['temperature'], total_kappa, strict=True):
            writer.writerow([temp, *kappa] + ([kappa[0] * thickness_m, kappa[1] * thickness_m,
                                              kappa[0] * thickness_m / 3.35e-10, kappa[1] * thickness_m / 3.35e-10] if sheet else []))
    (destination / 'provenance.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (destination / 'validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    (destination / 'README.md').write_text(
        '# phono3py material data\n\n'
        f"Source job: `{job}`. Mesh: {mesh.tolist()}. Temperatures (K): {data['temperature'].tolist()}.\n\n"
        '`POSCAR` is the primitive cell matching the q coordinates. `kappa-fbz.hdf5` contains a complete mesh, '
        'unit weights, frequencies, velocities, and `gamma`, `gamma_N`, `gamma_U`. '
        f"Additional independent resistive tables: {', '.join(extras) or 'none'}. "
        'Linewidths are HWHM in THz; inverse lifetimes are `4*pi*gamma` in ps^-1. '
        '`rate_N` and `rate_U` contain these converted inverse lifetimes; `rate_R` includes U plus the extra tables.\n\n'
        '`scattering_summary.csv` contains heat-capacity-weighted N/U inverse lifetimes '
        'over modes above the 1e-4 THz cutoff; these averages alone do not establish hydrodynamic transport.\n\n'
        '`bulk_kappa_RTA.csv` and `kappa` use three-phonon rates plus the listed extra tables. '
        '`kappa_phono3py` preserves the original phono3py conductivity, before custom impurity/defect additions. '
        'No additional boundary, electron or four-phonon rates are inferred. This is neither LBTE nor Callaway conductivity. '
        'No mesh-convergence study has been performed.\n\n'
        + (f'The slab height is {thickness_m / 1e-10:.8g} Angstrom. Raw conductivity uses that height. '
           'Sheet thermal conductance is kappa times this height; the additional 3.35 Angstrom conductivity '
           'is an explicitly chosen thickness convention, not graphite conductivity.\n\n' if sheet else '')
        + 'PhonoMC: set `io.material_folder` to this directory and `scattering.model="callaway"`; '
        f"choose a background temperature within [{data['temperature'].min():g}, {data['temperature'].max():g}] K. "
        '`frequency_raw_phono3py` preserves the original frequencies. In the transport `frequency` dataset, '
        'values below the phono3py 1e-4 THz cutoff are set to zero so PhonoMC uses the same mode set. '
        'This only excludes the three Gamma acoustic modes; it does not modify force constants.\n'
    )
    print(json.dumps(validation, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    export(args.job, args.destination)
