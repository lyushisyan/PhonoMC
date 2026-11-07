#!/usr/bin/env python3
"""Write selected Gaussian source fragments and previews, without running transport."""
import argparse
import json
import math
from pathlib import Path
import tomllib

from generate_matched_fets import preview


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--power-w',type=float,default=1e-6,
                        help='Reference total lattice-heating power per FET (W)')
    args = parser.parse_args()
    if not math.isfinite(args.power_w) or args.power_w <= 0:
        parser.error('Power must be finite and positive')
    out = Path(__file__).resolve().parents[1]/'model'/'matched_fets'
    geometry = json.loads((out/'geometry.json').read_text())
    devices = {name:d['boxes'] for name,d in geometry['devices'].items()}
    length = min(6.,geometry['parameters']['channel_length'])
    report = dict(selected_profile='gaussian',transverse_profile='uniform',
                  solver_axis_cutoff_sigma=2,transport_simulated=False,
                  power_is_reference_value=True,devices={})
    for name,parts in devices.items():
        channels = parts['channel']
        xmax = channels[0][1]
        xmin = xmax-length
        lo = [xmin,min(b[2] for b in channels),min(b[4] for b in channels)]
        hi = [xmax,max(b[3] for b in channels),max(b[5] for b in channels)]
        center = [(a+b)/2 for a,b in zip(lo,hi)]
        sigma = [length/4,0.,0.]
        body = f'''# Heat-source section only; merge into a complete simulation input.
# Total power is a reference value per FET, not per GAAFET layer.
# X is Gaussian with the solver's +/-2 sigma cutoff; Y/Z are uniform.
# Only active material cells are heated; GAAFET gaps are excluded.
[heat_source]
enabled = true
profile = "gaussian"
center = {json.dumps(center)}
sigma = {json.dumps(sigma)}
min = {json.dumps(lo)}
max = {json.dumps(hi)}
total_power = {args.power_w:.12g}
spectrum = "thermal"
time_profile = "constant"
start_time = 0
'''
        path = out/(name.lower()+'_heat_source.toml')
        path.write_text(body)
        source = tomllib.loads(path.read_text())['heat_source']
        assert math.isclose(source['total_power'],args.power_w,rel_tol=1e-11)
        # Check sampled power sharing on a 1 nm grid; this is source arithmetic,
        # not a local solver run. Current geometry has integer cell boundaries.
        per_layer = []
        for box in channels:
            weight = 0.
            for ix in range(math.floor(xmin),math.ceil(xmax)):
                x = ix+.5
                if xmin <= x <= xmax:
                    weight += math.exp(-.5*((x-center[0])/sigma[0])**2)
            weight *= (box[3]-box[2])*(box[5]-box[4])
            per_layer.append(weight)
        fractions = [w/sum(per_layer) for w in per_layer]
        assert all(math.isclose(f,1/len(channels)) for f in fractions)
        report['devices'][name] = dict(configuration=source,
                    sampled_layer_power_fractions=fractions,
                    spatial_weighted_volume_nm3=sum(per_layer))
    assert math.isclose(report['devices']['FinFET']['spatial_weighted_volume_nm3'],
                        report['devices']['GAAFET']['spatial_weighted_volume_nm3'])
    (out/'selected_heat_source.json').write_text(json.dumps(report,indent=2)+'\n')
    preview(out,devices,geometry['parameters'],heat_length=length,heat_profile='gaussian')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
