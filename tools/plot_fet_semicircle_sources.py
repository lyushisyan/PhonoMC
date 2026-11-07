#!/usr/bin/env python3
"""Illustrate proposed semicircular heat sources; no solver input is changed."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge, Patch


def main():
    out = Path(__file__).resolve().parents[1] / 'model' / 'matched_fets'
    data = json.loads((out/'geometry.json').read_text())
    colors = dict(substrate='#b5bdc9',source='#e5aa57',drain='#e5aa57',channel='#368ac1')
    red = '#d74737'
    fig, axes = plt.subplots(1,2,figsize=(12,6.5))
    fig.suptitle('Proposed top-localized semicircular heat sources',fontsize=18,y=.96)
    report = dict(description='Illustrative X-Z half disks extruded uniformly across channel width',
                  simulation_enabled=False,units='nm',devices={})
    for ax,(name,device) in zip(axes,data['devices'].items()):
        parts = device['boxes']
        sources = []
        for part,boxes in parts.items():
            for x0,x1,y0,y1,z0,z1 in boxes:
                ax.add_patch(Rectangle((x0,z0),x1-x0,z1-z0,facecolor=colors[part],edgecolor='#475569',lw=.8))
                if part in ('source','drain'):
                    ax.text((x0+x1)/2,(z0+z1)/2,part.title(),ha='center',va='center',rotation=90)
                if part == 'channel':
                    radius = min(3.,(x1-x0)/2,z1-z0)
                    cx = x1-radius
                    ax.add_patch(Wedge((cx,z1),radius,180,360,facecolor=red,edgecolor='#9c2828',lw=.8))
                    label = 'P' if name=='FinFET' else f'P/{len(boxes)}'
                    ax.annotate(label,xy=(cx,z1-radius/2),xytext=(cx-11,z1+4),color='#9c2828',
                                fontsize=11,ha='center',arrowprops=dict(arrowstyle='->',color='#9c2828',lw=1.1))
                    sources.append(dict(center_xz_nm=[cx,z1],radius_nm=radius,y_range_nm=[y0,y1],
                                        z_half='below_center',volume_nm3=.5*3.141592653589793*radius**2*(y1-y0),
                                        fraction_of_device_power=1/len(boxes)))
        report['devices'][name]=sources
        ax.text(26,10,'Substrate',ha='center',color='#374151')
        ax.set(xlim=(0,52),ylim=(0,57),xlabel='X (nm)',ylabel='Z (nm)',title=name+' — slice at Y = 12 nm')
        ax.set_aspect('equal')
        ax.set_xticks([0,10,36,42,52])
        ax.spines[['top','right']].set_visible(False)
    fig.legend(handles=[Patch(color=colors['channel'],label='Channel'),Patch(color=red,label='Assumed volumetric heat source')],
               loc='lower center',bbox_to_anchor=(.5,.075),ncol=2,frameon=False)
    fig.text(.5,.055,'Radius = 3 nm; diameter spans X = 36–42 nm; uniform across each channel width.',ha='center',fontsize=10)
    fig.text(.5,.02,'Illustrative source shape, not an electron-transport prediction. Same total power P per FET.',ha='center',fontsize=10,color='#64748b')
    fig.tight_layout(rect=(0,.13,1,.93))
    fig.savefig(out/'heat_sources_semicircle.png',dpi=200)
    fig.savefig(out/'heat_sources_semicircle.svg')
    (out/'heat_sources_semicircle.json').write_text(json.dumps(report,indent=2)+'\n')
    plt.close(fig)


if __name__ == '__main__':
    main()
