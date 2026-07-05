#!/usr/bin/env python3
"""
viz_tls_heatmap.py
Render a TLS heatmap PNG with filtered polygon outlines overlaid.

Usage:
    python viz_tls_heatmap.py SID [SID2 ...]  [--out-dir DIR] [--cc N]

Reads:
    outputs/tls_masks/{SID}/images/{SID}_asap_hooknettls_heat1.tif
    outputs/tls_masks/{SID}/images/filtered/{SID}_asap_hooknettls_tls_filtered.json
    outputs/tls_masks/{SID}/images/filtered/{SID}_asap_hooknettls_gc_filtered.json

Writes:
    {OUT_DIR}/{SID}_tls_viz.png
"""
import sys
import json
import argparse
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection
import tifffile

PROJECT = Path(__file__).resolve().parents[2]
OUT_BASE = PROJECT / 'outputs' / 'tls_masks'

def load_heatmap(sid):
    p = OUT_BASE / sid / 'images' / f'{sid}_asap_hooknettls_heat1.tif'
    if not p.exists():
        raise FileNotFoundError(p)
    with tifffile.TiffFile(str(p)) as tf:
        # Read lowest-resolution level for speed
        series = tf.series[0]
        levels = series.levels
        # Pick a level that's manageable (~2000 px on longest side)
        for lvl in levels:
            page = lvl.pages[0]
            h, w = page.shape[:2]
            if max(h, w) <= 4000:
                img = lvl.asarray()
                scale = page.shape[1] / levels[0].pages[0].shape[1]
                return img, scale
        # Fallback: full resolution
        img = levels[0].asarray()
        return img, 1.0

def load_polygons(sid, kind):
    p = OUT_BASE / sid / 'images' / 'filtered' / f'{sid}_asap_hooknettls_{kind}_filtered.json'
    if not p.exists():
        # try old layout
        p = OUT_BASE / sid / 'filtered' / f'{sid}_asap_hooknettls_{kind}_filtered.json'
    if not p.exists():
        return []
    raw = json.loads(p.read_text())
    feats = raw.get('features', raw) if isinstance(raw, dict) else raw
    polys = []
    for f in feats:
        coords = f.get('coordinates', [])
        if coords:
            polys.append(np.array(coords))
    return polys

def render(sid, out_dir):
    print(f'  {sid}', end='', flush=True)
    try:
        heatmap, scale = load_heatmap(sid)
    except FileNotFoundError as e:
        print(f'  SKIP (no heatmap): {e}')
        return

    tls_polys = load_polygons(sid, 'tls')
    gc_polys  = load_polygons(sid, 'gc')

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.imshow(heatmap, cmap='hot', vmin=0, vmax=255, aspect='equal')

    def add_polys(polys, color, label):
        patches = []
        for coords in polys:
            # Scale from full-res pixel coords to this pyramid level
            scaled = coords * scale
            patches.append(MplPolygon(scaled, closed=True))
        if patches:
            col = PatchCollection(patches, facecolor='none',
                                  edgecolor=color, linewidth=1.5, label=label)
            ax.add_collection(col)

    add_polys(tls_polys, '#00ff00', f'TLS (n={len(tls_polys)})')
    add_polys(gc_polys,  '#00aaff', f'GC  (n={len(gc_polys)})')

    tls_patch = mpatches.Patch(edgecolor='#00ff00', facecolor='none',
                                label=f'TLS (n={len(tls_polys)})', linewidth=1.5)
    gc_patch  = mpatches.Patch(edgecolor='#00aaff', facecolor='none',
                                label=f'GC  (n={len(gc_polys)})',  linewidth=1.5)
    ax.legend(handles=[tls_patch, gc_patch], loc='upper right', fontsize=9,
              framealpha=0.7)

    ax.set_title(f'{sid}  |  TLS={len(tls_polys)}  GC={len(gc_polys)}', fontsize=13)
    ax.axis('off')

    out_path = Path(out_dir) / f'{sid}_tls_viz.png'
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  -> {out_path}')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('sids', nargs='+')
    ap.add_argument('--out-dir', default=str(PROJECT / 'outputs' / 'viz'))
    args = ap.parse_args()

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    print(f'Writing PNGs to {args.out_dir}')
    for sid in args.sids:
        render(sid, args.out_dir)

if __name__ == '__main__':
    main()
