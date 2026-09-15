#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
checkvelsize - check that per-frame velocity mosaics share a common grid.

Replaces eyeballing `ls track-*/*_*/velocity/*.vx` for equal byte sizes, which
stopped working once the products became GeoTIFFs (compressed, so equal-sized
grids give unequal file sizes). Reads the grid from the sidecar that each
format already carries -- mosaicOffsets.vrt for tiff output, mosaicOffsets.*
.geodat for binary -- so it needs no marker files and works on products that
already exist, including a project part-way through the tiff migration.

setupS1Tracks --runVelstatsregions finds a common bounding box per
velocityStats range, rebuilds velocity/ on it and resizes velocity_nocull/ to
match, so frames are grouped by velocityStats range and velocity/ is compared
against velocity_nocull/ within each frame. Exits non-zero on any disagreement.

Run from the project root (holding track-*/) or inside a single track dir.

@author: ian
"""
import argparse
import glob
import os
import re
import sys
from collections import Counter, defaultdict

import utilities as u
from utilities import geodat

# Grids are compared to the nearest metre; 0.001 km absorbs the float noise of
# the km/metre round trip between the .geodat and GeoTransform conventions.
GRIDTOL = 0.001

VELDIRS = ['velocity', 'velocity_nocull']


def checkvelsizeArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\033[1mCheck per-frame velocity mosaics '
                    'share a grid\033[0m',
        epilog='Part of the mosaicworkflow package.')
    parser.add_argument('--tracks', nargs='+', metavar='track-N', default=None,
                        help='restrict to these tracks [all track-* dirs]')
    parser.add_argument('--velDirs', nargs='+', metavar='DIR', default=VELDIRS,
                        help=f'per-frame product dirs to check {VELDIRS}')
    parser.add_argument('--long', action='store_true',
                        help='list every frame, not just the mismatches')
    parser.add_argument('--missingIsError', action='store_true',
                        help='also exit non-zero when a frame has no product '
                        'in one of the velDirs [report only]')
    args = parser.parse_args()
    return args.tracks, args.velDirs, args.long, args.missingIsError


def getTrackDirs(tracks):
    '''
    Track dirs under the cwd, numeric-sorted. Running inside a single track dir
    (no track-* children) checks just that track.
    '''
    if tracks:
        trackDirs = [os.path.abspath(t) for t in tracks]
        missing = [d for d in trackDirs if not os.path.isdir(d)]
        if missing:
            u.myerror(f'Track directories not found: {missing}')
    else:
        trackDirs = glob.glob('track-*')
        if not trackDirs and re.search(r'track-\d+$', os.getcwd()):
            trackDirs = [os.getcwd()]
    return sorted(trackDirs,
                  key=lambda p: int(re.search(r'track-(\d+)',
                                              os.path.basename(p)).group(1)))


def getFrameRanges(trackDir):
    '''
    Burst ranges from velocityStats/<first>-<last>, which is what
    makevelstatsregions computes a common bounding box over. Returns a list of
    (name, first, last), or [] when the track has no velocityStats ranges.
    '''
    ranges = []
    statsDirs = glob.glob(os.path.join(trackDir, 'velocityStats', '*-*'))
    for vsd in sorted(statsDirs):
        name = os.path.basename(vsd)
        m = re.match(r'^(\d+)-(\d+)$', name)
        if m:
            ranges.append((name, int(m.group(1)), int(m.group(2))))
    return ranges


def frameRangeName(frameDir, ranges):
    ''' velocityStats range holding this <orbit>_<burst> dir, or
    unassigned '''
    m = re.match(r'^\d+_(\d+)$', os.path.basename(frameDir))
    if m is None:
        return None
    burst = int(m.group(1))
    for name, first, last in ranges:
        if first <= burst <= last:
            return name
    return 'unassigned'


def readGrid(prodDir):
    '''
    Grid of the velocity mosaic in prodDir, as (xs, ys, dx, dy, x0, y0) with
    x0/y0 the lower-left pixel centre in km. Prefers the tiff-mode VRT, falls
    back to the binary-mode .geodat sidecar, returns None if neither is there.

    Both are read through utilities.geodat so the two conventions (VRT
    GeoTransform is an upper-left pixel corner in metres, .geodat is a
    lower-left pixel centre in km) end up directly comparable -- that
    conversion lives in geodat.readGeodatFromTiff, not here.
    '''
    g = geodat(verbose=False)
    vrt = os.path.join(prodDir, 'mosaicOffsets.vrt')
    if os.path.exists(vrt):
        g.readGeodatFromTiff(vrt)
    else:
        sidecars = sorted(glob.glob(os.path.join(prodDir,
                                                 'mosaicOffsets.*.geodat')))
        if not sidecars:
            return None
        g.readGeodat(sidecars[0])
    return (g.xs, g.ys, g.dx, g.dy, g.x0, g.y0)


def gridKey(grid):
    ''' Hashable grid, rounded so float noise doesn't split equal grids '''
    xs, ys, dx, dy, x0, y0 = grid
    return (xs, ys, round(dx, 3), round(dy, 3),
            round(x0 / GRIDTOL), round(y0 / GRIDTOL))


def gridString(grid):
    xs, ys, dx, dy, x0, y0 = grid
    return f'{xs:>5} x {ys:<5} {dx:g}m  ({x0:9.1f},{y0:9.1f})'


def collectFrames(trackDir, velDirs):
    '''
    Every <orbit>_<burst> dir in the track, as a list of
    (frame, rangeName, {velDir: grid or None}). Frames with no product in any
    velDir are left out entirely -- they were never mosaicked, which is not
    what this tool is looking for.
    '''
    ranges = getFrameRanges(trackDir)
    frames = []
    for frameDir in sorted(glob.glob(os.path.join(trackDir, '*_*'))):
        if not os.path.isdir(frameDir):
            continue
        rangeName = frameRangeName(frameDir, ranges)
        if rangeName is None:
            continue
        grids = {}
        for velDir in velDirs:
            prodDir = os.path.join(frameDir, velDir)
            grids[velDir] = (readGrid(prodDir)
                             if os.path.isdir(prodDir) else None)
        if any(g is not None for g in grids.values()):
            frames.append((os.path.basename(frameDir), rangeName, grids))
    return frames


def regionsMode():
    ''' project.yaml velocityStatsRegions ('track' or 'frameBins') from the
    cwd (project root) or its parent (inside a track dir). '''
    for cand in ('project.yaml', '../project.yaml'):
        if os.path.exists(cand):
            try:
                import yaml
                with open(cand) as fp:
                    proj = yaml.safe_load(fp) or {}
                return 'track' if str(proj.get('velocityStatsRegions', 'frameBins')).lower() == 'track' else 'frameBins'
            except Exception:
                return 'frameBins'
    return 'frameBins'


def reportRangeTrackMode(track, rangeName, frames, velDirs, showLong):
    '''
    Track (master-grid) mode: products are meant to differ in size, so the test
    is alignment, not equality -- every product must share the most common
    posting and have its lower-left pixel-centre origin on an integer km.
    Returns (nMismatch, nMissing).
    '''
    postings = Counter((round(g[2], 3), round(g[3], 3)) for _, grids in frames
                       for g in grids.values() if g is not None)
    common = postings.most_common(1)[0][0]
    nMismatch = nMissing = 0
    lines = []
    for frame, grids in frames:
        for velDir in velDirs:
            grid = grids[velDir]
            if grid is None:
                nMissing += 1
                lines.append((True, f'    {frame}  {velDir:<15} (no product)'))
                continue
            xs, ys, dx, dy, x0, y0 = grid
            offGrid = (abs(x0 - round(x0)) > GRIDTOL or abs(y0 - round(y0)) > GRIDTOL
                       or (round(dx, 3), round(dy, 3)) != common)
            if offGrid:
                nMismatch += 1
                lines.append((True, f'    {frame}  {velDir:<15} {gridString(grid)}  <- off master grid'))
            else:
                lines.append((False, f'    {frame}  {velDir:<15} {gridString(grid)}'))
    verdict = 'OK' if nMismatch == 0 else f'OFF-GRID ({nMismatch})'
    print(f'{track:<10} {rangeName:<12} {len(frames):>4} frames  posting {common[0]:g}m  {verdict}')
    for flagged, line in lines:
        if flagged or showLong:
            print(line)
    return nMismatch, nMissing


def reportRange(track, rangeName, frames, velDirs, showLong):
    if regionsMode() == 'track':
        return reportRangeTrackMode(track, rangeName, frames, velDirs, showLong)
    '''
    Report one velocityStats range. The grid shared by the most products is
    taken as the expected one; anything else is a mismatch. Returns
    (nMismatch, nMissing).
    '''
    counts = Counter(gridKey(g) for _, grids in frames
                     for g in grids.values() if g is not None)
    common = counts.most_common(1)[0][0]
    example = next(g for _, grids in frames for g in grids.values()
                   if g is not None and gridKey(g) == common)

    nMismatch = 0
    nMissing = 0
    lines = []
    for frame, grids in frames:
        for velDir in velDirs:
            grid = grids[velDir]
            if grid is None:
                nMissing += 1
                lines.append((True, f'    {frame}  {velDir:<15} (no product)'))
            elif gridKey(grid) != common:
                nMismatch += 1
                lines.append((True, f'    {frame}  {velDir:<15} '
                                    f'{gridString(grid)}  <- differs'))
            else:
                lines.append((False, f'    {frame}  {velDir:<15} '
                                     f'{gridString(grid)}'))

    verdict = 'OK' if nMismatch == 0 else f'MISMATCH ({nMismatch})'
    print(f'{track:<10} {rangeName:<12} {len(frames):>4} frames  '
          f'{gridString(example)}  {verdict}')
    for flagged, line in lines:
        if flagged or showLong:
            print(line)
    return nMismatch, nMissing


def main():
    tracks, velDirs, showLong, missingIsError = checkvelsizeArgs()
    trackDirs = getTrackDirs(tracks)
    if not trackDirs:
        u.myerror('No track-* directories found -- run from the project root '
                  'or a track dir')

    totalMismatch, totalMissing, totalFrames, nRanges = 0, 0, 0, 0
    for trackDir in trackDirs:
        track = os.path.basename(trackDir)
        byRange = defaultdict(list)
        for frame, rangeName, grids in collectFrames(trackDir, velDirs):
            byRange[rangeName].append((frame, grids))
        if not byRange:
            print(f'{track:<10} no velocity products found')
            continue
        for rangeName in sorted(byRange):
            nMismatch, nMissing = reportRange(track, rangeName,
                                              byRange[rangeName], velDirs,
                                              showLong)
            totalMismatch += nMismatch
            totalMissing += nMissing
            totalFrames += len(byRange[rangeName])
            nRanges += 1

    print(f'\n{len(trackDirs)} track(s), {nRanges} range(s), '
          f'{totalFrames} frame(s): {totalMismatch} grid mismatch(es), '
          f'{totalMissing} missing product(s)')
    if totalMismatch or (missingIsError and totalMissing):
        sys.exit(1)


if __name__ == '__main__':
    main()
