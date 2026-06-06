#!/usr/bin/env python3
"""
initRAReference.py

Build the vr/va reference basemap for one velocityStats frame-range directory
(e.g. 0040-0049).  The basemap is placed in <frameRangeDir>/sims/ and is
consumed by velocityStats --doRA to initialise the outlier-rejection reference.

Algorithm
---------
1. Scan the track directory for all *_00n? frame subdirs that fall inside the
   frame range and collect their primary geodat geojson files.
2. Build a "super-geodat" that spans the combined azimuth extent of all found
   frames (plus ±10 km padding) and extends the range swath by ±5 km.
3. Clone the super-geodat and shift all times by deltaT days to create a
   zero-baseline secondary geodat.
4. Run simoffsets to produce synthetic range/azimuth offset files.
5. Write dummy (all-zero) az.est, rBaseline, and baseline.26x16 param files.
6. Write a mosaic3d inputFile and run mosaic3d -rOffsets -outputRA to produce
   sims/mosaicOffsets.vr and sims/mosaicOffsets.va.

Run from velocityStats/00n0-00n9/ or pass explicit paths.
"""

import argparse
import copy
import glob
import json
import math
import os
from datetime import datetime, timedelta
from subprocess import call

import sarfunc as s


# ---- helpers ----------------------------------------------------------------

def _parseNominalTime(timeStr):
    """Parse 'HH:MM:SS.ffffff' → datetime (arbitrary date 2000-01-01)."""
    base = '2000-01-01 ' + timeStr.strip()
    for fmt in ('%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S'):
        try:
            return datetime.strptime(base, fmt)
        except ValueError:
            pass
    raise ValueError(f'Cannot parse NominalTime: {timeStr!r}')


def _formatNominalTime(dt):
    """Format datetime back to 'HH:MM:SS.ffffff'."""
    return dt.strftime('%H:%M:%S.%f')


def _parseCorrectedTime(timeStr):
    """Parse 'HH MM SS.SSSSSSS' → datetime."""
    parts = timeStr.strip().split()
    hh, mm = int(parts[0]), int(parts[1])
    sec = float(parts[2])
    ss = int(sec)
    us = round((sec - ss) * 1e6)
    return datetime(2000, 1, 1, hh, mm, ss, us)


def _formatCorrectedTime(dt):
    """Format datetime back to 'HH MM SS.SSSSSSS'."""
    frac = dt.microsecond / 1e6
    return f'{dt.hour:2d} {dt.minute:2d} {dt.second + frac:.7f}'


def _secOfDay(dt):
    """Seconds since midnight for a datetime."""
    return dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond / 1e6


# ---- geodat building --------------------------------------------------------

def buildSuperGeodat(geodatFiles, simsDir):
    """
    Create a super-geodat from a list of primary geodat geojson files.

    Picks the frame with the largest MLAzimuthSize as the structural base,
    then extends its azimuth coverage to span all frames ±10 km and its range
    swath by ±5 km.  State vectors from all frames are merged and sorted.

    Returns the super-geodat dict (GeoJSON Feature properties).
    Writes to <simsDir>/super.geojson.
    """
    geojsons = []
    for f in geodatFiles:
        with open(f) as fp:
            geojsons.append(json.load(fp))

    props = [g['properties'] for g in geojsons]

    # Pick base: largest azimuth coverage
    base = max(props, key=lambda p: p['MLAzimuthSize'])
    result = copy.deepcopy(geojsons[max(range(len(props)),
                                        key=lambda i: props[i]['MLAzimuthSize'])])
    p = result['properties']

    nlr = p['NumberRangeLooks']
    nla = p['NumberAzimuthLooks']
    slpR = p['SLCRangePixelSize']   # m
    slpA = p['SLCAzimuthPixelSize']  # m
    prf  = p['PRF']                  # Hz
    mlPixR = slpR * nlr              # m per ML range pixel
    mlPixA = slpA * nla              # m per ML azimuth pixel
    pixPeriod = nla / prf            # seconds per ML azimuth line

    # ---- combined azimuth time span ----
    tStarts = []
    tEnds   = []
    for pp in props:
        tStart = _parseNominalTime(pp['NominalTime'])
        tEnd   = tStart + timedelta(seconds=pp['MLAzimuthSize'] * pixPeriod)
        tStarts.append(tStart)
        tEnds.append(tEnd)
    tSpanStart = min(tStarts)
    tSpanEnd   = max(tEnds)

    # ---- padding ----
    padAzSec = 10000.0 / mlPixA * pixPeriod   # seconds for 10 km
    padRangePx = round(5000.0 / mlPixR)        # ML pixels for 5 km

    newStart = tSpanStart - timedelta(seconds=padAzSec)
    newEnd   = tSpanEnd   + timedelta(seconds=padAzSec)
    newAzSize = round((newEnd - newStart).total_seconds() / pixPeriod)

    newNearRange = p['MLNearRange'] - 5000.0
    newFarRange  = p['MLFarRange']  + 5000.0
    newRangeSize = p['MLRangeSize'] + 2 * padRangePx
    newCenterRange = (newNearRange + newFarRange) / 2.0

    p['MLAzimuthSize'] = newAzSize
    p['MLRangeSize']   = newRangeSize
    p['MLNearRange']   = newNearRange
    p['MLFarRange']    = newFarRange
    p['MLCenterRange'] = newCenterRange
    p['NominalTime']   = _formatNominalTime(newStart)
    p['CorrectedTime'] = _formatCorrectedTime(newStart)

    # ---- merge state vectors from all frames ----
    svInterval = p.get('StateVectorInterval', 10.0)
    svAll = {}  # keyed by time (s-of-day) to avoid duplicates
    for pp in props:
        nSV = pp.get('NumberOfStateVectors', 0)
        t0  = pp.get('TimeOfFirstStateVector', 0.0)
        for i in range(1, nSV + 1):
            t = t0 + (i - 1) * svInterval
            if t not in svAll:
                svAll[t] = {
                    'pos': pp.get(f'SV_Pos_{i}'),
                    'vel': pp.get(f'SV_Vel_{i}'),
                }
    sortedTimes = sorted(svAll)
    # remove old SV keys from result
    oldN = p.get('NumberOfStateVectors', 0)
    for i in range(1, oldN + 1):
        p.pop(f'SV_Pos_{i}', None)
        p.pop(f'SV_Vel_{i}', None)
    # write merged SVs
    p['NumberOfStateVectors']  = len(sortedTimes)
    p['TimeOfFirstStateVector'] = sortedTimes[0] if sortedTimes else 0.0
    for idx, t in enumerate(sortedTimes, start=1):
        p[f'SV_Pos_{idx}'] = svAll[t]['pos']
        p[f'SV_Vel_{idx}'] = svAll[t]['vel']

    outPath = os.path.join(simsDir, 'super.geojson')
    with open(outPath, 'w') as fp:
        json.dump(result, fp, indent=2)
    print(f'Super-geodat written: {outPath}  '
          f'(az={newAzSize} rg={newRangeSize})')
    return result


def buildSecondaryGeodat(superGeojson, simsDir, deltaT=12):
    """
    Clone the super-geodat and shift all times by deltaT days.

    Because primary and secondary share the same orbit geometry the
    perpendicular baseline is zero, so all offset signal is velocity.
    State vector positions are shifted using linear extrapolation from the
    primary velocity vectors (sufficient for the short deltaT).
    """
    result = copy.deepcopy(superGeojson)
    p = result['properties']
    dt_sec = deltaT * 86400.0

    # Shift NominalTime and CorrectedTime
    tNom = _parseNominalTime(p['NominalTime'])
    tNom += timedelta(seconds=dt_sec)
    p['NominalTime']   = _formatNominalTime(tNom)
    p['CorrectedTime'] = _formatCorrectedTime(tNom)

    # Shift TimeOfFirstStateVector (mod 86400 to stay within day)
    p['TimeOfFirstStateVector'] = (p['TimeOfFirstStateVector'] + dt_sec) % 86400.0

    # Propagate SV positions using linear v·Δt
    nSV = p.get('NumberOfStateVectors', 0)
    for i in range(1, nSV + 1):
        pos = p[f'SV_Pos_{i}']
        vel = p[f'SV_Vel_{i}']
        if pos is not None and vel is not None:
            p[f'SV_Pos_{i}'] = [pos[j] + vel[j] * dt_sec for j in range(3)]

    outPath = os.path.join(simsDir, 'super.secondary.geojson')
    with open(outPath, 'w') as fp:
        json.dump(result, fp, indent=2)
    print(f'Secondary geodat written: {outPath}  (deltaT={deltaT} days)')
    return result


# ---- dummy param files ------------------------------------------------------

def createDummyParams(simsDir, superGeojson):
    """
    Write zero-polynomial az.est, rBaseline, and baseline.26x16 param files.

    Because the secondary geodat has zero baseline (same orbit, shifted time)
    all polynomial corrections are exactly zero.
    """
    p = superGeojson['properties']
    H   = p['SpaceCraftAltitude'] / 1000.0    # m → km
    Re  = (p['EarthRadiusMajor'] + p['EarthRadiusMinor']) / 2000.0  # mean, km
    RNear  = p['MLNearRange']   / 1000.0      # m → km
    Rc     = p['MLCenterRange'] / 1000.0
    nlr    = p['NumberRangeLooks']
    nla    = p['NumberAzimuthLooks']
    slpR   = p['SLCRangePixelSize']
    wl     = p['Wavelength']

    header = (
        f';;\n'
        f';; H, InSAR parameters  Re, RNear, Rc, nlr,nla\n'
        f';;\n'
        f'    {H:.3f}   {Re:.3f}   {RNear:.3f}  {Rc:.3f}'
        f'  {nlr} {nla} {slpR:.6f} {wl:.6f}\n'
        f';; Tide difference : none\n'
    )

    for fname in ('az.est', 'rBaseline'):
        with open(os.path.join(simsDir, fname), 'w') as fp:
            fp.write(header)
            fp.write(';\n;&\n;\n')
            fp.write('0.0 0.0 0.0 0.0\n')

    # baseline.26x16: two lines of zero baseline data
    baselineHdr = (
        f';;\n'
        f';; InSAR parameters H, Re, RNear, Rc, nlr,nla\n'
        f';;\n'
        f' {H:.3f}   {Re:.3f}  {RNear:.3f} {Rc:.3f}  {nlr} {nla} {slpR:.6f}\n'
        f';; Tide difference : none\n'
        f';\n'
        f'; Number of lines of baseline data\n'
        f';\n'
        f' 2\n'
        f';\n'
        f'; First flattening baseline\n'
        f';\n'
        f'0.0 0.0 0.0 0.0\n'
        f';\n'
        f'; Second flattening baseline\n'
        f';\n'
        f'0.0 0.0 0.0 0.0\n'
    )
    with open(os.path.join(simsDir, 'baseline.26x16'), 'w') as fp:
        fp.write(baselineHdr)
    print(f'Dummy param files written in {simsDir}')


# ---- simoffsets + mosaic3d --------------------------------------------------

def runSimoffsets(simsDir, region, dem, ompThreads=4):
    """Run simoffsets inside simsDir to produce synthetic offset files."""
    cmd = (
        f'cd {simsDir}; '
        f'simoffsets -region {region} -syncDat '
        f'-geodatFile super.geojson '
        f'-secondGeodatFile super.secondary.geojson '
        f'-offsetsDat offsets.dat '
        f'-azOffsets offsets.da '
        f'-dem {dem} '
        f'--ompThreads {ompThreads}'
    )
    print(cmd)
    call(cmd, shell=True)


def createInputFile(simsDir, deltaT=12):
    """Write the mosaic3d inputFile inside simsDir."""
    absDir = os.path.abspath(simsDir)
    content = (
        f'0 0 0 0 0.2 0.2\n'
        f';\n'
        f'1\n'
        f';\n'
        f'nophase {absDir}/super.geojson {absDir}/baseline.26x16 '
        f'{deltaT} 1. '
        f'{absDir}/offsets.da {absDir}/az.est '
        f'{absDir}/offsets.dr {absDir}/rBaseline\n'
    )
    path = os.path.join(simsDir, 'inputFile')
    with open(path, 'w') as fp:
        fp.write(content)
    print(f'mosaic3d inputFile written: {path}')


def runMosaic3d(simsDir, dem):
    """Run mosaic3d in RA speckle-track mode to produce mosaicOffsets.vr/va."""
    cmd = (
        f'cd {simsDir}; '
        f'mosaic3d -noVh -no3d -rOffsets -SVConst -outputRA -center '
        f'inputFile {dem} mosaicOffsets'
    )
    print(cmd)
    call(cmd, shell=True)


# ---- top-level entry point --------------------------------------------------

def initRAReference(frameRangeDir, trackDir, region, dem,
                    deltaT=12, ompThreads=4):
    """
    Build the vr/va reference basemap for one frame-range directory.

    Parameters
    ----------
    frameRangeDir : str  e.g. '0040-0049' or absolute path
    trackDir      : str  directory containing *_00?? frame subdirs
    region        : str  'greenland' or 'antarctica'
    dem           : str  path to DEM file
    deltaT        : int  repeat-pass interval in days (default 12)
    ompThreads    : int  OpenMP threads for simoffsets/siminsar
    """
    # resolve frame range
    dirName = os.path.basename(os.path.abspath(frameRangeDir))
    frame1, frame2 = (int(x) for x in dirName.split('-'))
    n = frame1 // 10    # tens digit group
    pattern = os.path.join(trackDir, f'*_{n:02d}??')
    frameDirs = sorted(glob.glob(pattern))
    if not frameDirs:
        raise RuntimeError(f'No frame dirs matching {pattern}')

    # collect geodat files for frames in range
    geodatFiles = []
    for fd in frameDirs:
        frameNum = int(fd.split('_')[-1])
        if frame1 <= frameNum <= frame2:
            candidates = glob.glob(os.path.join(fd, 'geodat*.geojson'))
            primary = [c for c in candidates if 'secondary' not in c]
            if primary:
                geodatFiles.append(primary[0])
    if not geodatFiles:
        raise RuntimeError(f'No primary geodat geojson files found for frames '
                           f'{frame1}-{frame2} in {trackDir}')
    print(f'Found {len(geodatFiles)} frame geodats for range {frame1}-{frame2}')

    # create sims subdir
    simsDir = os.path.join(os.path.abspath(frameRangeDir), 'sims')
    os.makedirs(simsDir, exist_ok=True)

    superGeo   = buildSuperGeodat(geodatFiles, simsDir)
    buildSecondaryGeodat(superGeo, simsDir, deltaT=deltaT)
    runSimoffsets(simsDir, region, dem, ompThreads=ompThreads)
    createDummyParams(simsDir, superGeo)
    createInputFile(simsDir, deltaT=deltaT)
    runMosaic3d(simsDir, dem)

    vrPath = os.path.join(simsDir, 'mosaicOffsets.vr')
    if os.path.exists(vrPath):
        print(f'Reference basemap ready: {vrPath}')
    else:
        print(f'WARNING: expected output not found: {vrPath}')


def main():
    parser = argparse.ArgumentParser(
        description='Build the vr/va reference basemap for a velocityStats '
                    'frame-range directory.',
        epilog='Part of the mosaicworkflow package.')
    parser.add_argument('--frameRangeDir', default='.',
                        help='Frame-range directory e.g. 0040-0049 [.]')
    parser.add_argument('--trackDir', default='../..',
                        help='Track directory containing *_00?? frame subdirs '
                             '[../.. from velocityStats/00n0-00n9/]')
    parser.add_argument('--region', default=None,
                        help='Region: greenland or antarctica [auto-detect]')
    parser.add_argument('--dem', default=None,
                        help='DEM file [region default]')
    parser.add_argument('--deltaT', type=int, default=12,
                        help='Repeat-pass interval in days [12]')
    parser.add_argument('--ompThreads', type=int, default=4,
                        help='OpenMP threads for simoffsets [4]')
    parser.add_argument('--initializeReference', action='store_true',
                        help='Rebuild even if basemap already exists')
    args = parser.parse_args()

    vrPath = os.path.join(os.path.abspath(args.frameRangeDir),
                          'sims', 'mosaicOffsets.vr')
    if os.path.exists(vrPath) and not args.initializeReference:
        print(f'Basemap already exists: {vrPath}  (use --initializeReference to rebuild)')
        return

    # resolve region and dem from project.yaml / sarfunc defaults if not given
    region = args.region
    dem    = args.dem
    if region is None or dem is None:
        # walk up to find project.yaml
        yamlPath = None
        cur = os.path.abspath(args.frameRangeDir)
        for _ in range(5):
            candidate = os.path.join(cur, 'project.yaml')
            if os.path.exists(candidate):
                yamlPath = candidate
                break
            cur = os.path.dirname(cur)
        if yamlPath:
            import yaml
            with open(yamlPath) as fp:
                proj = yaml.safe_load(fp) or {}
            if region is None:
                region = proj.get('region', 'greenland')
                # strip path if it's a full path to a region file
                if os.sep in str(region):
                    region = 'greenland'
        if region is None:
            region = 'greenland'
        if dem is None:
            myRegion = s.defaultRegionDefs(region)
            dem = myRegion.dem()

    initRAReference(args.frameRangeDir, args.trackDir, region, dem,
                    deltaT=args.deltaT, ompThreads=args.ompThreads)


if __name__ == '__main__':
    main()
