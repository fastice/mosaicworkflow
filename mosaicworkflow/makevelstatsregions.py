#!/usr/bin/env python3
"""
makevelstatsregions.py

Run from a track-N/ directory.  For each velocityStats/F1-F2 frame-range
sub-directory, finds all velocity thumbnail outputs whose frame number falls
within [F1, F2], computes their combined geographic bounding box (in km), and
writes or updates a tiepoints/vel_thumb_header_F1dashF2 file with the resolved
resolution line.

Supports both binary (mosaicOffsets.vx.geodat) and tiff (mosaicOffsets.vrt)
output from mosaic3d.  The VRT is tried first; the geodat is the fallback.

If a header file already exists its existing lines are preserved and only the
resolution line is replaced.
"""

import argparse
import os
import yaml
import utilities as u
from osgeo import gdal

gdal.UseExceptions()


def _boundsFromVrt(vrtPath):
    """
    Return (xMin, yMin, xMax, yMax, dx, dy) in km from a GDAL-readable VRT/tiff.
    Coordinates are pixel-centred polar stereographic metres converted to km,
    matching the convention used by geodat boundsInKm().
    """
    ds = gdal.Open(vrtPath)
    gt = ds.GetGeoTransform()   # (xUL, dx, 0, yUL, 0, -dy)  all in metres
    ncols = ds.RasterXSize
    nrows = ds.RasterYSize
    ds = None
    dx = gt[1] / 1000.0
    dy = abs(gt[5]) / 1000.0
    # shift pixel-edge origin to pixel-centre by half a pixel on each side
    xMin = (gt[0] + 0.5 * gt[1]) / 1000.0
    xMax = (gt[0] + (ncols - 0.5) * gt[1]) / 1000.0
    yMax = (gt[3] + 0.5 * gt[5]) / 1000.0       # gt[5] < 0, moves inward
    yMin = (gt[3] + (nrows - 0.5) * gt[5]) / 1000.0
    return xMin, yMin, xMax, yMax, dx, dy


def _boundsFromVelDir(velDir):
    """
    Return (xMin, yMin, xMax, yMax, dx, dy) in km for a velocity directory.
    Tries mosaicOffsets.vrt first, falls back to mosaicOffsets.vx.geodat.
    Returns None if neither is found.
    """
    vrtPath = os.path.join(velDir, 'mosaicOffsets.vrt')
    if os.path.exists(vrtPath):
        x1, y1, x2, y2, dx, dy = _boundsFromVrt(vrtPath)
        return x1, y1, x2, y2, dx, dy

    geodatPath = os.path.join(velDir, 'mosaicOffsets.vx.geodat')
    if os.path.exists(geodatPath):
        geo = u.geodat()
        geo.readGeodat(geodatPath)
        x1, y1, x2, y2 = geo.boundsInKm()
        dx, dy = geo.pixSizeInKm()
        return x1, y1, x2, y2, dx, dy

    return None


def writeVelThumb(frameDir, x0, y0, dx, dy, xs, ys, baseDir, dem=None):
    """Write or update a vel_thumb_header file for the given frame-range dir."""
    d = os.getcwd()
    frameRange = frameDir.replace('-', 'dash')
    filename = os.path.join(baseDir, 'tiepoints', f'vel_thumb_header_{frameRange}')

    existingLines = None
    if os.path.exists(filename):
        with open(filename) as fp:
            existingLines = fp.readlines()

    with open(filename, 'w') as fT:
        if existingLines is None:
            if dem is None:
                dem = '/Volumes/insar7/ian/gimp/gimp2/270m/dem.gimp2.270m'
            u.pushd(baseDir)
            pwd = os.getcwd()
            u.popd()
            print(f'DEM       = {dem}', file=fT)
            print(f'track_root= {pwd}', file=fT)
            print(f'tie_dir   = {pwd}/tiepoints', file=fT)
            print(f'vel_dir = {pwd}/velocity', file=fT)
            print('extra_flags ="-no3d"', file=fT)
        else:
            for line in existingLines:
                if 'resolution' not in line:
                    print(line, file=fT, end='')
        print(f'resolution ="', x0, y0, xs, ys, dx, dy, '"', file=fT)


def main():
    parser = argparse.ArgumentParser(
        description='Find common geographic regions for velocityStats and write '
                    'vel_thumb_header files.  Run from the track-N/ directory.',
        epilog='Part of the mosaicworkflow package.')
    parser.add_argument('-dem', '--dem', default=None,
                        help='DEM path to embed in new header files '
                             '[auto-detected from project.yaml]')
    parser.add_argument('-TSX', '--TSX', action='store_true',
                        help='TSX mode: set baseDir to ../ and disable 10 km padding')
    args = parser.parse_args()

    dem = args.dem
    baseDir = '../' if args.TSX else './'
    pad = 0 if args.TSX else 10

    # Locate project.yaml one level up from the track directory
    projectDir = os.path.dirname(os.getcwd())
    framePattern = '*'
    for yamlName in ('project.yaml', 'sensor.yaml'):
        yamlPath = os.path.join(projectDir, yamlName)
        if os.path.isfile(yamlPath):
            if yamlName == 'sensor.yaml':
                u.mywarning('sensor.yaml found but project.yaml expected — '
                            'rename to project.yaml')
            with open(yamlPath) as f:
                sensorData = yaml.safe_load(f)
            framePattern = sensorData.get('framePattern', '*')
            if dem is None:
                regionPath = sensorData.get('regionFile') or sensorData.get('region')
                if regionPath:
                    for candidate in [regionPath, regionPath + '.yaml',
                                      os.path.join(projectDir,
                                                   os.path.basename(regionPath) + '.yaml')]:
                        if os.path.exists(candidate):
                            with open(candidate) as rf:
                                dem = yaml.safe_load(rf).get('dem')
                            break
            break

    framePrefix = framePattern.split('?')[0] if '?' in framePattern else ''

    velFiles = u.dols('ls -d *_*/velocity')
    frameDirs = u.dols('ls -d velocityStats/*-*')
    print(velFiles)
    print(frameDirs)

    for frameDir in frameDirs:
        frameDir = frameDir.split('/')[-1]
        frame1, frame2 = (int(s[len(framePrefix):]) for s in frameDir.split('-'))
        print(frame1, frame2)

        xMin = yMin = 1e8
        xMax = yMax = -1e8
        lastDxDy = None

        for velFile in velFiles:
            d = velFile.split('/')[0]
            parts = d.split('_')
            if len(parts) < 2:
                print(f'Skipping {velFile}: unexpected directory format')
                continue
            frameStr = parts[1][len(framePrefix):]
            if not frameStr.lstrip('-').isdigit():
                print(f'Skipping {velFile}: cannot parse frame from "{parts[1]}"')
                continue
            frame = int(frameStr)
            if os.path.exists(os.path.join(d, 'Exclude')):
                print(f'\033[1m--- Skipping {velFile} because Exclude file found\033[0m')
                continue
            if frame1 <= frame <= frame2:
                bounds = _boundsFromVelDir(velFile)
                if bounds is None:
                    print(f'Skipping {velFile}: no mosaicOffsets.vrt or .vx.geodat')
                    continue
                x1, y1, x2, y2, dx, dy = bounds
                print(x1, y1, x2, y2)
                if x1 < xMin:
                    xMin = x1
                if y1 < yMin:
                    yMin = y1
                if x2 > xMax:
                    xMax = x2
                if y2 > yMax:
                    yMax = y2
                lastDxDy = (dx, dy)

        if lastDxDy is None:
            print(f'No velocity files found for frame range {frameDir}, skipping')
            continue

        print(xMin, yMin, xMax, yMax)
        dx, dy = lastDxDy
        x0 = round(xMin) - pad
        y0 = round(yMin) - pad
        xs = round(xMax + pad - x0) + pad
        ys = round(yMax + pad - y0) + pad
        print('resolution  = "', x0, y0, xs, ys, dx, dy, '"')
        writeVelThumb(frameDir, x0, y0, dx, dy, xs, ys, baseDir, dem=dem)


if __name__ == '__main__':
    main()
