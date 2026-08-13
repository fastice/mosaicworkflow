#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
grepnisar - Find status of all NISAR pairs for a year.

NISAR-specific analog of grepdate/greptops. Unlike the S1 workflow, each NISAR
product directory (<orbit>_<frame>) is already a pair: it carries both a
reference geodat (geodat*.geojson) and a secondary geodat (geodat*.secondary.
geojson), so the temporal baseline is intrinsic to the directory and no
cross-directory look-ahead is needed. Run inside a track-<N> directory.

By default only the virtual frames (<orbit>_0000, <orbit>_0001, ...) are listed;
--raw switches to the raw acquisition frames (<orbit>_35, <orbit>_36, ...).

@author: ian
"""
import argparse
import glob
import json
import os
from datetime import datetime, timedelta
from operator import itemgetter

import utilities as u
import yaml


def parseFrames(frame):
    if '*' in frame:
        # all frames
        return 0, 1000000
    try:
        frames = frame.split(',')
        if len(frames) == 1:
            return int(frames[0]), int(frames[0])
        else:
            return int(frames[0]), int(frames[1])
    except Exception:
        u.myerror(f'Could not parse frames: {frame}. Must be NNN or MMM,NNN')


def parseDate(date, dt):
    try:
        return datetime.strptime(date, "%Y-%m-%d") + dt
    except Exception:
        u.myerror(f'Could not parse date: {date}')


def grepnisarArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\033[1mFind status of all NISAR pairs for a year\033[0m',
        epilog='Part of the mosaicworkflow package.')
    parser.add_argument('year', type=int, nargs='?', default=None,
                        help='year [optional; if omitted use --firstdate/'
                        '--lastdate, default 2025-09-01 .. today]')
    parser.add_argument('--firstdate', type=str, default=None,
                        help='first date in range YYYY-MM-DD (overrides year)')
    parser.add_argument('--lastdate', type=str, default=None,
                        help='last date in range YYYY-MM-DD (overrides year)')
    parser.add_argument('--frame', type=str, default='*',
                        help='frame (for single) frame1,frame2 (for range)'
                        ' [*, all frames]')
    parser.add_argument('--raw', action='store_true', default=False,
                        help='list the raw acquisition frames instead of the '
                        'virtual frames (_0000)')
    parser.add_argument('--all', action='store_true', default=False,
                        help='add a column listing the raw frames that make up '
                        'each virtual frame (e.g. 35-40), read from frames.txt')
    args = parser.parse_args()
    frame1, frame2 = parseFrames(args.frame)
    year = args.year
    # Defaults for the date range: a year if one is given, else 2025-09-01..today
    if year is not None:
        if year < 2008 or year > 2040:
            u.myerror('invalid year {0:d}'.format(year))
        defaultFirst = datetime(year, 1, 1)
        defaultLast = datetime(year, 12, 31, 23, 59, 59)
    else:
        defaultFirst = datetime(2025, 9, 1)
        defaultLast = datetime.now()
    # process first date
    if args.firstdate is None:
        date1 = defaultFirst
    else:
        date1 = parseDate(args.firstdate, timedelta(seconds=0))
    # process second date
    if args.lastdate is None:
        date2 = defaultLast
    else:
        date2 = parseDate(args.lastdate, timedelta(seconds=86399.9))
    return date1, date2, frame1, frame2, args.raw, args.all


def getFramePattern():
    '''
    Read the virtual-frame glob (framePattern) from ../project.yaml, the same
    key read by makeframetie.py/autoclean.py/autocleanNISAR.py/updateDateDB.py.
    Running inside a track-<N> dir, project.yaml is one level up. Defaults to
    '00??' if the file or key is absent.
    '''
    projectFile = '../project.yaml'
    if os.path.exists(projectFile):
        with open(projectFile) as f:
            proj = yaml.safe_load(f)
        framePattern = proj.get('framePattern', None)
        if framePattern:
            return str(framePattern)
    return '00??'


def parseGeojsonDate(properties):
    '''Combine Date + NominalTime from a geodat geojson into a datetime.'''
    date = properties.get('Date')
    if date is None:
        return None
    time = properties.get('NominalTime', '0:0:0')
    for fmt in ('%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S'):
        try:
            return datetime.strptime(f'{date} {time}', fmt)
        except Exception:
            continue
    return None


def getDate(product, secondary=False):
    '''
    Get the reference (secondary=False) or secondary (secondary=True) date for
    a product directory from its geodat geojson.
    '''
    if secondary:
        files = sorted(glob.glob(f'{product}/geodat*.secondary.geojson'))
    else:
        files = [f for f in sorted(glob.glob(f'{product}/geodat*.geojson'))
                 if not f.endswith('.secondary.geojson')]
    if len(files) == 0:
        return None
    try:
        with open(files[0], 'r') as fp:
            properties = json.load(fp)['properties']
        return parseGeojsonDate(properties)
    except Exception:
        return None


def getProducts(date1, date2, frame1, frame2, raw, framePattern):
    '''
    Collect the orbit_frame product directories (virtual or raw), filtered by
    frame and reference-date range, with reference and secondary dates parsed.
    Virtual frames are those matching *_<framePattern> (from project.yaml); raw
    frames are the remaining orbit_frame directories.
    '''
    virtualDirs = set(glob.glob(f'*_{framePattern}'))
    if raw:
        allDirs = []
        for pattern in ['????_*', '?????_*', '???_*']:
            allDirs += glob.glob(pattern)
        productDirs = set(allDirs) - virtualDirs
    else:
        productDirs = virtualDirs
    products = []
    for product in sorted(productDirs):
        parts = product.split('_')
        if len(parts) != 2:
            continue
        try:
            orbit, frame = int(parts[0]), int(parts[1])
        except Exception:
            # not an orbit_frame directory, so continue
            continue
        # Filter by frame
        if frame < frame1 or frame > frame2:
            continue
        # Filter by reference date
        date = getDate(product, secondary=False)
        if date is None or date < date1 or date > date2:
            continue
        products.append({'dir': product, 'orbit': orbit, 'frame': frame,
                         'date': date,
                         'secondary': getDate(product, secondary=True)})
    return sorted(products, key=itemgetter('date'))


def slcStatus(product):
    d = product['dir']
    if os.path.exists(f'{d}/{d}.slc'):
        return '+'
    return '-'


def processingStatus(product):
    ''' grepdate's status scheme (x{fast}{thumbs}); the offsets test also
    accepts azimuth.offsets.tif since virtual frames carry the .tif form. '''
    d = product['dir']
    if os.path.exists(f'{d}/azimuth.offsets') or \
            os.path.exists(f'{d}/azimuth.offsets.tif'):
        thumbs = '-'
        if os.path.exists(f'{d}/velocity'):
            thumbs = 'v'
        fast = '.'
        if os.path.exists(f'{d}/fast/azimuth.offsets.fast'):
            fast = 'f'
        return f'x{fast}{thumbs}'
    elif os.path.exists(f'{d}/runboth'):
        return '.--'
    else:
        return 'o--'


def frameRange(product):
    '''
    Compact range of the constituent raw frames from the virtual frame's
    frames.txt (e.g. "35-40", or "35-37,39-40" if there is a gap). Empty
    string if frames.txt is absent (e.g. a raw frame directory).
    '''
    path = f'{product["dir"]}/frames.txt'
    if not os.path.exists(path):
        return ''
    with open(path) as fp:
        nums = sorted(int(x) for x in fp.read().split())
    if len(nums) == 0:
        return ''
    # Compress consecutive frame numbers into runs
    runs = []
    start = prev = nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
            continue
        runs.append((start, prev))
        start = prev = n
    runs.append((start, prev))
    return ','.join(f'{a}-{b}' if a != b else f'{a}' for a, b in runs)


def pairInterval(product):
    ''' Temporal baseline (days) from reference to secondary date. '''
    if product['secondary'] is None:
        return 0
    dt = (product['secondary'] - product['date'] + timedelta(hours=12)).days
    return dt


def processProduct(product, showAll=False):
    ''' Print the summary line for one pair directory. '''
    status = processingStatus(product)
    slc = slcStatus(product)
    dt = pairInterval(product)
    orbFrame = product['dir']
    frames = f'  {frameRange(product)}' if showAll else ''
    print(f'{status} {dt:3}d {orbFrame:>12} '
          f' {product["date"].strftime("%-d-%^b-%Y"):>11}  {slc}  '
          f'{product["date"].strftime("%-j"):>3}{frames}')


def outputSummary(products, showAll=False):
    for product in products:
        processProduct(product, showAll)


def main():
    #
    # Get command line args
    date1, date2, frame1, frame2, raw, showAll = grepnisarArgs()
    #
    framePattern = getFramePattern()
    products = getProducts(date1, date2, frame1, frame2, raw, framePattern)
    #
    # Header (no ascending/descending detection for NISAR)
    kind = 'raw' if raw else 'virtual'
    print(f'\n****  NISAR {kind} frames  ({len(products)} pairs)\n')
    #
    outputSummary(products, showAll)


if __name__ == '__main__':
    main()
