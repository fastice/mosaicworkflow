#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Update the sarDB date database used by the image/velocity mosaic workflows.

The database maps acquisition dates to orbit-frame product directories.
Run this from the datesDB/ directory that sits inside an S1 or NISAR project
root.  On first use the database is created from scratch; on subsequent runs
new directories are merged in, existing entries are left unchanged.

If ``../project.yaml`` exists and contains a ``framePattern`` key, that glob
pattern is used to restrict which orbit-frame directories are scanned (e.g.
``'00??'`` matches ``1667_0000``, ``1667_0052``, etc.).  Explicit
``--firstFrame``/``--lastFrame`` arguments always take priority over the yaml
pattern.

Scanning modes
--------------
  Default (no flag)     Scans for ``../track-*/*_{pattern}/*.pow``
                        Standard S1 binary power images.
  --vrt                 Scans for ``../track-*/*_{pattern}/*.pow.vrt``
                        S1 or NISAR power images referenced via VRT sidecars
                        (e.g. older greenlandProject NISAR products).
  --geojson             Scans for ``../track-*/*_{pattern}/geodat*.geojson``
                        NISAR virtual-frame products that have a GeoJSON
                        geodat but no .pow file (e.g. newGreenlandProject).

  In all modes ``{pattern}`` defaults to ``*`` (pow/vrt) or ``0000``
  (geojson) when no framePattern is found in project.yaml.

Usage examples
--------------
  # S1 binary power images
  cd /path/to/S1project/datesDB
  updateDateDB.py --sensor S1

  # NISAR virtual frames (geojson geodat, no .pow); framePattern read from yaml
  cd /path/to/NISARproject/datesDB
  updateDateDB.py --sensor NISAR80 --geojson

  # Subset of frames (overrides framePattern from yaml)
  updateDateDB.py --sensor S1 --firstFrame 50 --lastFrame 60
"""

import glob
import os
import argparse
import yaml
import mosaicfunc as mosf
import sarfunc as s


def updateDateDBArgs():
    '''Parse command-line arguments for updateDateDB.'''
    parser = argparse.ArgumentParser(
        description='\033[1mUpdate the date database that maps orbit-frame '
                    'product directories to acquisition dates.\033[0m',
        epilog='Run from the datesDB/ directory inside the project root.',
        allow_abbrev=False)
    parser.add_argument('--DBname', type=str, default='dateDataBase',
                        help='Filename for the date database [dateDataBase]')
    parser.add_argument('--sensor', type=str, default='S1',
                        help='Sensor name matching a sensors/*.yaml file '
                             '(e.g. S1, NISAR80) [S1]')
    parser.add_argument('--firstFrame', type=int, default=-1,
                        help='First frame number to include; '
                             '-1 means include all [-1]')
    parser.add_argument('--lastFrame', type=int, default=-1,
                        help='Last frame number to include; '
                             '-1 means include all [-1]')
    parser.add_argument('--vrt', action='store_true', default=False,
                        help='Scan for *.pow.vrt instead of *.pow '
                             '(S1/NISAR with VRT power sidecars)')
    parser.add_argument('--geojson', action='store_true', default=False,
                        help='Scan for geodat*.geojson in orbit-frame dirs '
                             '(NISAR virtual-frame products without .pow)')
    args = parser.parse_args()
    if args.vrt and args.geojson:
        import utilities as u
        u.myerror('--vrt and --geojson are mutually exclusive')
    myArgs = {'DBname': args.DBname,
              'sensor': s.sensorDefinitions(args.sensor),
              'firstFrame': args.firstFrame,
              'lastFrame': args.lastFrame,
              'vrt': args.vrt,
              'geojson': args.geojson}
    return myArgs


def readProjectFramePattern():
    '''
    Read framePattern from ../project.yaml if the file exists.

    Returns
    -------
    str or None
        The framePattern glob string (e.g. ``'00??'``), or None if
        ../project.yaml is absent or does not contain the key.
    '''
    projectFile = '../project.yaml'
    if not os.path.exists(projectFile):
        return None
    with open(projectFile) as f:
        proj = yaml.safe_load(f)
    framePattern = proj.get('framePattern', None)
    if framePattern:
        print(f'Found framePattern: {framePattern!r} in {projectFile}')
    return framePattern


def getGeojsonDirs(frames=None, framePattern=None):
    '''
    Return directories that contain a geodat*.geojson file (NISAR virtual
    frames without .pow files).  Assumes running one level below project root
    (i.e. inside datesDB/).

    Priority for the frame glob:
      1. Explicit ``frames`` list (from --firstFrame/--lastFrame).
      2. ``framePattern`` from ../project.yaml.
      3. Default: ``'0000'`` (virtual-frame directories only).

    Parameters
    ----------
    frames : list of str or None
        Explicit frame suffixes (e.g. ``['0000', '0052']``).
    framePattern : str or None
        Glob pattern for the frame suffix (e.g. ``'00??'``).

    Returns
    -------
    dirs : list of str
        Unique directory paths sorted alphabetically.
    '''
    if frames is not None:
        geodats = []
        for frame in frames:
            geodats += sorted(
                glob.glob(f'../track-*/*_{frame}/geodat*.geojson'))
    elif framePattern is not None:
        geodats = sorted(
            glob.glob(f'../track-*/*_{framePattern}/geodat*.geojson'))
    else:
        geodats = sorted(glob.glob('../track-*/*_0000/geodat*.geojson'))
    return list(dict.fromkeys(os.path.dirname(g) for g in geodats))


def main():
    '''Update the sarDB date database with newly available product directories.'''
    myArgs = updateDateDBArgs()
    print(myArgs)
    # Check for framePattern in project.yaml
    framePattern = readProjectFramePattern()
    # Read (or initialise) the database
    myDB = mosf.sarDB()
    print('Reading existing database...')
    myDB.readDB('.', dataBaseName=myArgs['DBname'])
    # Resolve explicit frame range (takes priority over framePattern)
    frames = None
    if myArgs['firstFrame'] >= 0 and myArgs['lastFrame'] > 0:
        frames = list(range(myArgs['firstFrame'], myArgs['lastFrame'] + 1))
        if myArgs['vrt']:
            frames = [f'000{x}' for x in frames]
    # Discover product directories
    print('Discovering product directories...')
    if myArgs['geojson']:
        myDirs = getGeojsonDirs(frames=frames, framePattern=framePattern)
    elif framePattern is not None and frames is None:
        # pow/vrt mode: use framePattern glob from project.yaml
        suffix = '.vrt' if myArgs['vrt'] else ''
        pows = sorted(
            glob.glob(f'../track-*/*_{framePattern}/*.pow{suffix}'))
        if not pows:
            pows = sorted(
                glob.glob(f'../track-*/*_{framePattern}/*.cor{suffix}'))
        myDirs = list(dict.fromkeys(os.path.dirname(p) for p in pows))
    else:
        myDirs = mosf.getPowerDirs(vrt=myArgs['vrt'], frames=frames)
    print(f'Found {len(myDirs)} directories to process.')
    # Merge new entries into the database
    print('Updating database...')
    myDB.updateDB('.', myDirs, myArgs['sensor'])
    # Write updated database back to disk
    print('Saving database...')
    myDB.saveDB('.', dataBaseName=myArgs['DBname'])
    print('Done.')


if __name__ == '__main__':
    main()
