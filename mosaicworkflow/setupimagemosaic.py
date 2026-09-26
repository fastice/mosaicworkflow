#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 13:42:20 2020

@author: ian
"""

import sarfunc as s
import argparse
import utilities as u
import numpy as np
import os
from datetime import datetime, timedelta
from subprocess import call
import mosaicfunc as mosf
import threading
import shapefile
import multiprocessing
import time
import glob
import shlex
from osgeo import gdal
from PIL import Image
from mosaicworkflow.setupquarters import addLabel, messageLines

currentVersion = 4
subVersion = 0

suffixDict = {'uncalibrated': ['image'], 'calibrated': ['sigma0', 'gamma0'],
              'corr': ['corr']}
noDataDict = {'uncalibrated': 0, 'calibrated': -30.0, 'corr': -2.0e9}
# geomosaic -int16 stores calibrated dB as Int16 dB*100 (scale 0.01), so -30 dB -> -3000
INT16_NODATA = -3000


def calSuffixes(calFlag, calOutput='both'):
    '''
    Products written for calFlag; for calibrated output calOutput selects
    sigma0, gamma0, or both (geomosaic -calOutput).
    '''
    if calFlag == 'calibrated' and calOutput in ('sigma0', 'gamma0'):
        return [calOutput]
    return suffixDict[calFlag]


def setupImageMosaicArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\033[1mSetup and run an image mosaics product \033[0m',
        epilog='Notes:\nPart of the mosaicworkflow package.', allow_abbrev='False')
    parser.add_argument('--firstdate', type=str, default='2000-01-01',
                        help='Use central dates >= first date [2000-01-01]')
    parser.add_argument('--lastdate', type=str, default=None,
                        help='Lastdate specify in place of nDays [None]')
    parser.add_argument('--nDays', type=int, default=6,
                        help='Use nDays interval starting at firstdate [6]')
    parser.add_argument('--noReprocess', action='store_true', default=False,
                        help='No proc, just reformat data in release'
                        '- has no effect on image data')
    parser.add_argument('--noMosaic', action='store_true', default=False,
                        help='Do not rerun mosaics - reassemble from pieces')
    parser.add_argument('--noClean', action='store_true', default=False,
                        help='Save the tiles in the the pieces directory')
    parser.add_argument('--setupOnly', action='store_true', default=False,
                        help='Do all steps but run - can be used to clean')
    parser.add_argument('--calibrate', action='store_true', default=False,
                        help='Make an S1 calibrated image product')
    parser.add_argument('--corr', action='store_true', default=False,
                        help='Make a NISAR correlation (coherence) mosaic')
    parser.add_argument('--check', action='store_true', default=False,
                        help='Setup command, but do not run')
    parser.add_argument('--sensor', type=str, default='S1',
                        help='Sensor (S1)')
    parser.add_argument('--prefix', type=str, default='GL_S1bks_mosaic',
                        help='Prefix for product name [GL_S1bks_mosaic]')
    parser.add_argument('--template', type=str, default='mosaic.template.yaml',
                        help='template that defines mosaic')
    parser.add_argument('--ascendingOnly', action='store_true', default=False,
                        help='Only include ascending data')
    parser.add_argument('--descendingOnly', action='store_true', default=False,
                        help='Only include descending data')
    parser.add_argument('--removePad', type=int, default=5,
                        help='Remove first and last n cols')
    parser.add_argument('--averageAll', action='store_true', default=False,
                        help='Average all scenes in date range; omit orbit '
                        'priority flag from geomosaic')
    parser.add_argument('--min', action='store_true', default=False,
                        help='Mosaic = pixel-wise minimum across all inputs')
    parser.add_argument('--max', action='store_true', default=False,
                        help='Mosaic = pixel-wise maximum across all inputs')
    parser.add_argument('--int16', action='store_true', default=False,
                        help='Calibrated only: write sigma0/gamma0 as Int16 '
                        'dB*100 (scale 0.01 in the file; lossless, ~2x '
                        'smaller). Template key: int16')
    parser.add_argument('--calOutput', type=str, default=None,
                        choices=['sigma0', 'gamma0', 'both'],
                        help='Calibrated only: products to make [both]. '
                        'Template key: calOutput')
    parser.add_argument('--featherDistance', type=float, default=None,
                        help='Feather length in metres (converted to output '
                        'pixels for geomosaic -fl); implies --averageAll '
                        'since geomosaic cannot feather with orbit priority '
                        '[0]. Template key: featherDistance')
    parser.add_argument('--gcov', type=str, default=None,
                        help='geomosaic -gcov yaml listing NISAR GCOV HDF5 '
                        'inputs (supplement or replace range/Doppler '
                        'products). Template key: gcov')
    parser.add_argument('--rangeSelect', type=str, default=None,
                        choices=['none', 'near', 'far'],
                        help='Keep only inputs whose incidence angle is within '
                        'angleTolerance of the per-pixel minimum (near) or '
                        'maximum (far) [none]. Template key: rangeSelect')
    parser.add_argument('--angleTolerance', type=float, default=None,
                        help='Tolerance in degrees for rangeSelect [1.0]. '
                        'Template key: angleTolerance')
    parser.add_argument('--minIncidence', type=float, default=None,
                        help='Drop data below this incidence angle (deg). Hard '
                        'gate: not rescued by the fallback. Template key: minIncidence')
    parser.add_argument('--maxIncidence', type=float, default=None,
                        help='Drop data above this incidence angle (deg). '
                        'Template key: maxIncidence')
    parser.add_argument('--angleRamp', action='store_true', default=False,
                        help='Fade each input out linearly over angleTolerance '
                        'instead of a hard cut. Template key: angleRamp')
    parser.add_argument('--angleStride', type=int, default=None,
                        help='Output pixels per cell in the rangeSelect '
                        'pre-pass [10]. Template key: angleStride')
    parser.add_argument('--nParallel', type=int, default=None,
                        help='Tiles run concurrently [min(nTiles, cpus/2)]. '
                        'Template key: nParallel')
    parser.add_argument('--ompThreads', type=int, default=None,
                        help='geomosaic -ompThreads per tile [geomosaic '
                        'default]. Template key: ompThreads')
    parser.add_argument('--message', type=str, default=None,
                        help='Free-text message stored in the release COG '
                        'metadata (key "message"). Template key: message')
    parser.add_argument('--cycleTag', type=str, default=None,
                        help='NISAR cycle label appended to prefix, e.g. '
                        'cycle016 (zero-padded; set by makeimagemosaics in '
                        'nisarCycle mode)')
    args = parser.parse_args()
    #
    if args.descendingOnly and args.ascendingOnly:
        u.myerror('descendingOnly and ascendingOnly flags mutually exclusive')
    #
    dem, smoothL, DBPaths, tracks, prefix, geodat, removePad, passType, _ =\
        mosf.mosaicTempArgs(args.template)
    if args.removePad != 5:
        removePad = args.removePad
    #
    firstDate = datetime.strptime(args.firstdate, "%Y-%m-%d")
    if args.lastdate is None:
        lastDate = firstDate + timedelta(days=args.nDays - 1)
    else:
        lastDate = datetime.strptime(args.lastdate, "%Y-%m-%d")
    #
    if args.corr:
        calFlag = 'corr'
    elif args.calibrate:
        calFlag = 'calibrated'
    else:
        calFlag = 'uncalibrated'
    prodPrefix = args.prefix
    if prefix is not None:
        prodPrefix = prefix
    if args.cycleTag is not None:
        prodPrefix = f'{prodPrefix}_{args.cycleTag}'
    #
    # Some args are hardwired for now but can be updated later
    if tracks is None:
        tracks = ['26', '74', '90', '112', '141', '170', '83']
    # CLI flags win; template passType applies only when neither flag is set
    ascendingOnly = args.ascendingOnly
    descendingOnly = args.descendingOnly
    if not ascendingOnly and not descendingOnly:
        if passType == 'ascending':
            ascendingOnly = True
        elif passType == 'descending':
            descendingOnly = True
    # geoMosaicMode: CLI flag wins over template key
    templateDict = mosf.readTemplate(args.template)
    geoMosaicMode = templateDict.get('geoMosaicMode', None)
    if args.min:
        geoMosaicMode = 'min'
    elif args.max:
        geoMosaicMode = 'max'
    # Make sure no reprocessing, including mosaic3d
    if args.noReprocess:
        args.noMosaic = True
    # geomosaic output options: CLI wins over template keys
    int16 = args.int16 or bool(templateDict.get('int16', False))
    calOutput = args.calOutput or templateDict.get('calOutput', 'both')
    if calOutput not in ('sigma0', 'gamma0', 'both'):
        u.myerror(f'calOutput must be sigma0, gamma0, or both, not {calOutput}')
    featherDistance = args.featherDistance if args.featherDistance is not None \
        else float(templateDict.get('featherDistance', 0))
    gcov = args.gcov or templateDict.get('gcov', None)
    message = args.message or templateDict.get('message', None)
    nParallel = args.nParallel or templateDict.get('nParallel', None)
    rangeSelect = args.rangeSelect or templateDict.get('rangeSelect', 'none')
    if rangeSelect not in ('none', 'near', 'far'):
        u.myerror(f'rangeSelect must be none, near or far, not {rangeSelect}')
    angleTolerance = args.angleTolerance if args.angleTolerance is not None \
        else templateDict.get('angleTolerance', None)
    angleStride = args.angleStride or templateDict.get('angleStride', None)
    angleRamp = args.angleRamp or bool(templateDict.get('angleRamp', False))
    minIncidence = args.minIncidence if args.minIncidence is not None \
        else templateDict.get('minIncidence', None)
    maxIncidence = args.maxIncidence if args.maxIncidence is not None \
        else templateDict.get('maxIncidence', None)
    if rangeSelect == 'none' and (angleTolerance is not None or
                                  angleStride is not None or angleRamp):
        u.myerror('angleTolerance/angleStride require rangeSelect near or far')
    ompThreads = args.ompThreads or templateDict.get('ompThreads', None)
    if gcov is not None:
        gcov = os.path.abspath(gcov)  # geomosaic runs in the pieces directory
        if not os.path.exists(gcov):
            u.myerror(f'gcov yaml {gcov} does not exist')
    if (int16 or calOutput != 'both') and calFlag != 'calibrated':
        u.myerror('int16 and calOutput apply only to --calibrate mosaics')
    featherPixels = int(round(featherDistance / templateDict['dx'])) \
        if featherDistance > 0 else 0
    # min/max mode always averages over the full date range (no orbit priority)
    # geomosaic refuses -fl with orbit priority, so feathering also averages
    # rangeSelect also forces averaging: geomosaic refuses orbit priority with it, because the
    # priority modes key off the scale buffer that the selection's fallback test reads
    averageAll = args.averageAll or geoMosaicMode in ('min', 'max') or \
        featherPixels > 0 or rangeSelect in ('near', 'far')
    myArgs = {'prefix': prodPrefix,
              'firstDate': firstDate, 'lastDate': lastDate,
              'noReprocess': args.noReprocess, 'check': args.check,
              'calFlag': calFlag, 'sensor': args.sensor,
              'dateDBs': DBPaths, 'geodatFile': geodat,
              'tracks': tracks, 'noMosaic': args.noMosaic,
              'templateFile': args.template,
              'ascendingOnly': ascendingOnly,
              'descendingOnly': descendingOnly,
              'removePad': removePad,
              'dem': dem, 'smoothL': smoothL, 'cleanUp': not args.noClean,
              'runMosaics': not args.setupOnly,
              'averageAll': averageAll,
              'geoMosaicMode': geoMosaicMode,
              'int16': int16, 'calOutput': calOutput,
              'featherPixels': featherPixels, 'gcov': gcov,
              'message': message, 'nParallel': nParallel,
              'ompThreads': ompThreads, 'rangeSelect': rangeSelect,
              'angleTolerance': angleTolerance, 'angleStride': angleStride,
              'angleRamp': angleRamp, 'minIncidence': minIncidence,
              'maxIncidence': maxIncidence}
    return myArgs


def padInputFiles(inputFiles, padPixels, dx):
    '''
    Grow each tile's output grid by padPixels on every side (first line of a
    geomosaic input file: x0 y0 xSize ySize dx dy, in km). Feathering weights
    depend on the distance to data edges, which a tile cannot see beyond its
    own grid, so tiles are run with a margin of one feather length and
    cropped back in mergeImageTiles.
    '''
    pad = padPixels * dx / 1000.
    for inputFile in inputFiles:
        with open(inputFile) as fp:
            lines = fp.readlines()
        x0, y0, xs, ys, ddx, ddy = [float(v) for v in lines[0].split()]
        lines[0] = f' {x0 - pad:13.6f} {y0 - pad:13.6f} {xs + 2 * pad:13.6f} ' \
            f'{ys + 2 * pad:13.6f} {ddx:9.4f} {ddy:9.4f}\n'
        with open(inputFile, 'w') as fp:
            fp.writelines(lines)


def gcovFilesUsed(gcovYaml):
    '''
    The GCOV files geomosaic -gcov uses from a yaml: files: items and glob:
    matches (sorted), minus downloads in progress (.aria2 alongside) and
    granules superseded by a higher product counter -- the same rules as
    geoMosaic/gcovMosaic.c.
    '''
    files, inFiles = [], False
    with open(gcovYaml) as fp:
        for line in fp:
            line = line.split('#')[0].strip()
            if not line:
                continue
            if line.startswith('-'):
                if inFiles:
                    files.append(line[1:].split()[0])
                continue
            key, _, value = line.partition(':')
            inFiles = key.strip() == 'files'
            if key.strip() == 'glob':
                files += sorted(glob.glob(value.strip()))
    files = [f for f in files if not os.path.exists(f + '.aria2')]
    newest = {}
    for f in files:
        stem, counter = os.path.basename(f).rsplit('_', 1)
        if stem not in newest or counter > newest[stem][0]:
            newest[stem] = (counter, f)
    return sorted({v[1] for v in newest.values()})


def gcovFilesFromLogs(logDir='log'):
    '''
    GCOV files geomosaic actually mosaicked, from the tile stderr logs: each
    'GCOV <file>' line followed by a 'GCOV window' line (not a Skip). Returns
    None if no tile log mentions a GCOV. Run in the pieces directory.
    '''
    used, found = set(), False
    for log in glob.glob(f'{logDir}/stderr.*'):
        current = None
        with open(log, errors='replace') as fp:
            for line in fp:
                if line.startswith('GCOV /'):
                    current, found = line.split()[1], True
                elif line.startswith('GCOV window') and current is not None:
                    used.add(current)
                    current = None
                elif line.startswith('Skip'):
                    current = None
    return sorted(used) if found else None


def addGCOVShapes(shapeFile, gcovYaml, date1, date2, template):
    '''
    Append one record per GCOV used in the mosaic to the shapefile made by
    makeimageshapefile.py, which only knows range/Doppler products. The files
    are those the tiles actually read (gcovFilesFromLogs); without tile logs
    (e.g. --noMosaic after cleanup) the yaml is resolved instead and limited to
    the date range and the mosaic box. Same fields and conventions as
    makeimageshapefile: lon/lat polygon, TRACK, DATE (start time), ORBIT
    (absolute orbit), SAT (NISAR). Footprint is identification/boundingPolygon.
    '''
    import h5py
    from shapely import wkt as shapelyWkt
    from shapely.geometry import box
    from shapely.ops import transform
    from pyproj import Transformer
    toPS = Transformer.from_crs(4326, template['epsg'], always_xy=True).transform
    x0, y0 = template['xll'], template['yll']
    mosaicBox = box(x0, y0, x0 + template['sx'] * template['dx'],
                    y0 + template['sy'] * template['dy'])
    records = []
    fromLogs = gcovFilesFromLogs()
    for f in (fromLogs if fromLogs is not None else gcovFilesUsed(gcovYaml)):
        with h5py.File(f, 'r') as h:
            idn = h['science/LSAR/identification']
            start = datetime.strptime(
                idn['zeroDopplerStartTime'][()].decode()[:19],
                '%Y-%m-%dT%H:%M:%S')
            poly = shapelyWkt.loads(idn['boundingPolygon'][()].decode())
            if fromLogs is None and (start < date1 or start > date2 +
                                     timedelta(days=1) or not
                                     transform(toPS, poly).intersects(mosaicBox)):
                continue
            records.append(([[(x, y) for x, y, *_ in poly.exterior.coords]],
                            {'TRACK': str(int(idn['trackNumber'][()])),
                             'DATE': start.strftime('%Y-%m-%d %H:%M:%S'),
                             'ORBIT': str(int(idn['absoluteOrbitNumber'][()])),
                             'SAT': 'NISAR'}))
    # rewrite with the existing (range/Doppler) records plus the GCOV ones
    old = shapefile.Reader(f'{shapeFile}.shp')
    fields, oldShapes, oldRecords = old.fields[1:], old.shapes(), old.records()
    old.close()
    # explicit names: product names end in '.0', which pyshp would take as the extension
    w = shapefile.Writer(shp=f'{shapeFile}.shp', shx=f'{shapeFile}.shx',
                         dbf=f'{shapeFile}.dbf', shapeType=shapefile.POLYGON)
    for field in fields:
        w.field(*field)
    names = [field[0] for field in fields]
    for shape, record in zip(oldShapes, oldRecords):
        w.poly([shape.points])
        w.record(*record)
    for parts, rec in records:
        w.poly(parts)
        w.record(**{k: v for k, v in rec.items() if k in names})
    w.close()
    print(f'addGCOVShapes: added {len(records)} GCOV frames to {shapeFile}.shp '
          f'(from {"tile logs" if fromLogs is not None else gcovYaml})')
    return len(records)


def cropTiles(suffix, padPixels, cropDir='cropped'):
    '''
    Crop the padPixels margin (see padInputFiles) from each tile into cropDir,
    leaving the padded tiles in place so a reassembly (--noMosaic) crops again
    from the originals. Run in the pieces directory. Returns the glob for the
    cropped tiles.
    '''
    os.makedirs(cropDir, exist_ok=True)
    for tile in sorted(glob.glob(f'mosaic.*{suffix}.tif')):
        ds = gdal.Open(tile)
        nx, ny = ds.RasterXSize, ds.RasterYSize
        ds = None
        gdal.Translate(f'{cropDir}/{tile}', tile,
                       srcWin=[padPixels, padPixels, nx - 2 * padPixels,
                               ny - 2 * padPixels],
                       creationOptions=['COMPRESS=DEFLATE', 'TILED=YES'])
    return f'{cropDir}/*{suffix}.tif'


def burnMessage(jpgFile, lines, template):
    '''
    Burn message lines into a quicklook JPG with the velocity-browse credit
    code (setupquarters.addLabel: xLabel/yLabel placement, font auto-fit).
    Nothing is done without xLabel/yLabel in the template.
    '''
    if template is None or 'xLabel' not in template or 'yLabel' not in template:
        print('burnMessage: no xLabel/yLabel in template; message is in the '
              'COG metadata only')
        return
    gray = np.array(Image.open(jpgFile).convert('L'))
    # white text on the grey-scale backscatter image
    rgb = addLabel(np.dstack([gray, gray, gray]), template, labels=lines,
                   fill=255)
    Image.fromarray(rgb.astype(np.uint8), 'RGB').save(jpgFile, quality=99)


def geoMosaicExtraFlags(myArgs):
    '''
    geomosaic flags for the output options (int16, calOutput, feathering,
    GCOV inputs); empty when none are set, so existing runs are unchanged.
    '''
    flags = ''
    if myArgs.get('int16'):
        flags += ' -int16'
    if myArgs.get('calOutput', 'both') != 'both':
        flags += f' -calOutput {myArgs["calOutput"]}'
    if myArgs.get('featherPixels', 0) > 0:
        flags += f' -fl {myArgs["featherPixels"]}'
    if myArgs.get('gcov') is not None:
        flags += f' -gcov {myArgs["gcov"]}'
    if myArgs.get('ompThreads') is not None:
        flags += f' -ompThreads {myArgs["ompThreads"]}'
    if myArgs.get('rangeSelect', 'none') in ('near', 'far'):
        flags += f' -{myArgs["rangeSelect"]}Range'
        if myArgs.get('angleTolerance') is not None:
            flags += f' -angleTolerance {myArgs["angleTolerance"]}'
        if myArgs.get('angleStride') is not None:
            flags += f' -angleStride {myArgs["angleStride"]}'
        if myArgs.get('angleRamp'):
            flags += ' -angleRamp'
    for key, flag in (('minIncidence', '-minIncidence'), ('maxIncidence', '-maxIncidence')):
        if myArgs.get(key) is not None:
            flags += f' {flag} {myArgs[key]}'
    return flags


def runGeoMosaic(inputFile, date1, date2, dem, smoothL, calFlag, srsInfo,
                 removePad, averageAll=False, geoMosaicMode=None,
                 extraFlags=''):
    '''
    Run geomosaic command
    Parameters
    ----------
    inputFile : inputFile.
    date1 : first date.
    date2 : last date.
    dem: dem
    calFlag: 'callibrated' or 'uncalibrated'
    srsInfo: {'epsg': ..., 'wktFile': ....}
    averageAll : bool, optional
        When True, omit the orbit-priority flag so geomosaic averages all scenes.
    extraFlags : str, optional
        Additional geomosaic flags (see geoMosaicExtraFlags).
    Returns
    -------
    None.
    '''
    flags = {'uncalibrated': f'-removePad {removePad} {smoothL} -byteScale '
             '-BSlowerBound 0.53  -BSupperBound 2.4 -BSexponent .2 '
             '-BSscale 0.0011538463',
             'calibrated': f'-removePad {removePad} {smoothL} -S1Cal ',
             'corr': '-noPower'}[calFlag]
    #
    suffix = {'uncalibrated': '.image', 'calibrated': '', 'corr': '.corr'}[calFlag]
    # suffixes = suffixDict[calFlag]
    outDir = os.path.dirname(inputFile)
    outFile = os.path.basename(inputFile).replace('inputFile', 'mosaic')
    stdoutFile = f'{outDir}/log/{outFile.replace("mosaic","stdout")}'
    stderrFile = f'{outDir}/log/{outFile.replace("mosaic","stderr")}'
    # orbit priority flag: omitted when averaging all scenes
    orbitFlag = '' if averageAll else '-descending'
    # min/max mode flag
    minMaxFlag = f'-{geoMosaicMode}' if geoMosaicMode in ('min', 'max') else ''
    # assemble command
    command = f'cd {outDir}; ' \
        f'geomosaic -GTiff -center {orbitFlag} {minMaxFlag} {flags}{extraFlags} -xyDEM '\
        f'-date1 {date1.strftime("%m-%d-%Y")} ' \
        f'-date2 {date2.strftime("%m-%d-%Y")} '  \
        f'{os.path.basename(inputFile)} {dem} {outFile}{suffix}'
    # Call command
    with open(stdoutFile, 'w') as stdout:
        with open(stderrFile, 'w') as stderr:
            print(command, file=stdout)
            # executable='/bin/csh',
            call(command, shell=True, stdout=stdout, stderr=stderr)
            time.sleep(1)  # Add a delay to see if this fixes random seg faults
            # for suffix in suffixes:
            #     fpImageToTiff(f'{outDir}/{outFile}.{suffix}',
            #                   f'{outDir}/{outFile}.{suffix}',
            #                   srsInfo, calFlag=calFlag)
            #     # make browse
            #     time.sleep(1)


def setupGeoMosaic(inputFiles, myArgs, srsInfo):
    '''
    Setup threads to run a mosaicked image
    Parameters
    ----------
    inputFiles : List of inputFiles to generate submosaics.
    myArgs:  Dictionary of input parameters
    srsInfo: {'epsg': ..., 'wktFile': ....}
    Returns
    -------
    threads : list of threads to run to generate the mosaic.
    '''
    threads = []
    for inputFile in inputFiles:
        threads.append(threading.Thread(target=runGeoMosaic,
                                        args=[inputFile,
                                              myArgs['firstDate'],
                                              myArgs['lastDate'],
                                              myArgs['dem'],
                                              myArgs['smoothL'],
                                              myArgs['calFlag'],
                                              srsInfo,
                                              myArgs['removePad'],
                                              myArgs.get('averageAll', False)],
                                        kwargs={'geoMosaicMode': myArgs.get('geoMosaicMode'),
                                                'extraFlags': geoMosaicExtraFlags(myArgs)}))
    return threads


def fpImageToTiff(imageTile, tiffTile, srsInfo, calFlag='uncalibrated'):
    '''
    Apply standard image scaling to images and make a tiff. Scaling comes
    from maketopstif100.pro
    Parameters
    ----------
    imageTile : geomosaic output.
    tiffTile : name for the tif output.
    srsInfo: {'epsg': ..., 'wktFile': ....}
    Returns
    -------
    None.
    '''
    noData = noDataDict[calFlag]
    myImage = u.geoimage(geoType='scalar', verbose=False)
    # do uncal scaling
    if calFlag == 'uncalibrated':
        suffix = suffixDict['uncalibrated'][0]
        # remove suffix for inputs since geomosaic doesn't add for uncal
        myImage.readData(imageTile.replace(f'.{suffix}', ''))
        lb, ub = 0.53, 2.4
        kZero = myImage.x == 0
        myImage.x = (0.00015 / 0.13 * myImage.x)**0.2  # Scale
        myImage.x = (myImage.x - lb) * (1. / (ub - lb)) * 255.  # map in 0to255
        myImage.x = np.clip(myImage.x, 1, 255).astype(np.uint8)  # clip min 1
        myImage.x[kZero] = 0  # Force all noData values truly zero
    else:
        myImage.readData(imageTile)
    myImage.writeMyTiff(tiffTile, noDataDefault=noData, epsg=srsInfo['epsg'],
                        wktFile=srsInfo['wktFile'], computeStats=False)
    time.sleep(1)  # Add a delay to see if this fixes random seg faults


def makeImageSpatial(corners, outDir, spatialFile):
    '''
    Write spatial file
    Parameters
    ----------
    corners :lat/lon corners of image.
    outDir : main product dir.
    spatialFile : Name of spatial file to create (no path).
    Returns
    -------
    None.
    '''
    u.pushd(f'{outDir}/release')
    # assume makePremet error check
    fp = open(spatialFile, 'w')
    # corner points
    print(f'{corners["ll"]["lon"]:10.5f} {corners["ll"]["lat"]:10.5f}',
          file=fp)
    print(f'{corners["ul"]["lon"]:10.5f} {corners["ul"]["lat"]:10.5f}',
          file=fp)
    print(f'{corners["ur"]["lon"]:10.5f} {corners["ur"]["lat"]:10.5f}',
          file=fp)
    print(f'{corners["lr"]["lon"]:10.5f} {corners["lr"]["lat"]:10.5f}',
          file=fp)
    #
    fp.close()
    u.popd()


def imageSatTypes(shapeFileName):
    '''
    Read the sat types (S1A, S1B, S1C, S1D, NISAR) from the summary shape file.
    Assumes being called from pieces directory.
    Parameters
    ----------
    shapeFilesName : shapefile name (xxx.shp).
    Returns
    -------
    sTypes : {'S1A': True/False, 'S1B': True/False, 'S1C': True/False,
              'S1D': True/False, 'NISAR': True/False}.

    '''
    u.pushd('../release')
    # any records, record LS8 present
    sTypes = {'S1A': False, 'S1B': False, 'S1C': False, 'S1D': False,
              'NISAR': False}
    #
    current = shapefile.Reader(shapeFileName)
    myDict = dict(zip([x[0] for x in current.fields],
                      range(0, len(current.fields))))
    for myRecord in current.records():
        sTypes[myRecord[myDict['SAT']-1]] = True  # Not sure why -1 needed.
    current = []
    u.popd()
    return sTypes


def populateImagePremet(preMetData, firstDate, lastDate, validPlatforms={},
                        prodType=''):
    '''
    populate premet data fields

    Parameters
    ----------
    preMetData : list of prexisting or blank ([]) premet data.
    firstDate : first date in mosaic.
    lastDate : last date in mosaic.
    validPlatforms : , optional
        dict of platforms (e.g.,{"S1A": True, "S1B": False}. The default is {}.
    prodType : Not used for now, optional. The default is ''.
    Returns
    -------
    preMetData : list with preMetData.
    '''
    preMetData.append(s.premet({'VersionID_local': f'{currentVersion}'}))
    preMetData.append(s.premet(
        {'Begin_date': f'{firstDate.strftime("%Y-%m-%d")}'}))
    preMetData.append(s.premet(
        {'End_date': f'{lastDate.strftime("%Y-%m-%d")}'}))
    preMetData.append(s.premet({'Begin_time': '00:00:01.000'}))
    preMetData.append(s.premet({'End_time': '23:59:59.0000'}))
    #
    sensorShort = {'S1A': 'C-SAR', 'S1B': 'C-SAR',
                   'S1C': 'C-SAR', 'S1D': 'C-SAR',
                   'NISAR': 'L-SAR'}
    instShort = {'S1A': 'C-SAR', 'S1B': 'C-SAR',
                 'S1C': 'C-SAR', 'S1D': 'C-SAR',
                 'NISAR': 'NISAR'}
    instrument = {'S1A': 'Sentinel-1A', 'S1B': 'Sentinel-1B',
                  'S1C': 'Sentinel-1C', 'S1D': 'Sentinel-1D',
                  'NISAR': 'NISAR'}
    #
    sensor = ''
    for platform in validPlatforms.keys():
        if validPlatforms[platform]:
            myPlatform = {'AssociatedPlatformShortname': instrument[platform],
                          'AssociatedInstrumentShortname': instShort[platform],
                          'AssociatedSensorShortname': sensorShort[platform]}
            sensor += platform + ' '
            preMetData.append(
                s.premet(myPlatform,
                         container='AssociatedPlatformInstrumentSensor'))
    #
    dataSetID = {'DataSetId': 'MEaSUREs Greenland '
                 f'{prodType} image mosaics from {sensor.strip()}'}
    # force this to be first
    preMetData = [s.premet(dataSetID, container='Collection')] + preMetData
    return preMetData


def setupMosaicDir(date1, date2, prefix):
    '''
    Setup directory structure (middle date) and makes a pieces and release dir
    underneath
    Parameters
    ----------
    date1 : First date.
    date2 : Second Date.
    Returns
    -------
    dirName : Top level product dir.
    '''
    dirName = \
        f'{prefix}.{date1.strftime("%Y-%m-%d")}.{date2.strftime("%Y-%m-%d")}'
    mosf.mkMyDir(dirName, 'Product', verbose=True)
    mosf.mkMyDir(f'{dirName}/pieces', 'pieces', verbose=True)
    mosf.mkMyDir(f'{dirName}/pieces/log', 'log', verbose=True)
    mosf.mkMyDir(f'{dirName}/release', 'release', verbose=True)
    return dirName


def makeProdName(prefix, date1, date2):
    ''' Makek product name'''
    prodName = f'{prefix}_{date1.strftime("%d%b%y")}_' \
        f'{date2.strftime("%d%b%y")}_*_RESm_v{currentVersion:02d}.{subVersion}'
    return prodName


def makeShapePremet(myShapeName, sTypes, firstDate, lastDate):
    ''' write premet data for shapefile '''
    preMetFile = myShapeName.replace('.shp', '.premet')
    preMets = []
    preMets = populateImagePremet(preMets, firstDate, lastDate,
                                  validPlatforms=sTypes,
                                  prodType='shape file meta data for')
    prodName = myShapeName.split('/')[-1]
    # now write the file
    fp = open(preMetFile, 'w')
    # forcee file name to the beginning
    preMets = [preMets[0],
               s.premet({'ProducerGranuleId': prodName.split('.shp')[0]},
                        container='DataGranule')] + preMets[1:]
    for preMet in preMets:  # Output results
        preMet.printParams(fp=fp)
    fp.close()


def mergeImageTiles(prodDir, prefix, date1, date2, posting, corners,
                    cleanUp=True, noReprocess=False, calFlag='uncalibrated',
                    calOutput='both', int16=False, gcov=None, padPixels=0,
                    message=None, template=None):
    '''
    Merge the tifs with a vrt and write the final cloudoptimized geotiff
    Parameters
    ----------
    prodDir : Directory (YYYY-MM-DD) for product.
    date1: First date
    date2: Last date
    posting: x-posting (assume y the same)
    corners:  Corners used to make spo file for shape
    cleanUp : bool, optional
        Remove tiles to clean up pieces directory. The default is True.
    calFlag: str, optional
        'uncalibrated', 'uncalibrated'
    calOutput: str, optional
        sigma0, gamma0, or both (calibrated only)
    int16: bool, optional
        Tiles are Int16 dB*100 (geomosaic -int16): nodata -3000, and a
        predictor for the COG; the 0.01 scale is carried through.
    gcov: str, optional
        GCOV yaml used, if any; marks NISAR as a platform.
    padPixels: int, optional
        Feather margin to crop from each tile before assembly (padInputFiles).
    message: str or list of str, optional
        Stored as dataset metadata "message" in the release COGs and, when
        the template gives xLabel/yLabel, burned into the JPG quicklook like
        the velocity-browse credit block (setupquarters.addLabel).
    template: dict, optional
        Mosaic template (for xLabel/yLabel).
    Returns
    -------
    satTypes : dict indicating which sats used.
    '''
    # prodTypes = {'uncalibrated': ['uncalibrated'],
    #              'calibrated': ['calibrated', 'calibrated']}[calFlag]
    suffixes = calSuffixes(calFlag, calOutput)
    noData = INT16_NODATA if int16 else noDataDict[calFlag]
    predictor = '-co PREDICTOR=YES ' if int16 else ''
    lines = messageLines(message)
    messageMeta = f'-mo {shlex.quote("message=" + "; ".join(lines))} ' \
        if lines else ''
    jpgRes = '500'
    try:  # Use try so cleanUp won't occur if there is a problem
        u.pushd(directory=f'{prodDir}/pieces')
        #
        prodName = makeProdName(prefix, date1, date2)
        satTypes = {k: False for k in ['S1A', 'S1B', 'S1C', 'S1D', 'NISAR']}
        if calFlag == 'corr':
            # makeimageshapefile.py parses S1/TSX/CSK product names and errors with
            # "unsupported sensor" on a NISAR .nisar.cor line. main() already skips the premets
            # for corr for the same reason, so skip the shapefile too rather than fail the merge.
            satTypes['NISAR'] = True
        else:
            shapeFileName = prodName.replace('*', 'shape').replace('_RESm', '')
            command = f'makeimageshapefile.py inputFile.0.0 ' \
                f'../release/{shapeFileName}'
            print(command)
            call(command, shell=True)  # , executable='/bin/csh')
            if gcov is not None:
                # makeimageshapefile.py only knows range/Doppler inputs
                addGCOVShapes(f'../release/{shapeFileName}', gcov, date1, date2,
                              template)
            satTypes = imageSatTypes(f'{shapeFileName}.shp')
            if gcov is not None:
                satTypes['NISAR'] = True
            makeImageSpatial(corners, '..', f'{shapeFileName}.spo')
            makeShapePremet(f'../release/{shapeFileName}.shp', satTypes,
                            date1, date2)
        # build the vrt
        if not noReprocess:
            for suffix in suffixes:  # Loop on prod types
                # Remove prexisting
                if os.path.exists(f'{prodDir}.{suffix}.vrt'):
                    os.remove(f'{prodDir}.{suffix}.vrt')
                # generate vrt from tiffs (margins cropped if tiles were padded)
                tiles = cropTiles(suffix, padPixels) if padPixels > 0 \
                    else f'*{suffix}.tif'
                command = f'gdalbuildvrt -vrtnodata {noData} '\
                    f'{prodDir}.{suffix}.vrt {tiles}'
                call(command, shell=True)  # , executable='/bin/csh')
                # create the cloud optimized geo
                tifName = \
                    prodName.replace("*", suffix).replace("RES", str(posting))
                jpgName = prodName.replace("*", suffix).replace("RES", jpgRes)
                newTif = f'../release/{tifName}.tif'
                origVrt = f'{prodDir}.{suffix}.vrt'
                command = 'gdal_translate -of COG -co COMPRESS=DEFLATE ' \
                    f'{predictor}-co RESAMPLING=AVERAGE -co OVERVIEWS=IGNORE_EXISTING ' \
                    '-co GEOTIFF_VERSION=1.1 -co BIGTIFF=IF_SAFER -stats ' \
                    '-co NUM_THREADS=ALL_CPUS ' \
                    f'{messageMeta}' \
                    f'{origVrt} {newTif}'
                call(command, shell=True)  # , executable='/bin/csh')
                #
                # make quicklook
                # Int16 input would otherwise give a 12-bit JPEG few readers open
                byteFlag = '-ot Byte ' if int16 else ''
                call(f'gdal_translate -co "QUALITY=99" {byteFlag}-scale -of JPEG -r '
                     f'average -tr {jpgRes} {jpgRes} '
                     f'../release/{tifName}.tif ../release/{jpgName}.jpg',
                     shell=True)  # , executable='/bin/csh')
                if lines:
                    burnMessage(f'../release/{jpgName}.jpg', lines, template)
    except Exception:
        u.myerror('mergeImageTiles: Problem creating final mosaic')
    #
    if cleanUp:
        # remove all image files - leave inputs and vrts
        # dols returns '' (not []) when nothing matches
        toClean = list(u.dols('ls mosaic.*.*')) + glob.glob('cropped/mosaic.*.*')
        for mosaicFile in toClean:
            os.remove(mosaicFile)
    u.popd()
    #
    return satTypes


def makeImagePremets(outDir, prefix, firstDate, lastDate, sTypes, posting,
                     calFlag='uncalibrated', calOutput='both'):
    ''' write the premet file '''
    u.pushd(f'{outDir}/release')
    preMetFiles = []
    suffixes = calSuffixes(calFlag, calOutput)
    #
    for suffix in suffixes:
        prodType = f'{calFlag} {suffix.replace(".","")}'
        prodName = makeProdName(prefix, firstDate, lastDate)
        preMets = populateImagePremet([], firstDate, lastDate,
                                      validPlatforms=sTypes,
                                      prodType=prodType)
        # generate granule file name using prodName as template
        granuleName = \
            prodName.replace("*", suffix).replace("RES", str(posting))
        imageFile = granuleName + '.tif'
        #
        if not os.path.exists(imageFile):
            u.myerror(f'makeImagePremets: could not find release/{imageFile} '
                      'as a template')
        #
        GranuleID = s.premet({'ProducerGranuleId': granuleName},
                             container='DataGranule')
        preMetFile = f'{granuleName}.premet'
        #
        preMetsSuffix = [preMets[0], GranuleID] + preMets[1:]
        # now write the file
        print(f'Writing premet {preMetFile}')
        with open(preMetFile, 'w') as fp:
            for preMet in preMetsSuffix:
                preMet.printParams(fp=fp)
        preMetFiles.append(preMetFile)
    u.popd()
    return preMetFiles


def checkDateRange(myArgs, sarDB):
    '''
    Make sure the requested dates are in date range
    Parameters
    ----------
    myArgs : dict
        Parsed command line args.
    sarDB : sarDB object
        The date base.
    Returns
    -------
    inRange : Bool
        True if in range, False if not

    '''
    return myArgs['firstDate'] >= sarDB.minDate and \
        myArgs['lastDate'] <= (sarDB.maxDate+timedelta(days=6))


def checkGdalVersion():
    ''' Make sure running at least gdal 3.1 '''
    version = float('.'.join(gdal.VersionInfo().split('0')[0:2]))
    if version < 3.1:
        u.myerror('Gdal must be at least V3.1')
    print(f'Gdal version {version}')


def main():
    ''' Create an image mosaic '''
    # get args
    myArgs = setupImageMosaicArgs()
    checkGdalVersion()
    #
    myDB = mosf.sarDB()
    myDB.readDB(myArgs['dateDBs'])
    # myDB.printByDate()
    #
    if not checkDateRange(myArgs, myDB):
        u.mywarning('Requested range falls outside of DB date range:'
                    f'{myDB.minDate} {myDB.maxDate}\nContinuing with available data')
    # create single dict by date of products across multiple dbs with paths
    dateRangeAll = \
        myDB.extractByDateRange(myArgs['firstDate'],
                                myArgs['lastDate'],
                                ascendingOnly=myArgs['ascendingOnly'],
                                descendingOnly=myArgs['descendingOnly'])

    dateRangeTracks = dateRangeAll.selectTrackDateRangeDB(myArgs['tracks'])
    #
    prodDir = setupMosaicDir(myArgs['firstDate'], myArgs['lastDate'],
                             f'{myArgs["prefix"]}_{myArgs["calFlag"]}')
    #
    writeFiles = not myArgs['noReprocess']
    if myArgs['calFlag'] == 'corr':
        inputFiles, srsInfo, corners, posting = \
            mosf.makeSectionedCorrInputFiles(f'{prodDir}/pieces',
                                             myArgs['templateFile'],
                                             dateRangeTracks,
                                             myArgs['geodatFile'],
                                             writeFiles=writeFiles)
    else:
        inputFiles, srsInfo, corners, posting = \
            mosf.makeSectionedPowerInputFiles(f'{prodDir}/pieces',
                                              myArgs['templateFile'],
                                              dateRangeTracks,
                                              myArgs['geodatFile'],
                                              writeFiles=writeFiles)
    # u.myerror('fsdf')
    # Feathering needs each tile to see data edges beyond its own grid
    if writeFiles and myArgs['featherPixels'] > 0:
        padInputFiles(inputFiles, myArgs['featherPixels'], posting)
    # Determine max threads
    maxThreads = min(len(inputFiles), int(multiprocessing.cpu_count()/2))
    if myArgs['nParallel'] is not None:
        maxThreads = min(len(inputFiles), int(myArgs['nParallel']))
    # setup tiles
    threads = setupGeoMosaic(inputFiles, myArgs, srsInfo)
    # produce tiles
    if myArgs['runMosaics'] and not myArgs['noMosaic']:
        u.runMyThreads(threads, maxThreads, 'geomosaic')
    # merge tiles and produce mosaic in release dir
    print('Merging tiles...')
    satTypes = mergeImageTiles(prodDir, myArgs['prefix'], myArgs['firstDate'],
                               myArgs['lastDate'], posting, corners,
                               calFlag=myArgs['calFlag'],
                               noReprocess=myArgs['noReprocess'],
                               cleanUp=myArgs['cleanUp'],
                               calOutput=myArgs['calOutput'],
                               int16=myArgs['int16'], gcov=myArgs['gcov'],
                               padPixels=myArgs['featherPixels'],
                               message=myArgs['message'],
                               template=mosf.readTemplate(myArgs['templateFile']))
    # Shapefile, premet, and spatial files are S1-specific
    if myArgs['calFlag'] != 'corr':
        print('Making premets...')
        preMetFiles = makeImagePremets(prodDir, myArgs['prefix'],
                                       myArgs['firstDate'],
                                       myArgs['lastDate'], satTypes, posting,
                                       calFlag=myArgs['calFlag'],
                                       calOutput=myArgs['calOutput'])
        for preMetFile in preMetFiles:
            makeImageSpatial(corners, prodDir,
                             preMetFile.replace('premet', 'spo'))


if __name__ == "__main__":
    main()
