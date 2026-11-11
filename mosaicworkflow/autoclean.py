#!/usr/bin/env python3
"""
autoclean.py

Cycle through offset directories inside a track-N/ directory, compare each
per-segment velocity map against the velocityStats reference mean and sigma,
flag outlier pixels, map them back to range/azimuth coordinates, and write
bad-pixel index lists that cleanoff can apply.

XY mode (default):  a pixel is bad if vx OR vy deviates from the reference
    by more than sigThresh × sigma.  One list → badoffsets_auto.list.
RA mode (--doRA):   range and azimuth are tested independently.
    Range outliers  → badoffsets_auto.list.dr
    Azimuth outliers → badoffsets_auto.list.da

Run from the track-N/ directory.

Each product gets a summary appended to its velFile.report when cleaning
finishes: how much it had to begin with, how much was cleaned, how much is
left, the mean/median departure from the reference, the coverage of the
culled+interpolated offsets the product starts from, and a one-word status
(OK / SPARSE / NEARLY-EMPTY / EMPTY).  The file's first line is unchanged.

Note that products carrying an Exclude.pending marker are skipped, so their
masks stay frozen at whatever the last non-pending run produced.
"""
import argparse
import yaml
import utilities as u
import numpy as np
import re
import os
from subprocess import call
import scipy.ndimage as morph
import threading
import sarfunc as s
from datetime import datetime
import shapefile
import glob


try:
    from PIL import Image, ImageDraw
except Exception:
    print('could not open PIL - may cause problems')

debug = False


def getSensor():
    """Return the sensor based on the directory name."""
    cwd = os.getcwd()
    if 'TSX' in cwd:
        return 'TSX'
    elif 'Sentinel' in cwd:
        return 'S1'
    elif 'CSK' in cwd:
        return 'CSK'
    u.myerror('cannot parse sensor from path name')


def parseDateArg(dateStr):
    """Parse YYYY:MM:DD → datetime (used as argparse type)."""
    try:
        return datetime.strptime(dateStr, "%Y:%m:%d")
    except ValueError:
        raise argparse.ArgumentTypeError(
            f'Invalid date {dateStr!r}: expected YYYY:MM:DD')


def _readProjectYaml():
    """Return parsed ../project.yaml dict, or {} if not found."""
    for name in ('project.yaml', 'sensor.yaml'):
        path = os.path.join(os.path.dirname(os.getcwd()), name)
        if os.path.exists(path):
            try:
                with open(path) as f:
                    return yaml.safe_load(f) or {}
            except Exception:
                pass
    return {}


def getFramePrefix():
    """Return the non-wildcard prefix from framePattern (e.g. '00' for '00??')."""
    pattern = _readProjectYaml().get('framePattern', '*')
    return pattern.split('?')[0] if '?' in pattern else ''


def getVelocityStatsMode():
    """Return 'RA' or 'XY' from ../project.yaml velocityStatsMode key."""
    return _readProjectYaml().get('velocityStatsMode', 'XY').upper()


def autocleanProcessArgs():
    """Parse command-line arguments and return processing parameters."""
    parser = argparse.ArgumentParser(
        epilog='Part of the mosaicworkflow package.',
        description='Flag outlier offsets by comparison with velocityStats '
                    'reference.  Run from the track-N/ directory.')
    parser.add_argument('-nocull', '--nocull', action='store_true',
                        help='Use the nocull stats instead of culled [False]')
    parser.add_argument('-refresh', '--refresh', action='store_true',
                        help='Overwrite existing bad-offset lists [False]')
    parser.add_argument('-reportOnly', '--reportOnly', action='store_true',
                        help='Write velFile.report from the bad-offset lists '
                             'already on disk; compute and write nothing else. '
                             'Includes Exclude.pending products, which a normal '
                             'run skips [False]')
    parser.add_argument('-remove', '--remove', action='store_true',
                        help='Remove prior bad-offset lists and exit [False]')
    parser.add_argument('-applyall', '--applyall', action='store_true',
                        help='Run cleanoff in all dirs, even those not updated [False]')
    parser.add_argument('-noapply', '--noapply', action='store_true',
                        help='Do not run cleanoff in any directory [False]')
    parser.add_argument('-firstDate', '--firstDate', type=parseDateArg,
                        default=datetime(1900, 12, 31), metavar='YYYY:MM:DD',
                        help='Only process pairs on or after this date [1900:12:31]')
    parser.add_argument('-lastDate', '--lastDate', type=parseDateArg,
                        default=datetime(2100, 12, 31), metavar='YYYY:MM:DD',
                        help='Only process pairs on or before this date [2100:12:31]')
    parser.add_argument('-sigThresh', '--sigThresh', type=float, default=3.0,
                        help='Discard threshold in units of sigma [3.0]')
    parser.add_argument('-threads', '--threads', type=int, default=12,
                        help='Number of parallel threads [12]')
    parser.add_argument('-offDir', '--offDir', default='./*_*',
                        help='Glob pattern for offset directories [./*_*]')
    parser.add_argument('-dem', '--dem', default=None,
                        help='DEM file [auto-selected by date and region]')
    parser.add_argument('-region', '--region', default='greenland',
                        help='Region: greenland or antarctica [greenland]')
    parser.add_argument('-doRA', '--doRA', action='store_true', default=False,
                        help='RA mode: test vr and va independently [auto from project.yaml]')
    parser.add_argument('-doXY', '--doXY', action='store_true', default=False,
                        help='Force XY mode even when project.yaml says RA')
    args = parser.parse_args()

    if args.applyall and args.noapply:
        u.myerror('--applyall and --noapply are mutually exclusive')
    if args.reportOnly and (args.refresh or args.remove or args.applyall):
        u.myerror('--reportOnly cannot be combined with --refresh, --remove '
                  'or --applyall: it must not change anything')
    if args.firstDate > args.lastDate:
        u.myerror(f'firstDate ({args.firstDate}) > lastDate ({args.lastDate})')

    # resolve doRA: explicit flags override yaml
    if args.doXY:
        doRA = False
    elif args.doRA:
        doRA = True
    else:
        doRA = (getVelocityStatsMode() == 'RA')

    region = getRegion(args.region)
    print(region)
    regionDef = s.defaultRegionDefs(region)
    dems = getDems(args.dem, regionDef)

    return (args.nocull, args.refresh, args.applyall, args.noapply,
            args.sigThresh, dems, args.threads, args.offDir,
            args.firstDate, args.lastDate, args.remove, regionDef, doRA,
            args.reportOnly)


def getDems(dem, regionDef):
    if regionDef.name() == 'greenland' and dem is None:
        dems = ['/Volumes/insar7/ian/gimp/gimp1/270m/dem.gimp1.270m',
                '/Volumes/insar7/ian/gimp/gimp2/270m/dem.gimp2.270m']
    elif dem is None:
        dems = [regionDef.dem()] * 2
    else:
        dems = [dem] * 2
    return dems


def getRegion(region):
    """Return region from argument; override from ../region file if present."""
    regionFile = os.path.join(os.path.dirname(os.getcwd()), 'region')
    if os.path.exists(regionFile):
        with open(regionFile) as fp:
            regionFromFile = fp.readline().strip()
        if region != regionFromFile:
            u.mywarning(f'region arg ({region}) does not agree with '
                        f'file value ({regionFromFile}); using file value')
        return regionFromFile
    return region


def _velTiff(basePath, primaryExt):
    """Return True if the GeoTIFF variant exists, False for binary."""
    return os.path.exists(basePath + primaryExt + '.tif')


def getVelRef(noCull, epsg, wktFile, framePrefix='', doRA=False):
    """Get reference velocity (mean and sigma); return dicts keyed by frame."""
    geoType = 'velocityRA' if doRA else 'velocity'
    primaryExt = '.vr' if doRA else '.vx'
    velFile = ['velocity', 'velocity_nocull'][noCull]
    refs = sorted(glob.glob('velocityStats/*-*'))
    sigmas = {}
    means = {}
    ranges = []
    for ref in refs:
        basePath = f'{ref}/{velFile}'
        tiff = _velTiff(basePath, primaryExt)
        checkExt = '.tif' if tiff else ''
        if os.path.exists(basePath + primaryExt + checkExt):
            parts = ref.split('/')[-1].split('-')
            firstFrame = int(parts[0][len(framePrefix):])
            lastFrame  = int(parts[1][len(framePrefix):])
            mean  = u.geoimage(geoType=geoType, verbose=False)
            sigma = u.geoimage(geoType='error', verbose=False)
            mean.readData(basePath, epsg=epsg, wktFile=wktFile, tiff=tiff)
            sigma.readData(basePath, epsg=epsg, wktFile=wktFile, tiff=tiff)
            ranges.append(range(firstFrame, lastFrame + 1))
            for frame in range(firstFrame, lastFrame + 1):
                means[frame] = mean
                sigmas[frame] = sigma
    return means, sigmas, ranges


def getOffsetDirs(noCull, offDirRoot, keepPending=False):
    """Return offset dirs that have azimuth.offsets and no Exclude file.

    keepPending retains Exclude.pending products.  A normal run must skip them
    -- which is also why their masks freeze at whatever the last non-pending
    run produced -- but a report wants exactly those.
    """
    offs = u.globOffsetProducts(f'{offDirRoot}/azimuth.offsets')
    if len(offs) < 1:
        u.myerror('no *_*/azimuth.offsets: in track-XXX dir?')
    velFile = ['velocity', 'velocity_nocull'][noCull]
    offDirs = []
    for off in offs:
        d = off.replace('/azimuth.offsets', '')
        pending = os.path.exists(f'{d}/Exclude.pending')
        if (not os.path.exists(f'{d}/Exclude') and
                (keepPending or not pending) and
                os.path.exists(f'{d}/{velFile}')):
            offDirs.append(d)
        else:
            u.mywarning(f'Skipping {d}: Exclude={os.path.exists(d+"/Exclude")} '
                        f'ExcludePending={os.path.exists(d+"/Exclude.pending")} '
                        f'missingVelFile={not os.path.exists(d+"/"+velFile)}')
    return offDirs


def checkFrame(frame, ranges):
    """Raise NameError if frame is not in any of the ranges."""
    for r in ranges:
        if frame in r:
            return frame
    u.mywarning(f'\033[1m Frame = {frame}  out of range\033[0m')
    raise NameError()


def computeRA(lat, lon, dem, geodatFile):
    """Convert lat/lon to r/a via lltora binary."""
    imDir = '/'.join(geodatFile.split('/')[0:-1])
    templl = imDir + '/temp.ll'
    tempra = imDir + '/temp.ra'
    u.writeLLtoRAformat(lat, lon, tempfile=templl)
    command = f'lltora {geodatFile} {dem} {templl} {tempra} '
    with open('stdout', 'w') as fout, open('stderr', 'w') as ferr:
        call(command, shell=True, executable='/bin/csh', stdout=fout, stderr=ferr)
    r, a = u.readLLtoRA(tempfile=tempra)
    os.remove(templl)
    os.remove(tempra)
    return r, a


def runCleanoff(offDir):
    """Run cleanoff as a thread."""
    with open('stdout', 'w') as fout, open('stderr', 'w') as ferr:
        call(f'cd {offDir} ; cleanoff', shell=True, executable='/bin/csh',
             stdout=fout, stderr=ferr)


def readSigThresh(offDir, sigThresh):
    """Return sigThresh, overridden by offDir/sigThresh.override if present."""
    sigThreshUse = sigThresh
    try:
        overrideFile = f'{offDir}/sigThresh.override'
        if os.path.exists(overrideFile):
            with open(overrideFile) as fp:
                sigThreshUse = float(fp.readline().strip())
    except Exception:
        u.myerror(f'problem with sigThresh file at {offDir}')
    return sigThreshUse


def getVelAndFiles(offDir, epsg, wktFile, noCull, mean, sensorInfo, doRA=False):
    """Read velocity, geodat, and open the report file for offDir."""
    geoType = 'velocityRA' if doRA else 'velocity'
    primaryExt = '.vr' if doRA else '.vx'
    vel = u.geoimage(geoType=geoType, verbose=False)
    velFile = 'velocity_nocull/mosaicOffsets'
    try:
        fpReport = open(f'{offDir}/velFile.report', 'w')
    except Exception:
        u.myerror(f'could not open {offDir}/velFile.report: '
                  'verify velocity_nocull exists')
    basePath = f'{offDir}/{velFile}'
    tiff = _velTiff(basePath, primaryExt)
    checkExt = '.tif' if tiff else ''
    primaryFile = basePath + primaryExt + checkExt
    if not os.path.exists(primaryFile):
        u.mywarning(f'\n\033[1m--- Warning no velocity directory for '
                    f'{primaryFile} - consider running makeframetie.py - \033[0m')
        return False
    vel.readData(basePath, epsg=epsg, wktFile=wktFile, tiff=tiff)
    refAttr = 'vr' if doRA else 'vx'
    if getattr(vel, refAttr).shape != getattr(mean, refAttr).shape:
        print(getattr(vel, refAttr).shape, getattr(mean, refAttr).shape)
        u.mywarning(f'\033[1mvelFile = {velFile} shape does not match '
                    f'reference vel \033[0m')
        raise NameError()
    geodatName = f'geodat{sensorInfo["nlooksR"]}x{sensorInfo["nlooksA"]}.in'
    geodatFile = os.path.join(offDir, geodatName)
    if not os.path.exists(geodatFile):
        u.mywarning(f'\033[1m offDir = {offDir}  missing geodatFile \033[0m')
    return vel, geodatFile, fpReport


def makeSigThresh2D(vel, sigThresh, sigShape, doRA=False):
    """Return a float32 array of per-pixel threshold values."""
    nx, ny = vel.geo.sizeInPixels()
    sigThresh2D = Image.new('F', (nx, ny), sigThresh)
    if sigShape is not None:
        if not os.path.exists(sigShape):
            u.myerror(f'Shape file: {sigShape} does not exist')
        shape = shapefile.Reader(sigShape)
        for feature in shape.shapeRecords():
            xyPoly = np.array(feature.shape.points)
            scaleFactor = feature.record[1]
            if len(xyPoly.shape) == 2:
                xi, yi = np.rint(vel.geo.lltoImage(xyPoly[:, 1], xyPoly[:, 0]))
                poly = list(zip(xi, yi))
                ImageDraw.Draw(sigThresh2D).polygon(poly, outline=1,
                                                    fill=sigThresh * scaleFactor)
    return np.array(sigThresh2D, dtype='f4')


def makeInitialCullMask(vel, mean, doRA=False):
    """Return (maskb, goodpoints, kernel1) for the initial outlier pass."""
    refAttr = 'vr' if doRA else 'vx'
    v0 = getattr(vel, refAttr)
    m0 = getattr(mean, refAttr)
    krad = 5
    kernel1 = np.ones((krad, krad), dtype=bool)
    goodpoints = np.logical_and(np.isfinite(v0), np.isfinite(m0))
    maskb = np.zeros(v0.shape, dtype=np.byte)
    maskb[np.isfinite(v0)] = 1
    return maskb, goodpoints, kernel1


def cleanMaskEdges(maskb, kernel1):
    """Open the mask to find ragged edge pixels."""
    maskbOpen = morph.binary_dilation(
        morph.binary_erosion(maskb, structure=kernel1), structure=kernel1)
    edgePoints = np.logical_and(maskb > 0, maskbOpen < 1)
    maskbOpen[:] = 0
    maskbOpen[edgePoints] = 1
    return maskbOpen, edgePoints


def diffStatsFor(names, diffs, sigmas, goodpoints):
    """Per-component (name, meanAbs, medianAbs, signedMedian, biasInSigma).

    The signed median is the one that matters.  A frame whose baseline was fit
    from a handful of tiepoints comes out DISPLACED, not noisy: a couple of bad
    tiepoints bias the fit, the whole field sits off the reference by a near
    constant, and every pixel then trips the sigThresh*sigma test at once --
    cleaning strips the product for a fit error rather than a data error.
    Scatter and bias are indistinguishable in the absolute difference, which is
    all this report used to carry.
    """
    out = []
    if not np.any(goodpoints):
        return out
    for name, d, sig in zip(names, diffs, sigmas):
        v = d[goodpoints]
        med = float(np.median(v))
        s = float(np.median(sig[goodpoints]))
        out.append((name, float(np.average(np.abs(v))),
                    float(np.median(np.abs(v))), med,
                    med / s if s else float('nan')))
    return out


def makeDiffMask(vel, goodpoints, fpReport, mean, sigma, edgePoints,
                 sigThresh2D, kernel1):
    """XY mode: return one closed mask combining vx and vy outliers + edges."""
    dx = np.zeros(vel.vx.shape)
    dy = np.zeros(vel.vx.shape)
    dx[goodpoints] = np.abs(vel.vx[goodpoints] - mean.vx[goodpoints])
    dy[goodpoints] = np.abs(vel.vy[goodpoints] - mean.vy[goodpoints])
    print(np.average(dx), np.average(dy), file=fpReport)
    if debug:
        u.writeImage('testsig', sigThresh2D, 'f4')
        u.writeImage('testsig.dx', dx, 'f4')
        u.writeImage('testsig.dy', dy, 'f4')
        u.writeImage('testsig.ex', sigma.ex, 'f4')
        u.writeImage('testsig.ey', sigma.ey, 'f4')
    badpoints = np.logical_or(dx > sigThresh2D * sigma.ex,
                              dy > sigThresh2D * sigma.ey)
    badpoints = np.logical_or(badpoints, edgePoints)
    mask = np.zeros(vel.vx.shape, dtype=np.byte)
    mask[badpoints] = 1
    sx = np.zeros(vel.vx.shape)
    sy = np.zeros(vel.vx.shape)
    sx[goodpoints] = vel.vx[goodpoints] - mean.vx[goodpoints]
    sy[goodpoints] = vel.vy[goodpoints] - mean.vy[goodpoints]
    stats = diffStatsFor(('dx', 'dy'), (sx, sy), (sigma.ex, sigma.ey),
                         goodpoints)
    return morph.binary_erosion(
        morph.binary_dilation(mask, structure=kernel1),
        structure=kernel1), stats


def makeDiffMasksRA(vel, goodpoints, fpReport, mean, sigma, edgePoints,
                    sigThresh2D, kernel1):
    """RA mode: return (maskR, maskA) — independent closed masks for vr and va."""
    dr = np.zeros(vel.vr.shape)
    da = np.zeros(vel.va.shape)
    dr[goodpoints] = np.abs(vel.vr[goodpoints] - mean.vr[goodpoints])
    da[goodpoints] = np.abs(vel.va[goodpoints] - mean.va[goodpoints])
    print(np.average(dr), np.average(da), file=fpReport)
    if debug:
        u.writeImage('testsig.dr', dr, 'f4')
        u.writeImage('testsig.da', da, 'f4')
        u.writeImage('testsig.er', sigma.er, 'f4')
        u.writeImage('testsig.ea', sigma.ea, 'f4')

    def _closeMask(bad):
        mask = np.zeros(bad.shape, dtype=np.byte)
        mask[bad] = 1
        return morph.binary_erosion(
            morph.binary_dilation(mask, structure=kernel1), structure=kernel1)

    badR = np.logical_or(dr > sigThresh2D * sigma.er, edgePoints)
    badA = np.logical_or(da > sigThresh2D * sigma.ea, edgePoints)
    sr = np.zeros(vel.vr.shape)
    sa = np.zeros(vel.va.shape)
    sr[goodpoints] = vel.vr[goodpoints] - mean.vr[goodpoints]
    sa[goodpoints] = vel.va[goodpoints] - mean.va[goodpoints]
    stats = diffStatsFor(('dr', 'da'), (sr, sa), (sigma.er, sigma.ea),
                         goodpoints)
    return _closeMask(badR), _closeMask(badA), stats


def offsetsLoss(offDir, badIndex):
    """What the bad-point list costs in the OFFSETS domain, as
    {component: (original, cleaned, left, gridPx)}.

    This is the number that matters.  A product starts from its culled+
    interpolated .dr/.da rasters and what mosaic3d finally reads is those
    minus the bad-point list, so cleaning is measured against `original` --
    not against the velocity grid, which spans the whole velocityStats box and
    makes every single frame look like a rounding error.

    `badIndex` maps a component to the flat indices written for it (in XY mode
    the one list applies to both).  The indices are intersected with the valid
    pixels, because the 3x3 dilation that builds them also lands on no-data
    pixels -- so the list is routinely longer than the data it removes.
    Returns {} if nothing can be read.
    """
    try:
        from osgeo import gdal
    except Exception:
        return {}
    out = {}
    for comp in ['dr', 'da']:
        files = [f for f in glob.glob(f'{offDir}/*.offsets.cull.interp.{comp}.vrt')
                 if 'register' not in os.path.basename(f)]
        if len(files) < 1:
            continue
        try:
            # keep the dataset referenced: a band outlives its dataset only
            # as long as something holds the dataset, else gdal throws on use
            ds = gdal.Open(files[0])
            band = ds.GetRasterBand(1)
            data = band.ReadAsArray()
            noData = band.GetNoDataValue()
            good = np.isfinite(data) if noData is None else data != noData
            original = int(np.count_nonzero(good))
            cleaned = 0
            idx = badIndex.get(comp)
            if idx is not None and len(idx) > 0:
                flat = good.reshape(-1)
                inRange = idx[(idx >= 0) & (idx < flat.size)]
                cleaned = int(np.count_nonzero(flat[inRange]))
            out[comp] = (original, cleaned, original - cleaned, int(data.size))
        except Exception as e:
            u.mywarning(f'could not measure {files[0]} for the report: {e}')
    return out


def readExistingList(offDir, doRA):
    """The bad-offset indices already on disk, as {component: array}."""
    out = {}
    pairs = [('dr', '.dr'), ('da', '.da')] if doRA else [('dr', ''), ('da', '')]
    for comp, suffix in pairs:
        f = f'{offDir}/badoffsets_auto.list{suffix}'
        if os.path.exists(f):
            try:
                out[comp] = np.fromfile(f, dtype=np.int32)
            except Exception as e:
                u.mywarning(f'could not read {f}: {e}')
    return out


def readFitTiepoints(offDir):
    """(nUsed, sigma) from the range baseline fit, or (None, None)."""
    for name in ('rBaseline.deltabp.yaml', 'rBaseline.yaml'):
        f = f'{offDir}/motion/{name}'
        if not os.path.exists(f):
            continue
        try:
            t = open(f).read()
            n = re.search(r'^nTiepointsUsed:\s*(\d+)', t, re.M)
            g = re.search(r'^sigma:\s*([-\d.eE+]+)', t, re.M)
            return (int(n.group(1)) if n else None,
                    float(g.group(1)) if g else None)
        except Exception:
            return None, None
    return None, None


def cleanStatus(left, original):
    """One-word verdict on what survives, as a fraction of what the product
    originally had -- never of the grid."""
    if left < 1:
        return 'EMPTY'
    if not original:
        return 'UNKNOWN'
    pct = 100. * left / original
    if pct < 1.:
        return 'NEARLY-EMPTY'
    if pct < 25.:
        return 'SPARSE'
    return 'OK'


def writeCleanReport(fpReport, offDir, doRA, sigThresh, counts, components,
                     diffStats, offLoss, reportOnly=False):
    """Append a per-product summary of what cleaning removed.

    The file's first line is left as it was (the mean |dx|/|dy|), so anything
    reading it by position still works; this block goes after it.
    """
    valid, withRef = counts['valid'], counts['withRef']

    def pct(n, d):
        return f'{100. * n / d:5.1f}%' if d else '  n/a'

    print(f'\n# autoclean summary  {datetime.now():%Y-%m-%d %H:%M:%S}',
          file=fpReport)
    print(f'product          : {os.path.basename(os.path.normpath(offDir))}',
          file=fpReport)
    print(f'mode             : {"RA" if doRA else "XY"}'
          f'{"  (reportOnly: measured from the lists on disk)" if reportOnly else ""}',
          file=fpReport)
    print(f'sigThresh        : {sigThresh}', file=fpReport)

    # what cleaning costs, against the cull.interp product it starts from
    fractions = []
    for comp, (orig, cleaned, left, gridPx) in sorted(offLoss.items()):
        fractions.append((left, orig))
        print(f'offsets .{comp}      : original {orig:>9d} px valid '
              f'({pct(orig, gridPx)} of the {gridPx} px offsets grid)',
              file=fpReport)
        print(f'  cleaned        : {cleaned:>9d} ({pct(cleaned, orig)} of '
              'original)', file=fpReport)
        print(f'  left           : {left:>9d} ({pct(left, orig)} of original)',
              file=fpReport)
    # why those points went: the decision itself is made on the velocity grid
    print('why              :', file=fpReport)
    for name, _, _, nList in components:
        print(f'  badList {name:<7s}: {nList:>9d} indices written (3x3 '
              'dilated, so more than it removes)', file=fpReport)
    if valid < 1:
        print('  velocity       : EMPTY -- no solution for this frame, so no '
              'bias can be measured; any list on disk predates that',
              file=fpReport)
    else:
        print(f'  noReference    : {valid - withRef:>9d} of {valid} solved px '
              '-- no reference, cannot be checked', file=fpReport)
        print(f'  raggedEdge     : {counts["edge"]:>9d} of {valid} solved px '
              '-- always removed', file=fpReport)
    worstBias = 0.
    for name, meanAbs, medAbs, signed, biasSig in diffStats:
        print(f'  absDiff {name:<6s} : mean {meanAbs:9.3f}  median {medAbs:9.3f}'
              '  (m/yr vs reference)', file=fpReport)
        print(f'  bias    {name:<6s} : median signed {signed:9.3f} m/yr '
              f'({biasSig:+5.1f} sigma)', file=fpReport)
        if np.isfinite(biasSig):
            worstBias = max(worstBias, abs(biasSig))
    nTie, fitSigma = readFitTiepoints(offDir)
    if nTie is not None:
        print(f'  rangeFit       : {nTie} tiepoints, sigma {fitSigma}',
              file=fpReport)
    if worstBias >= 1.:
        print(f'  WARNING        : the field is displaced {worstBias:.1f} sigma '
              'from the reference -- a fit error, not a data error, biases every '
              'pixel past sigThresh at once and strips the product',
              file=fpReport)

    if fractions:
        left, orig = min(fractions, key=lambda p: p[0] / p[1] if p[1] else 0.)
        status = cleanStatus(left, orig)
    elif reportOnly:
        # legacy product: no cull.interp vrt to measure against, and reportOnly
        # computes no mask of its own, so there is nothing to judge -- saying
        # EMPTY here would libel every pre-tiff product
        print('  offsetsLoss    : not measurable (no cull.interp vrt -- legacy '
              'product)', file=fpReport)
        status = 'UNKNOWN'
    else:  # offsets rasters unreadable -- fall back to the velocity grid
        left, orig = min(((lf, valid) for _, _, lf, _ in components),
                         default=(0, valid))
        status = cleanStatus(left, orig)
    print(f'status           : {status}', file=fpReport)


def computeBadRA(mask1, vel, geodatFile, dem, offDir):
    """Return (r, a) range/azimuth coordinates of bad pixels."""
    yi, xi = np.where(mask1)
    if len(xi) == 0:
        print(f'No Bad Points for: {offDir}')
        return [None], [None]
    x, y = vel.xx[xi], vel.yy[yi]
    lat, lon = vel.geo.xykmtoll(x, y)
    r, a = computeRA(lat, lon, dem, geodatFile)
    return r, a


def makeRAMask(offDir, r, a):
    """Build a dilated range/azimuth index mask from bad-point coordinates."""
    fileroot = u.globOffsetProducts(f'{offDir}/*.interp.da')
    if len(fileroot) < 1:
        u.myerror(f'offDir = {offDir} missing offset files')
    fileroot = fileroot[0]
    vrtFile = fileroot.replace('.interp', '').replace('.cull.da', '.vrt')
    if os.path.exists(vrtFile):
        off = u.offsets(fileRoot=fileroot, vrtFile=vrtFile, verbose=False)
    else:
        datFile = fileroot.replace('.interp', '').replace('.cull.da', '.dat')
        datFile = datFile.split('/')[-1]
        off = u.offsets(fileRoot=fileroot, datFile=datFile, verbose=False)
    off.readOffsetsDat()
    ro, ao, index = off.slpRAtoOffsetCoords(r, a)
    ramask = np.zeros((off.na, off.nr), dtype=np.byte)
    ramask[ao, ro] = 1
    if debug:
        u.writeImage('testra0.mask', ramask, 'u1')
    kernelra = np.ones((3, 3), dtype=bool)
    ramask = morph.binary_dilation(ramask, structure=kernelra)
    if debug:
        u.writeImage('testra1.mask', ramask, 'u1')
    return np.int32(np.flatnonzero(ramask))


def writeIndexList(offDir, raIndex, suffix=''):
    """Write bad-pixel indices to badoffsets_auto.list[suffix]."""
    listfile = f'{offDir}/badoffsets_auto.list{suffix}'
    with open(listfile, 'w') as fp:
        raIndex.tofile(fp)


def _writeMaskIfBad(mask, vel, geodatFile, dem, offDir, suffix):
    """Helper: convert mask → r/a → RA index → write list.

    Returns the indices written, or None if there were no bad points -- the
    report intersects them with the valid offsets pixels to say what cleaning
    actually costs.
    """
    r, a = computeBadRA(mask, vel, geodatFile, dem, offDir)
    if r[0] is not None:
        raIndex = makeRAMask(offDir, r, a)
        writeIndexList(offDir, raIndex, suffix=suffix)
        return raIndex
    return None


def processOffsets(offDir, noCull, mean, sigma, sigThresh, sigShape, dem, epsg,
                   wktFile, sensorInfo, doRA=False, reportOnly=False):
    """Compute and write bad-offset mask(s) for one directory.

    With reportOnly, nothing is computed or written but the report: the bad
    point lists already on disk are read back and measured, so the numbers
    describe the cleaning the products actually carry.
    """
    result = getVelAndFiles(offDir, epsg, wktFile, noCull, mean, sensorInfo, doRA)
    if result is False:
        return False
    vel, geodatFile, fpReport = result
    sigThresh2D = makeSigThresh2D(vel, sigThresh, sigShape, doRA)
    maskb, goodpoints, kernel1 = makeInitialCullMask(vel, mean, doRA)
    maskbClose, edgePoints = cleanMaskEdges(maskb, kernel1)

    counts = {'total': int(maskb.size),
              'valid': int(np.count_nonzero(maskb)),
              'withRef': int(np.count_nonzero(goodpoints)),
              'edge': int(np.count_nonzero(edgePoints))}

    # count both against the valid pixels only: the closed mask spreads past
    # the data edge under dilation, so a raw mask.sum() can exceed the number
    # of pixels the product ever had.  cleaned + left == goodAtStart this way.
    def _cleaned(mask):
        return int(np.count_nonzero(np.logical_and(maskb > 0, mask > 0)))

    def _left(mask):
        return int(np.count_nonzero(np.logical_and(maskb > 0, mask < 1)))

    if reportOnly:
        badIndex = readExistingList(offDir, doRA)
        diffStats = diffStatsFor(
            ('dr', 'da') if doRA else ('dx', 'dy'),
            [getattr(vel, a) - getattr(mean, a) for a in
             (('vr', 'va') if doRA else ('vx', 'vy'))],
            (sigma.er, sigma.ea) if doRA else (sigma.ex, sigma.ey),
            goodpoints)
        components = [(c, 0, 0, len(v)) for c, v in sorted(badIndex.items())]
        if not doRA:        # one list covers the product; report it once
            components = components[:1]
            if components:
                components = [('offsets', 0, 0, components[0][3])]
    elif doRA:
        maskR, maskA, diffStats = makeDiffMasksRA(vel, goodpoints, fpReport,
                                                  mean, sigma, edgePoints,
                                                  sigThresh2D, kernel1)
        if debug:
            u.writeImage('testBadR', maskR, 'u1')
            u.writeImage('testBadA', maskA, 'u1')
        idxR = _writeMaskIfBad(maskR, vel, geodatFile, dem, offDir, suffix='.dr')
        idxA = _writeMaskIfBad(maskA, vel, geodatFile, dem, offDir, suffix='.da')
        badIndex = {'dr': idxR, 'da': idxA}
        components = [('range', _cleaned(maskR), _left(maskR),
                       0 if idxR is None else len(idxR)),
                      ('azimuth', _cleaned(maskA), _left(maskA),
                       0 if idxA is None else len(idxA))]
    else:
        mask1, diffStats = makeDiffMask(vel, goodpoints, fpReport, mean, sigma,
                                        edgePoints, sigThresh2D, kernel1)
        if debug:
            u.writeImage('testBad1', mask1, 'u1')
        raIndex = None
        r, a = computeBadRA(mask1, vel, geodatFile, dem, offDir)
        if r[0] is not None:
            raIndex = makeRAMask(offDir, r, a)
            writeIndexList(offDir, raIndex)
        # one list covers the whole product in XY mode, so it applies to both
        badIndex = {'dr': raIndex, 'da': raIndex}
        components = [('offsets', _cleaned(mask1), _left(mask1),
                       0 if raIndex is None else len(raIndex))]

    try:
        writeCleanReport(fpReport, offDir, doRA, sigThresh, counts, components,
                         diffStats, offsetsLoss(offDir, badIndex),
                         reportOnly=reportOnly)
    except Exception as e:
        u.mywarning(f'could not write the autoclean summary for {offDir}: {e}')

    fpReport.close()
    return True


def getDem(dems, myDate):
    """Return dems[0] (gimp1) for dates before 2015, dems[1] (gimp2) otherwise."""
    return dems[1] if myDate.year >= 2015 else dems[0]


def _listExists(offDir, doRA):
    """Return True if the expected bad-offset list file(s) exist."""
    if doRA:
        return (os.path.exists(f'{offDir}/badoffsets_auto.list.dr') and
                os.path.exists(f'{offDir}/badoffsets_auto.list.da'))
    return os.path.exists(f'{offDir}/badoffsets_auto.list')


def _removeList(offDir, doRA):
    """Remove existing bad-offset list file(s)."""
    suffixes = ['.dr', '.da'] if doRA else ['']
    for suf in suffixes:
        f = f'{offDir}/badoffsets_auto.list{suf}'
        if os.path.exists(f):
            print(f'removing {f}')
            os.remove(f)


def main():
    """Flag outlier offsets by comparison with velocityStats reference."""
    noCull, refresh, applyFlag, noapply, sigThresh, dems, nThreads, \
        offDirRoot, firstDate, lastDate, remove, regionDef, doRA, \
        reportOnly = autocleanProcessArgs()

    framePrefix = getFramePrefix()

    offDirs = getOffsetDirs(noCull, offDirRoot, keepPending=reportOnly)
    print(f'Number of products culled {len(offDirs)}')
    sensor = getSensor()
    print(sensor)
    sarDef = s.sensorDefinitions(sensor)
    sensorInfo = sarDef.SAR
    epsg, wktFile = regionDef.epsg(), regionDef.wktFile()
    sigShape = regionDef.sigmaShape()
    means, sigmas, ranges = getVelRef(noCull, epsg, wktFile,
                                      framePrefix=framePrefix, doRA=doRA)

    threads, news, myDates = [], [], []
    skip = 0
    for offDir in offDirs:
        geoFile = sorted(glob.glob(f'{offDir}/geodat*x*.in'))[0]
        myGeo = u.geodatrxa(file=geoFile)
        myDates.append(myGeo.datetime)
        if myGeo.datetime < firstDate or myGeo.datetime > lastDate:
            news.append(False)
            continue
        new = False
        if remove:
            _removeList(offDir, doRA)
            news.append(False)
            continue
        try:
            frameStr = offDir.split('_')[-1][len(framePrefix):]
            frame = checkFrame(int(frameStr), ranges)
            if reportOnly:
                if not _listExists(offDir, doRA):
                    skip += 1          # nothing cleaned yet: nothing to report
                else:
                    threads.append(threading.Thread(
                        target=processOffsets,
                        args=[offDir, noCull, means[frame], sigmas[frame],
                              readSigThresh(offDir, sigThresh), sigShape,
                              getDem(dems, myGeo.datetime), epsg, wktFile,
                              sensorInfo, doRA, True]))
            elif refresh or not _listExists(offDir, doRA):
                new = True
                sigThreshUse = readSigThresh(offDir, sigThresh)
                dem = getDem(dems, myGeo.datetime)
                thread = threading.Thread(
                    target=processOffsets,
                    args=[offDir, noCull, means[frame], sigmas[frame],
                          sigThreshUse, sigShape, dem, epsg, wktFile,
                          sensorInfo, doRA])
                threads.append(thread)
                if noapply:
                    new = False
            else:
                skip += 1
            news.append(new)
        except NameError:
            u.mywarning('\033[1;31mError processing \033[0m ' + offDir)

    if skip > 0:
        if reportOnly:
            u.myalert(f'Skipped {skip} products with no bad-offset list')
        else:
            u.myalert(f'Skipped {skip} products because --refresh not set and '
                      'bad-offset list(s) exist')

    if len(threads) > 0:
        u.runMyThreads(threads, nThreads,
                       'Reports' if reportOnly else 'Making Masks')

    if reportOnly:
        print(f'reportOnly: wrote {len(threads)} velFile.report; '
              'nothing else changed')
        return

    if remove and not applyFlag:
        print('remove flag set: no further action')
        return

    threads = []
    for offDir, new, myDate in zip(offDirs, news, myDates):
        try:
            if myDate < firstDate or myDate > lastDate:
                continue
            if applyFlag or (new and
                             _listExists(offDir, doRA) and
                             os.path.exists(offDir + '/cleanoff')):
                thread = threading.Thread(target=runCleanoff, args=[offDir])
                threads.append(thread)
        except NameError:
            u.mywarning('Error processing ' + offDir)

    if len(threads) > 0:
        u.runMyThreads(threads, nThreads, 'Run cleanoff')


if __name__ == '__main__':
    main()
