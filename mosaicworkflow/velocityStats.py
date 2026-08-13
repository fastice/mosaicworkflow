#!/usr/bin/env python3
import argparse
import utilities as u
import numpy as np
import os
import re
import shapefile
import sarfunc as s
from osgeo import gdal

gdal.UseExceptions()

try:
    from PIL import Image, ImageDraw
except Exception:
    print('couldnot open PIL - may cause problems')

from shutil import copyfile


def getVelocityStatsMode():
    """Read velocityStatsMode from ../../project.yaml if present."""
    for levels in ('../../project.yaml', '../../../project.yaml'):
        yamlPath = os.path.join(os.getcwd(), levels)
        if os.path.exists(yamlPath):
            try:
                import yaml
                with open(yamlPath) as fp:
                    proj = yaml.safe_load(fp) or {}
                return proj.get('velocityStatsMode', 'XY').upper()
            except Exception:
                pass
    return 'XY'


def velocityStatsProcessArgs():
    """Parse command-line arguments and return processing parameters."""
    parser = argparse.ArgumentParser(
        description='Compute velocity statistics for culling.  '
                    'Run in a velocityStats/F1-F2 directory.',
        epilog='Part of the mosaicworkflow package.')
    parser.add_argument('-nocull', '--nocull', action='store_true',
                        help='Use velocity_nocull directories [False]')
    parser.add_argument('-sumVz', '--sumVz', action='store_true',
                        help='Sum only the Vz channel [False]')
    parser.add_argument('-velmap', '--velmap', default=None, metavar='FILE',
                        help='Base velocity map file [region default]')
    parser.add_argument('-region', '--region', default=None, metavar='NAME',
                        help='Region: greenland or antarctica [auto]')
    parser.add_argument('-regionFile', '--regionFile', default=None, metavar='FILE',
                        help='Project-specific region yaml (dem/velMap/sigmaShape); '
                             'takes precedence over -region when both given')
    parser.add_argument('-doRA', '--doRA', action='store_true', default=False,
                        help='Process vr/va components instead of vx/vy')
    parser.add_argument('-doXY', '--doXY', action='store_true', default=False,
                        help='Force vx/vy mode even when project.yaml says RA')
    parser.add_argument('-initializeReference', '--initializeReference',
                        action='store_true',
                        help='(Re)build sims/mosaicOffsets.vr basemap')
    args = parser.parse_args()
    # resolve doRA: explicit flags override yaml
    if args.doXY:
        doRA = False
    elif args.doRA:
        doRA = True
    else:
        doRA = (getVelocityStatsMode() == 'RA')
    return (args.nocull, args.velmap, args.sumVz, args.region, args.regionFile, doRA,
            args.initializeReference)


def detectFileFormat(basePath, doRA=False):
    """Return True if tiff velocity files are present, False for binary."""
    primarySuffix = '.vr.tif' if doRA else '.vx.tif'
    binarySuffix  = '.vr'     if doRA else '.vx'
    if os.path.exists(basePath + primarySuffix):
        return True
    if os.path.exists(basePath + binarySuffix):
        return False
    u.myerror(f'No velocity files found at {basePath}')


def getFullMap(velMapFile):
    """ Read processed velocity map as basemap to use when little new data """
    if velMapFile is None:
        return None
    velMap = u.geoimage(geoType='velocity', verbose=False)
    isTiff = velMapFile.lower().endswith(('.tif', '.vrt'))
    velMap.readData(velMapFile, tiff=isTiff)
    velMap.setupInterp()
    return velMap


def getRABasemap(frameDir):
    """Load the vr/va reference basemap from sims/mosaicOffsets."""
    simsBase = os.path.join(frameDir, 'sims', 'mosaicOffsets')
    tiff = detectFileFormat(simsBase, doRA=True)
    velMap = u.geoimage(geoType='velocityRA', verbose=False)
    velMap.readData(simsBase, tiff=tiff)
    velMap.setupInterp()
    return velMap


def parseVelFilePath(velFile):
    ''' parse vel file path to get orbit frame'''
    d = velFile.split("/")[1]
    try:
        orbit, frame = (int(s) for s in d.split('_'))
    except Exception:
        orbit, frame = -1, 1
        u.mywarning(f'Could no parse velFile: {velFile}')
    if os.path.exists('../'+d+'/Exclude'):
        u.mywarning(f'Skipping {velFile} because Exclude file found')
        frame = -1
    elif os.path.exists('../'+d+'/Exclude.pending'):
        # Soft exclude: keep the frame. Its velocity map is a blank (all
        # no-data) mosaic -- tieScript keeps no-solution segs in the input file
        # and mosaic3d's sigma<0 skip gates contribute nothing -- so it adds
        # zero to the stats (findKeepers finds no finite pixels) but supplies a
        # real grid for setupFirstVel, letting velocityStats produce output for
        # ranges where every frame is pending instead of falling back to the
        # thumb-header no-data path (or nothing at all).
        u.mywarning(f'Including {velFile} (Exclude.pending -- blank map, '
                    f'grid only)')
    return orbit, frame


def setupFirstVel(vel, velFile, frameDir, velMap, doRA=False, tiffMode=False,
                  epsg=None, wktFile=None):
    ''' declare and init variables based on first image shape.
    Returns (mu1, mu2, var1, var2, numAvg, velBase, geo) where
    mu1/mu2 are the component accumulators and geo is the geodat object. '''
    cfg1, cfg2 = ('vr', 'va') if doRA else ('vx', 'vy')
    c1 = getattr(vel, cfg1)
    mux, muy = np.zeros(c1.shape), np.zeros(c1.shape)
    varx, vary = np.zeros(c1.shape), np.zeros(c1.shape)
    numAvg = np.zeros(c1.shape)

    # Build velBase on the same grid as the input velocity product
    geoType = 'velocityRA' if doRA else 'velocity'
    velBase = u.geoimage(geoType=geoType, verbose=False)

    if tiffMode:
        velBase.getGeoFile(velFile + '/mosaicOffsets',
                           'greenland', tiff=True)
    else:
        geodatSuffix = '.vr.geodat' if doRA else '.vx.geodat'
        velBase.readGeodat(velFile + '/mosaicOffsets' + geodatSuffix)
    velBase.xyCoordinates()

    xps = np.zeros(c1.shape)
    yps = np.zeros(c1.shape)
    for i in range(len(velBase.yy)):
        xps[i, :] = velBase.xx
    for i in range(len(velBase.xx)):
        yps[:, i] = velBase.yy

    if velMap is not None:
        interp_result = velMap.interpGeo(xps, yps)
        if doRA:
            # velMap is the region's vx/vy reference (RA mode has no vr/va
            # basemap to interpolate -- see main()). Only speed magnitude is
            # meaningful here (rotation-invariant, same value in any basis);
            # vr/va stay NaN so normStats()'s mean blend (iFewData/iNoData)
            # stays untouched -- the mean remains pure data-derived, only the
            # sigma classification bins (iSlow/iMid/iFast/iExtreme) use
            # velBase.v below.
            setattr(velBase, cfg1, np.full(c1.shape, np.nan))
            setattr(velBase, cfg2, np.full(c1.shape, np.nan))
        else:
            setattr(velBase, cfg1, interp_result[0])
            setattr(velBase, cfg2, interp_result[1])
        velBase.v = interp_result[2]
    else:
        setattr(velBase, cfg1, np.full(c1.shape, np.nan))
        setattr(velBase, cfg2, np.full(c1.shape, np.nan))
        velBase.v = np.full(c1.shape, -1.0)

    velBase.v[np.isnan(velBase.v)] = -1.0

    baseOut = frameDir + ('/mosaicBaseRA' if doRA else '/mosaicBase')
    writeMosaicBase(velBase, baseOut, epsg=epsg, wktFile=wktFile)

    velBase.v[velBase.v < 0] = 1000000
    return mux, muy, varx, vary, numAvg, velBase, velBase.geo


def findKeepers(vel, velBase, mean, sigma, sigMask, doRA=False):
    ''' find good values; flag outliers '''
    cfg1, cfg2 = ('vr', 'va') if doRA else ('vx', 'vy')
    eSfx1, eSfx2 = ('er', 'ea') if doRA else ('ex', 'ey')
    if sigma is None:
        vel.v[np.isnan(vel.v)] = -1
        jOutliers = np.logical_and(
            np.logical_and(vel.v > velBase.v * 3, velBase.v > 100),
            sigMask < 1.01)
    else:
        iNoData = np.logical_or(np.isnan(getattr(vel, cfg1)),
                                np.isnan(getattr(mean, cfg1)))
        ex = np.absolute(getattr(vel, cfg1) - getattr(mean, cfg1))
        ey = np.absolute(getattr(vel, cfg2) - getattr(mean, cfg2))
        ex[iNoData], ey[iNoData] = 0, 0
        jOutliers = np.logical_and(
            np.logical_or(ex > 3 * getattr(sigma, eSfx1),
                          ey > 3 * getattr(sigma, eSfx2)),
            sigMask < 1.01)
    getattr(vel, cfg1).__setitem__(jOutliers, np.nan)
    iGood = np.isfinite(getattr(vel, cfg1))
    return iGood


def sumStats(vel, mux, muy, varx, vary, numAvg, iGood, doRA=False):
    ''' accumulate stats '''
    cfg1, cfg2 = ('vr', 'va') if doRA else ('vx', 'vy')
    c1 = getattr(vel, cfg1)
    c2 = getattr(vel, cfg2)
    mux[iGood]    += c1[iGood]
    muy[iGood]    += c2[iGood]
    varx[iGood]   += c1[iGood] ** 2
    vary[iGood]   += c2[iGood] ** 2
    numAvg[iGood] += 1


def normStats(velBase, numAvg, mux, muy, varx, vary, froot, velFileSave,
              sigMask, doRA=False):
    ''' compute normed stats '''
    cfg1, cfg2 = ('vr', 'va') if doRA else ('vx', 'vy')
    rawSfx1 = '.raw.vr' if doRA else '.raw.vx'
    rawSfx2 = '.raw.va' if doRA else '.raw.vy'

    iData, iNoData, iSomeData, iFewData, iSlow, iMid, iFast, iExtreme, \
        baseGood, countThreshF = computePointCounts(velBase, numAvg, doRA=doRA)

    mux[iSomeData] /= numAvg[iSomeData]
    muy[iSomeData] /= numAvg[iSomeData]

    u.writeImage(froot + rawSfx1, mux, '>f4')
    u.writeImage(froot + rawSfx2, muy, '>f4')
    geodatSrc = 'mosaicOffsets' + ('.vr.geodat' if doRA else '.vx.geodat')
    if os.path.exists(velFileSave + '/' + geodatSrc):
        copyfile(velFileSave + '/' + geodatSrc, froot + rawSfx1 + '.geodat')
        copyfile(velFileSave + '/' + geodatSrc, froot + rawSfx2 + '.geodat')

    # blend few-data regions with basemap
    w1 = numAvg[iFewData] / countThreshF
    w2 = (countThreshF - numAvg[iFewData]) / countThreshF
    base1 = getattr(velBase, cfg1)
    base2 = getattr(velBase, cfg2)
    mux[iFewData] = w2 * base1[iFewData] + w1 * mux[iFewData]
    muy[iFewData] = w2 * base2[iFewData] + w1 * muy[iFewData]
    mux[iNoData]  = base1[iNoData]
    muy[iNoData]  = base2[iNoData]

    shapeTmp = mux.shape
    sigx, sigy = np.zeros(shapeTmp), np.zeros(shapeTmp)
    nScale = np.zeros(shapeTmp)
    nScale[iData] = 1. / (numAvg[iData] - 1.0)
    sigx[:], sigy[:] = 18, 18
    sigx[baseGood] = 15. + .2 * np.absolute(base1[baseGood])
    sigy[baseGood] = 15. + .2 * np.absolute(base2[baseGood])
    sigx[iMid] += 15
    sigy[iMid] += 15
    sigx[iData] = np.sqrt(nScale[iData] * varx[iData] - mux[iData] ** 2)
    sigy[iData] = np.sqrt(nScale[iData] * vary[iData] - muy[iData] ** 2)
    sigx[iMid] = np.minimum(sigx[iMid], 0.5 * np.abs(mux[iMid]))
    sigy[iMid] = np.minimum(sigy[iMid], 0.5 * np.abs(muy[iMid]))
    sigx[iMid] = np.maximum(sigx[iMid], velBase.v[iMid] * 0.075 + 7.5)
    sigy[iMid] = np.maximum(sigy[iMid], velBase.v[iMid] * 0.075 + 7.5)
    sigx[iSlow] = np.minimum(sigx[iSlow], 30)
    sigy[iSlow] = np.minimum(sigy[iSlow], 30)
    sigx[iSomeData] = np.maximum(sigx[iSomeData], 5)
    sigy[iSomeData] = np.maximum(sigy[iSomeData], 5)
    sigx[iFast] = np.maximum(sigx[iFast], 75)
    sigy[iFast] = np.maximum(sigy[iFast], 75)
    sigx[iExtreme] = np.maximum(sigx[iExtreme], 75 + 0.05 * velBase.v[iExtreme])
    sigy[iExtreme] = np.maximum(sigy[iExtreme], 75 + 0.05 * velBase.v[iExtreme])
    iHiVar = np.logical_and(sigMask > 1.01, iData)
    if np.max(iHiVar):
        sigx[iHiVar] = np.sqrt(nScale[iHiVar] * varx[iHiVar] - mux[iHiVar] ** 2)
        sigy[iHiVar] = np.sqrt(nScale[iHiVar] * vary[iHiVar] - muy[iHiVar] ** 2)
    return sigx, sigy


def computePointCounts(velBase, numAvg, doRA=False):
    ''' compute counts '''
    cfg1, cfg2 = ('vr', 'va') if doRA else ('vx', 'vy')
    try:
        countThresh  = 7
        countThreshF = float(countThresh)
        iData = numAvg >= countThresh
        base1 = getattr(velBase, cfg1)
        base2 = getattr(velBase, cfg2)
        baseGood = np.logical_and(np.absolute(base1) < 30000,
                                  np.absolute(base2) < 30000)
        iFewData = np.logical_and(
            np.logical_and(numAvg <= countThresh, numAvg > 0), baseGood)
        iFast    = np.greater(velBase.v, 300)
        iExtreme = np.logical_and(np.greater(velBase.v, 1500), baseGood)
        iMid     = np.logical_and(velBase.v < 300, velBase.v > 80)
        iSlow    = np.logical_and(velBase.v <= 80, velBase.v > -0.0001)
        iNoData  = np.logical_and(numAvg < 1, np.isfinite(base1))
        iSomeData = numAvg >= 1
    except NameError:
        u.myerror("Error: Run in velocityStats ?")
    return (iData, iNoData, iSomeData, iFewData, iSlow, iMid, iFast, iExtreme,
            baseGood, countThreshF)


def referenceSigmaFloor(vComp, speed):
    ''' Reference-velocity sigma floor (m/yr) for ONE velocity component,
    reproducing normStats()'s reference-regime (no accumulated data) threshold
    model: the flat 18 base, the component-linear 15 + 0.2*|vComp| term, and the
    speed-classified mid/slow/fast/extreme clamps. This is the piecewise model
    that has driven velocityStats' XY-mode sigma for a long time.

    Factored out so autocleanNISAR can apply the SAME model per-frame in
    range/azimuth: in RA mode velocityStats has no vr/va basemap on its map grid
    (see setupFirstVel), so its map-grid er/ea keep only the isotropic speed-bin
    part -- the component-linear term never fires. autocleanNISAR recovers the
    range/azimuth reference velocity components from the synthetic offsets and
    calls this to rebuild the full, component-aware floor.

    vComp / speed are m/yr arrays (speed = ground-speed magnitude, used only for
    the bin classification). The data-driven branches in normStats (measured std
    where navg>=7, the navg>=1 max-5 clamp, high-variance polygons) are
    deliberately omitted here -- those depend on accumulated data, which the
    caller (autocleanNISAR) supplies instead by blending in the velocityStats
    data sigma with the navg weight. Kept in deliberate lock-step with
    normStats()'s inline model: any change to that floor must be mirrored here. '''
    vComp = np.asarray(vComp, dtype=float)
    speed = np.asarray(speed, dtype=float)
    sig = np.full(speed.shape, 18.0)
    baseGood = np.abs(vComp) < 30000
    sig[baseGood] = 15. + 0.2 * np.abs(vComp[baseGood])
    iMid = np.logical_and(speed < 300, speed > 80)
    iSlow = np.logical_and(speed <= 80, speed > -1e-4)
    iFast = speed > 300
    iExtreme = np.logical_and(speed > 1500, baseGood)
    sig[iMid] += 15
    sig[iMid] = np.minimum(sig[iMid], 0.5 * np.abs(vComp[iMid]))
    sig[iMid] = np.maximum(sig[iMid], 0.075 * speed[iMid] + 7.5)
    sig[iSlow] = np.minimum(sig[iSlow], 30)
    sig[iFast] = np.maximum(sig[iFast], 75)
    sig[iExtreme] = np.maximum(sig[iExtreme], 75 + 0.05 * speed[iExtreme])
    return sig


def writeStats(mux, muy, sigx, sigy, numAvg, froot, geo, doRA=False,
              epsg=None, wktFile=None, computeStats=True):
    ''' write stats results as GeoTIFF + VRT (mean pair, sigma pair, navg) '''
    meanType  = 'velocityRA' if doRA else 'velocity'
    sigmaType = 'errorRA'    if doRA else 'error'
    cfg1, cfg2   = ('vr', 'va') if doRA else ('vx', 'vy')
    eCfg1, eCfg2 = ('er', 'ea') if doRA else ('ex', 'ey')

    mean = u.geoimage(geoType=meanType, verbose=False)
    mean.geo = geo
    setattr(mean, cfg1, mux.astype('f4'))
    setattr(mean, cfg2, muy.astype('f4'))
    mean.writeMyTiff(froot, epsg=epsg, wktFile=wktFile, noV=True,
                     computeStats=computeStats)
    mean.writeMyVrt(froot)

    sigma = u.geoimage(geoType=sigmaType, verbose=False)
    sigma.geo = geo
    setattr(sigma, eCfg1, sigx.astype('f4'))
    setattr(sigma, eCfg2, sigy.astype('f4'))
    sigma.writeMyTiff(froot, epsg=epsg, wktFile=wktFile, noV=True,
                      computeStats=computeStats)
    sigma.writeMyVrt(froot, vrtFile=froot + '.err.vrt')

    navg = u.geoimage(geoType='scalar', verbose=False)
    navg.geo = geo
    navg.x = numAvg.astype('f4')
    navg.writeMyTiff(froot + '.navg', epsg=epsg, wktFile=wktFile,
                     computeStats=computeStats)


def writeMosaicBase(velBase, baseOut, epsg=None, wktFile=None):
    ''' Write the few-data-blend reference basemap as a GeoTIFF+VRT combo
    (per-component .tif + multiband .vrt, same convention as writeStats()),
    instead of the legacy flat-binary + .geodat pair writeData() produced.
    mosaicBase/mosaicBaseRA is a QC-only artifact (not read downstream), so this
    is purely a format change. noV=True: only the vr/va (or vx/vy) components,
    matching writeData()'s old output. computeStats=False: the RA no-data
    basemap is all-noData (vr/va NaN -> -2e9), which would crash GDAL
    GetStatistics. Like writeData(), writeMyTiff() replaces NaN with the same
    -2e9 (_NO_DATA['.vr']/['.va']) in place, so normStats()'s later
    base1/base2 = velBase.vr/.va reads are unchanged. '''
    velBase.writeMyTiff(baseOut, epsg=epsg, wktFile=wktFile, noV=True,
                        computeStats=False)
    velBase.writeMyVrt(baseOut)


def _findFrameDirsInRange(frame1, frame2):
    ''' Return frame dirs in .. (e.g. ../3391_0020) where frame is in [frame1, frame2]. '''
    result = []
    for name in sorted(os.listdir('..')):
        parts = name.split('_')
        if len(parts) != 2:
            continue
        try:
            frame = int(parts[1])
        except ValueError:
            continue
        if frame1 <= frame <= frame2:
            result.append(os.path.join('..', name))
    return result


def _gridFromThumbHeader(frameDir):
    ''' Parse the polar-stereographic map grid for this frame range from the
    tiepoints thumb header's "resolution" line -- "x0 y0 xSizeKm ySizeKm dxKm
    dyKm" (all km, the standard mosaic3d region argument, identical to the grid
    the mosaic velocity for this range would have been produced on). Returns
    (x0Km, y0Km, xs, ys, dxM, dyM) or None if the header/line is missing.

    This is the correct grid for the no-data-path velBase: the earlier version
    took the grid from an excluded frame's range.offsets.tif, whose geotransform
    is in radar (range/azimuth) pixel coordinates, not PS map coordinates -- so
    the no-data velocityStats output landed at a nonsense location that
    autocleanNISAR's PS-footprint interpolation never overlapped, silently
    flagging every offset pixel as bad. '''
    velStatsDir = os.path.dirname(os.path.abspath(frameDir))
    trackDir = os.path.dirname(velStatsDir)
    tag = os.path.basename(frameDir).replace('-', 'dash')
    header = os.path.join(trackDir, 'tiepoints', f'vel_thumb_header_{tag}')
    if not os.path.exists(header):
        return None
    with open(header) as fp:
        text = fp.read()
    m = re.search(r'resolution\s*=\s*"?\s*'
                  r'([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+'
                  r'([-\d.]+)\s+([-\d.]+)', text)
    if m is None:
        return None
    x0, y0, xSizeKm, ySizeKm, dxKm, dyKm = (float(v) for v in m.groups())
    # Guard the "0 0 0 0" auto-size sentinel (or any unfilled/degenerate line):
    # makevelstatsregions.py rewrites these resolutions during
    # --runVelstatsregions, and an un-filled one would otherwise yield a 0-pixel
    # grid. Refuse it so the caller warns + writes noVelocityFiles.txt instead.
    if dxKm <= 0 or dyKm <= 0 or xSizeKm <= 0 or ySizeKm <= 0:
        return None
    xs = int(round(xSizeKm / dxKm))
    ys = int(round(ySizeKm / dyKm))
    if xs <= 0 or ys <= 0:
        return None
    return x0, y0, xs, ys, dxKm * 1000., dyKm * 1000.


def velFileLoop(velFiles, frame1, frame2, frameDir, velMap, noCull, myRegion,
                useSig=False, doVz=False, doRA=False):
    ''' loop over velocity files to accumulate stats, then normalise and save '''
    first = True
    froot = [frameDir + '/velocity', frameDir + '/velocity_nocull'][noCull]
    geoType = 'velocityRA' if doRA else 'velocity'
    errType = 'errorRA'    if doRA else 'error'
    eSfx1   = '.er'        if doRA else '.ex'
    eSfx2   = '.ea'        if doRA else '.ey'

    if useSig:
        print(froot)
        sigma = u.geoimage(geoType=errType, verbose=False)
        mean  = u.geoimage(geoType=geoType, verbose=False)
        # read back the tiff+vrt this frame-range's own first velFileLoop()
        # pass just wrote via writeStats()
        sigma.readData(froot, geoType=errType, tiff=True)
        mean.readData(froot, tiff=True)
    else:
        mean, sigma = None, None

    tiffMode = None  # detect from first real file

    for velFile in velFiles:
        orbit, frame = parseVelFilePath(velFile)
        if frame >= frame1 and frame <= frame2:
            basePath = velFile + '/mosaicOffsets'
            if tiffMode is None:
                tiffMode = detectFileFormat(basePath, doRA=doRA)
            vel = u.geoimage(geoType=geoType, verbose=False)
            vel.readData(basePath, tiff=tiffMode)
            if doVz:
                vZ = u.geoimage(geoType='scalar', verbose=False)
                vZ.readData(velFile + '/mosaicOffsets.vz')

            if first:
                velFileSave = velFile
                mux, muy, varx, vary, numAvg, velBase, outGeo = \
                    setupFirstVel(vel, velFile, frameDir, velMap,
                                  doRA=doRA, tiffMode=tiffMode,
                                  epsg=myRegion.epsg(), wktFile=myRegion.wktFile())
                shapeSave = getattr(vel, 'vr' if doRA else 'vx').shape
                sigShape  = myRegion.sigmaShape()
                sigMask   = makeSigMask(velBase, 1.0, sigShape)
                first = False
                if doVz:
                    muz = np.zeros(mux.shape)

            c1 = getattr(vel, 'vr' if doRA else 'vx')
            if shapeSave != c1.shape:
                # A frame whose mosaic wasn't regenerated to the current
                # velstats-region extent (e.g. its tie/mosaic step was skipped,
                # so it kept a stale size from a prior extent) can't be
                # accumulated into the common-grid arrays. Skip it with a loud,
                # actionable warning instead of aborting the whole track's
                # velStats -- regenerate that frame's velocity mosaic to pick it
                # up. (Was a fatal myerror, which let one stale frame kill the
                # entire track.)
                u.mywarning(f'velocityStats.py: SKIPPING {velFile}: shape '
                            f'{c1.shape} != velstats-region grid {shapeSave} '
                            f'-- stale mosaic, regenerate this frame to include it')
                continue

            iGood = findKeepers(vel, velBase, mean, sigma, sigMask, doRA=doRA)
            sumStats(vel, mux, muy, varx, vary, numAvg, iGood, doRA=doRA)
            if doVz:
                muz[iGood] += vZ.x[iGood]

    if first:
        explainFile = os.path.join(frameDir, 'noVelocityFiles.txt')
        if useSig:
            u.mywarning(f'velocityStats: no velocity files for frames '
                        f'{frame1}-{frame2}; skipping sigma-refinement pass')
            return False
        # Build a default all-no-data output on the SAME polar-stereographic map
        # grid the mosaic velocity for this frame range would have used, read
        # from the tiepoints thumb header's "resolution" line (see
        # _gridFromThumbHeader()). candidateDirs (excluded-but-valid frames with
        # geodats present) still gates whether we emit output at all vs. just
        # warn -- but the grid itself no longer comes from those frames' radar
        # geometry.
        candidateDirs = _findFrameDirsInRange(frame1, frame2)
        grid = _gridFromThumbHeader(frameDir)
        if not candidateDirs or grid is None:
            reason = ('no reference frame geodats' if not candidateDirs
                      else 'no tiepoints thumb-header grid')
            u.mywarning(f'velocityStats: no velocity files and {reason} for frames '
                        f'{frame1}-{frame2} in {frameDir}')
            with open(explainFile, 'w') as fp:
                fp.write(f'No velocity files for frames {frame1}-{frame2}, and '
                         f'{reason} to build a default no-data grid.\n')
            return False
        x0, y0, xs, ys, dxM, dyM = grid
        u.mywarning(f'velocityStats: no velocity files for frames {frame1}-{frame2}; '
                    f'writing all-no-data output on the tiepoints thumb-header PS '
                    f'grid ({xs}x{ys} @ {dxM / 1000.:g} km)')
        domain = 'antarctica' if myRegion.epsg() == 3031 else 'greenland'
        cfg1, cfg2 = ('vr', 'va') if doRA else ('vx', 'vy')
        velBase = u.geoimage(geoType=geoType, verbose=False)
        velBase.geo = u.geodat(x0=x0, y0=y0, xs=xs, ys=ys, dx=dxM, dy=dyM,
                               domain=domain, verbose=False)
        velBase.xyCoordinates()
        shape = (ys, xs)
        mux, muy = np.zeros(shape), np.zeros(shape)
        varx, vary = np.zeros(shape), np.zeros(shape)
        numAvg = np.zeros(shape)
        xps, yps = np.zeros(shape), np.zeros(shape)
        for i in range(len(velBase.yy)):
            xps[i, :] = velBase.xx
        for i in range(len(velBase.xx)):
            yps[:, i] = velBase.yy
        if velMap is not None:
            interp_result = velMap.interpGeo(xps, yps)
            if doRA:
                setattr(velBase, cfg1, np.full(shape, np.nan))
                setattr(velBase, cfg2, np.full(shape, np.nan))
            else:
                setattr(velBase, cfg1, interp_result[0])
                setattr(velBase, cfg2, interp_result[1])
            velBase.v = interp_result[2]
        else:
            setattr(velBase, cfg1, np.full(shape, np.nan))
            setattr(velBase, cfg2, np.full(shape, np.nan))
            velBase.v = np.full(shape, -1.0)
        velBase.v[np.isnan(velBase.v)] = -1.0
        baseOut = frameDir + ('/mosaicBaseRA' if doRA else '/mosaicBase')
        writeMosaicBase(velBase, baseOut, epsg=myRegion.epsg(),
                        wktFile=myRegion.wktFile())
        velBase.v[velBase.v < 0] = 1000000
        sigShape = myRegion.sigmaShape()
        sigMask = makeSigMask(velBase, 1.0, sigShape)
        outGeo = velBase.geo
        velFileSave = candidateDirs[0]
        with open(explainFile, 'w') as fp:
            fp.write(f'No velocity files found for frames {frame1}-{frame2}. '
                     f'All-no-data output written on the tiepoints thumb-header '
                     f'PS grid ({xs}x{ys}, {dxM / 1000.:g} km posting).\n')
        # fall through to normStats/writeStats with numAvg=0 everywhere;
        # computeStats=False because all pixels will be noData (-2e9)
        _noDataWrite = True
    else:
        _noDataWrite = False
        # A prior pass may have left a noVelocityFiles.txt explanation from
        # when this range had no velocity files -- now that it does, remove
        # the stale marker so it can't mislead.
        staleExplain = os.path.join(frameDir, 'noVelocityFiles.txt')
        if os.path.exists(staleExplain):
            os.remove(staleExplain)
    sigx, sigy = normStats(velBase, numAvg, mux, muy, varx, vary, froot,
                           velFileSave, sigMask, doRA=doRA)
    writeStats(mux, muy, sigx, sigy, numAvg, froot, outGeo, doRA=doRA,
               epsg=myRegion.epsg(), wktFile=myRegion.wktFile(),
               computeStats=not _noDataWrite)

    if doVz:
        iGood = numAvg > 1
        muz[iGood] /= numAvg[iGood]
        u.writeImage(froot + '.vz', muz, '>f4')
        copyfile(froot + ('.vr.geodat' if doRA else '.vx.geodat'),
                 froot + '.vz.geodat')


def makeSigMask(vel, sigThresh, sigShape):
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


def getRegion(velFiles):
    ''' get region from "region" file at top level or infer from epsg'''
    if os.path.exists('../../region'):
        with open('../../region', 'r') as fp:
            return fp.readline().strip()
    myPath = os.path.dirname(velFiles[0])
    geodats = u.dols(f'ls {myPath}/geodat*x*.in')
    if len(geodats) > 0:
        geodat = u.geodatrxa(file=geodats[0])
        if geodat.isSouth():
            return 'antarctica'
        else:
            return 'greenland'
    u.myError('getRegion: error parsing region')


def main():
    """ Compute stats for a stack of velocity files """
    noCull, velMapFile, doVz, region, regionFile, doRA, initRef = velocityStatsProcessArgs()

    frameDirs = u.dols("ls -d *-*")
    velFiles = [u.dols('ls -d ../*_*/velocity'),
                u.dols('ls -d ../*_*/velocity_nocull')][noCull]
    if len(velFiles) < 1:
        u.myerror('No files found; in velocityStats directory?')

    if region is None and regionFile is None:
        region = getRegion(velFiles)
    print(f'region {regionFile if regionFile else region}')
    myRegion = s.defaultRegionDefs(region, regionFile=regionFile)

    if doRA:
        print('Mode: RA (vr/va)')
        # velocity.vr/.va/.er/.ea's MEAN is pure data-derived from accumulated
        # real frames -- no basemap blended in (the regional-reference prior
        # for the mean is applied externally, in autoclean's per-frame
        # comparison, via offsets.velocity). The region velocity map is still
        # loaded here because normStats()'s sigma classification bins
        # (iSlow/iMid/iFast/iExtreme) need a real reference speed magnitude
        # (velBase.v) -- see setupFirstVel(), which uses only .v from this
        # map, not vr/va.
        if velMapFile is None:
            velMapFile = myRegion.velMap()
        print(f'velMap {velMapFile}')
        velMap = getFullMap(velMapFile)
    else:
        print('Mode: XY (vx/vy)')
        if velMapFile is None:
            velMapFile = myRegion.velMap()
        print(f'velMap {velMapFile}')
        velMap = getFullMap(velMapFile)

    for frameDir in frameDirs:
        frame1, frame2 = (int(x) for x in frameDir.split('-'))
        print(frame1, frame2)

        if velFileLoop(velFiles, frame1, frame2, frameDir, velMap, noCull,
                       myRegion, doRA=doRA) is not False:
            velFileLoop(velFiles, frame1, frame2, frameDir, velMap, noCull,
                        myRegion, useSig=True, doVz=doVz, doRA=doRA)


if __name__ == '__main__':
    main()
