#!/usr/bin/env python3
"""
autocleanNISAR.py

Per-frame outlier-flagging step for NISAR. Compares the real, measured
range/azimuth offsets for one frame against an expected offset blended from:
  (a) offsets.velocity.vrt -- the region-reference-map simulation already
      computed for this frame by the normal ROFFtoGrimp/SetupNISAR pipeline
      (native radar-pixel-grid, no extra siminsar run needed here), and
  (b) velocityStats' own accumulated, pure-data-derived velocity.vr/.va,
      weighted per-pixel by velocity.navg (more accumulated data -> more
      weight on the data-derived estimate, less on the region reference).

Everything is compared in offset-pixel units (not velocity units). Sigma
comes from velocity.er/.ea alone (unblended) -- velocityStats' normStats()
already produces a sigma that transitions from a region-speed-based floor at
low data counts to a real data-derived sigma at high counts, so no second
sigma source is needed.

Writes range.offsets.good.tif / azimuth.offsets.good.tif (byte, 1=valid,
0=not valid) on the same native pixel grid as
range.offsets.vrt/azimuth.offsets.vrt, then embeds each as a GDAL per-dataset
mask band (GMF_PER_DATASET) directly in the corresponding .vrt's own XML --
this polarity matches GDAL's native MaskBand convention (nonzero=valid), so
no inversion is needed. mosaic3d's readOffsets.c honors this mask
automatically (disable with -noMask).

--allFrames also writes <trackDir>/autocleanSummary.yaml: per-frame percent-
good/bad-pixel-count/pre-baseline-fit-bias stats for range and azimuth, plus
a track-level mean-percent-good aggregate.

Run from a track-N/ directory: autocleanNISAR.py <orbit>_<frame>
"""
import argparse
import glob
import os
import threading
import traceback
import xml.etree.ElementTree as ET
import numpy as np
import yaml
import utilities as u
import sarfunc as s
from osgeo import gdal
from mosaicworkflow.simoffsets import groundToSlantRangeResolution
from mosaicworkflow.velocityStats import referenceSigmaFloor

gdal.UseExceptions()

COUNT_THRESH = 7.0  # matches velocityStats.py normStats()'s countThresh


def resolveRegionDef(trackDir):
    ''' Resolve this project's region definition from <trackDir>/../project.yaml
    (same discovery convention makevelstatsregions.py uses). Needed so the
    u.geoimage().readData(..., epsg=..., wktFile=...) calls in evaluateFrame()
    get the correct polar-stereographic projection -- readData()'s epsg=None
    default silently resolves to Greenland (EPSG:3413) regardless of the
    project, which produces nonsense (all-NaN) interpolation results for any
    other region (e.g. Antarctica, EPSG:3031) if left unset. '''
    projectDir = os.path.dirname(os.path.abspath(trackDir))
    projectYaml = os.path.join(projectDir, 'project.yaml')
    if not os.path.exists(projectYaml):
        u.myerror(f'autocleanNISAR: could not find {projectYaml} to resolve region/epsg')
    with open(projectYaml) as f:
        proj = yaml.safe_load(f) or {}
    regionPath = proj.get('region') or proj.get('regionFile')
    if not regionPath:
        u.myerror(f'autocleanNISAR: {projectYaml} has no region/regionFile key')
    return s.defaultRegionDefs(None, regionFile=regionPath)


def findVelStatsDir(frameNum, framePattern='00??', trackDir=None):
    ''' Return the velocityStats/<f1>-<f2> range-dir name containing frameNum.
    When trackDir is given, scan the ranges that actually exist there (so a
    single track-wide range such as 0000-0099 -- project.yaml
    velocityStatsRegions: track -- is found, preferring one that already holds
    a velocity composite); otherwise, or when no range contains frameNum, fall
    back to the tens-digit grouping of
    setupNISARTracks.py:_vel_stats_dirs_for_track(). '''
    if trackDir is not None:
        hits = []
        for vsd in sorted(glob.glob(os.path.join(trackDir, 'velocityStats', '*-*'))):
            try:
                f1, f2 = (int(x) for x in os.path.basename(vsd).split('-'))
            except ValueError:
                continue
            if f1 <= frameNum <= f2:
                hits.append(vsd)
        if hits:
            withData = [h for h in hits if os.path.exists(os.path.join(h, 'velocity.vr.tif'))]
            return os.path.basename((withData or hits)[0])
    prefix = framePattern.split('?')[0]
    tailLen = framePattern.count('?') - 1
    frameStr = f'{frameNum:0{len(prefix) + framePattern.count("?")}d}'
    mDigit = frameStr[len(prefix)]
    start = prefix + mDigit + '0' * tailLen
    end = prefix + mDigit + '9' * tailLen
    return f'{start}-{end}'


def readOffsetPair(vrtFile, meanBand, sigmaBand):
    ''' Read a range.offsets.vrt/azimuth.offsets.vrt-style file (mean + per-
    pixel measurement sigma bands) into a u.offsets object. '''
    off = u.offsets(vrtFile=vrtFile, verbose=False)
    off.readVrt({meanBand: 'off', sigmaBand: 'sigma'})
    return off


def readIonCorrection(frameDir, rangeOff):
    ''' Return the ionosphere range correction array (already in the same
    offset-pixel units as range.offsets.vrt's own RangeOffsets band -- see
    estimateIonosphere.py's ".offset" variant), or None if range.offsets.vrt
    has no ionosphereRangeOffsetCorrection metadata tag (matches
    makeMaster.py:find_ion_correction()'s resolution convention). Ionosphere
    is a range/dispersive effect only -- no azimuth equivalent. '''
    ionName = rangeOff.meta.get('ionosphereRangeOffsetCorrection')
    if not ionName:
        return None
    ionVrt = os.path.join(frameDir, ionName)
    if not os.path.exists(ionVrt):
        return None
    ion = u.offsets(vrtFile=ionVrt, verbose=False)
    ion.readVrt({'ionosphereCorrection': 'ionCorrection'})
    corr = ion.ionCorrection
    if corr.shape == rangeOff.off.shape:
        return corr
    # The merged ionosphere product can be SHORTER than the merged offsets when
    # a constituent frame has no ionosphere estimate at all (a 10-second
    # dual-bandwidth partial frame whose chosen granule carries no unwrapped
    # phase: Antarctica80 track-119/3678_123, track-120/3679_126, track-104/
    # 4701_123). Both rasters share the radar grid and its geotransform, so
    # place the correction by its azimuth-time offset and leave the uncovered
    # rows at zero (no correction) rather than dying on the shape mismatch and
    # leaving the whole virtual frame without an autoclean mask or a tie.
    try:
        gtR = gdal.Open(rangeOff.vrtFile).GetGeoTransform()
        gtI = gdal.Open(ionVrt).GetGeoTransform()
        colOff = int(round((gtI[0] - gtR[0]) / gtR[1]))
        rowOff = int(round((gtI[3] - gtR[3]) / gtR[5]))
    except Exception as e:
        u.mywarning(f'readIonCorrection: {ionVrt} shape {corr.shape} != offsets '
                    f'{rangeOff.off.shape} and no usable geotransform ({e}); '
                    f'ignoring the ionosphere correction for this frame')
        return None
    full = np.zeros(rangeOff.off.shape, dtype=corr.dtype)
    r0, c0 = max(rowOff, 0), max(colOff, 0)
    r1 = min(rowOff + corr.shape[0], full.shape[0])
    c1 = min(colOff + corr.shape[1], full.shape[1])
    if r1 <= r0 or c1 <= c0:
        u.mywarning(f'readIonCorrection: {ionVrt} does not overlap the offsets grid; '
                    f'ignoring the ionosphere correction for this frame')
        return None
    full[r0:r1, c0:c1] = corr[r0 - rowOff:r1 - rowOff, c0 - colOff:c1 - colOff]
    u.mywarning(f'readIonCorrection: {os.path.basename(ionVrt)} covers rows {r0}-{r1 - 1} '
                f'of {full.shape[0]} (a constituent frame has no ionosphere estimate); '
                f'uncovered rows get zero correction')
    return full


def writeGoodMask(templateVrt, goodArray, outTif):
    ''' Write a byte GeoTIFF (1=valid, 0=not valid) on the same pixel grid
    as templateVrt. '''
    src = gdal.Open(templateVrt)
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(outTif, src.RasterXSize, src.RasterYSize, 1,
                       gdal.GDT_Byte, options=['COMPRESS=DEFLATE'])
    ds.SetGeoTransform(src.GetGeoTransform())
    ds.SetProjection(src.GetProjection())
    ds.GetRasterBand(1).WriteArray(goodArray.astype(np.uint8))
    ds.FlushCache()
    ds, src = None, None


def embedMaskBand(vrtFile, maskTif):
    ''' Embed maskTif (byte GeoTIFF, 1=valid/0=not, same grid as vrtFile) as a
    GDAL per-dataset mask band (GMF_PER_DATASET) directly in vrtFile's own
    XML, so mosaic3d's readOffsets.c (GDALGetMaskBand/GDALGetMaskFlags) picks
    it up with no separate embedding tool. Idempotent -- replaces any
    existing <MaskBand> block from a prior autoclean run. '''
    tree = ET.parse(vrtFile)
    root = tree.getroot()
    for existing in root.findall('MaskBand'):
        root.remove(existing)
    maskBand = ET.SubElement(root, 'MaskBand')
    maskRasterBand = ET.SubElement(maskBand, 'VRTRasterBand', {'dataType': 'Byte'})
    simpleSource = ET.SubElement(maskRasterBand, 'SimpleSource')
    sourceFilename = ET.SubElement(simpleSource, 'SourceFilename', {'relativeToVRT': '1'})
    sourceFilename.text = os.path.basename(maskTif)
    ET.SubElement(simpleSource, 'SourceBand').text = '1'
    tree.write(vrtFile)


def evaluateFrame(frameDir, trackDir='.', sigThresh=3.0, framePattern='00??'):
    '''
    Compute range/azimuth valid-pixel masks for one real frame directory
    (e.g. track-N/<orbit>_<frame>). Returns (goodR, goodA, biasR, biasA) --
    boolean arrays plus the median pre-baseline-fit bias (px) removed from
    each component -- and writes range.offsets.good.tif/azimuth.offsets.good.tif
    in frameDir.
    '''
    orbit, frameNum = (int(x) for x in os.path.basename(frameDir).split('_'))
    velStatsDir = os.path.join(trackDir, 'velocityStats',
                               findVelStatsDir(frameNum, framePattern, trackDir))

    rangeVrt = os.path.join(frameDir, 'range.offsets.vrt')
    azVrt = os.path.join(frameDir, 'azimuth.offsets.vrt')
    velVrt = os.path.join(frameDir, 'offsets.velocity.vrt')
    geomVrt = os.path.join(frameDir, 'offsets.geom.vrt')
    llVrt = os.path.join(frameDir, 'offsets.velocity.ll.vrt')

    rangeOff = readOffsetPair(rangeVrt, 'RangeOffsets', 'RangeSigma')
    azOff = readOffsetPair(azVrt, 'AzimuthOffsets', 'AzimuthSigma')

    geodatFile = os.path.join(frameDir, rangeOff.meta['geo1'])
    rangeOff.geodatrxa = u.geodatrxa(file=geodatFile, echo=False)
    rangeOff.slpRg, rangeOff.slpAz = rangeOff.geodatrxa.singleLookResolution()

    # range.offsets.vrt's own RangeOffsets band is not ionosphere-corrected --
    # add the correction (already in the same offset-pixel units) before
    # comparing against the expected offset, or the ionosphere signal itself
    # would show up as a spurious range-only outlier.
    ionCorrection = readIonCorrection(frameDir, rangeOff)
    rangeMeasured = rangeOff.off if ionCorrection is None \
        else rangeOff.off + ionCorrection

    ll = u.offsets(vrtFile=llVrt, verbose=False)
    ll.readVrt({'lat': 'lat', 'lon': 'lon'})

    # target attribute names avoid readVrt()'s special {rgOff,azOff,SigmaR,
    # SigmaA} set -- those trigger a file-path-recording lookup (r.files[i])
    # that indexes out of range for a pixel-interleaved single-tif VRT like
    # this one (one .tif backs both bands, not one .tif per band).
    velSim = u.offsets(vrtFile=velVrt, verbose=False)
    velSim.readVrt({'RangeOffsets': 'rgOffSim', 'AzimuthOffsets': 'azOffSim'})
    deltaT = float(velSim.meta['deltaT'])

    # offsets.velocity = geom (zero-velocity static baseline) + velocity*scale
    # (mirrors simoffsets.py:main()'s own "dr = dr0 + vr*scaleR" construction).
    # The velStats-driven term below must add this same static baseline back
    # in -- omitting it means comparing a measurement dominated by a large
    # geometric offset (often tens to hundreds of pixels) against a pure
    # velocity-only correction, which is wrong by exactly that geometric term.
    geom = u.offsets(vrtFile=geomVrt, verbose=False)
    geom.readVrt({'RangeOffsets': 'rgGeom', 'AzimuthOffsets': 'azGeom'})

    # scale factors: range varies with incidence angle across the swath;
    # azimuth does not (no incidence-angle dependence).
    groundToSlant = groundToSlantRangeResolution(rangeOff)
    scaleR = groundToSlant * deltaT / 365.
    scaleA = deltaT / 365. / rangeOff.slpAz

    # velocityStats reference: pure data-derived mean/sigma/navg. epsg/wktFile
    # must be passed explicitly -- readData()'s epsg=None default resolves to
    # Greenland regardless of the actual project (see resolveRegionDef()).
    regionDef = resolveRegionDef(trackDir)
    epsg, wktFile = regionDef.epsg(), regionDef.wktFile()
    statsBase = os.path.join(velStatsDir, 'velocity')
    # Read only the window of the composite this frame touches (a track-wide
    # composite can be hundreds of Mpx; the frame needs a few percent of it),
    # and tolerate a missing composite (no products yet in this range) by
    # falling back to the reference-only blend below (wData = 0 everywhere)
    # instead of dying in readgeodat -- that failure used to drop every frame
    # of a new range from the master until someone rebuilt the stats.
    vr = va = er = ea = nAvgInterp = None
    if os.path.exists(statsBase + '.vr.tif'):
        probe = u.geoimage(geoType='velocityRA', verbose=False)
        probe.getGeoFile(statsBase, probe.getDomain(epsg), tiff=True,
                         wkt=probe.getWKT_PROJ(epsg, wktFile))
        xps, yps = probe.geo.lltoxykm(ll.lat, ll.lon)
        box = (np.nanmin(xps), np.nanmax(xps), np.nanmin(yps), np.nanmax(yps))
        mean = u.geoimage(geoType='velocityRA', verbose=False)
        sigma = u.geoimage(geoType='errorRA', verbose=False)
        navg = u.geoimage(geoType='scalar', verbose=False)
        if (mean.readDataWindow(statsBase, *box, tiff=True, epsg=epsg, wktFile=wktFile) is not None
                and sigma.readDataWindow(statsBase, *box, tiff=True, epsg=epsg, wktFile=wktFile) is not None
                and navg.readDataWindow(statsBase + '.navg', *box, tiff=True, epsg=epsg, wktFile=wktFile) is not None):
            mean.setupInterp()
            sigma.setupInterp()
            navg.setupInterp()
            xps, yps = mean.geo.lltoxykm(ll.lat, ll.lon)
            vr, va, _ = mean.interpGeo(xps, yps)
            er, ea, _ = sigma.interpGeo(xps, yps)
            nAvgInterp = navg.interpGeo(xps, yps)
        else:
            u.mywarning(f'autocleanNISAR: {frameDir} does not overlap the velocityStats '
                        f'composite in {velStatsDir} -- reference-only mask')
    else:
        u.mywarning(f'autocleanNISAR: no velocity composite in {velStatsDir} '
                    f'(velocityStats not run yet for this range?) -- reference-only mask')
    if vr is None:
        shape = ll.lat.shape
        vr, va, er, ea = (np.full(shape, np.nan) for _ in range(4))
        nAvgInterp = np.zeros(shape)

    nAvgSafe = np.nan_to_num(nAvgInterp, nan=0.)
    vrSafe = np.nan_to_num(vr, nan=0.)
    vaSafe = np.nan_to_num(va, nan=0.)
    wData = np.clip(nAvgSafe, 0, COUNT_THRESH) / COUNT_THRESH
    wRef = 1. - wData

    drExpected = wRef * velSim.rgOffSim + wData * (geom.rgGeom + vrSafe * scaleR)
    daExpected = wRef * velSim.azOffSim + wData * (geom.azGeom + vaSafe * scaleA)

    # Component-aware reference sigma floor. velocityStats' map-grid sigma
    # (er/ea) is isotropic in RA mode -- with no vr/va basemap on the map grid,
    # normStats()'s component-linear term (15 + 0.2*|v_comp|) never fires and
    # only the speed-bin floor survives. Recover the range/azimuth reference
    # velocity components from the synthetic offsets (invert dr = geom + v*scale,
    # the construction offsets.velocity uses) and rebuild the SAME piecewise
    # normStats floor per component via referenceSigmaFloor() -- so this matches
    # the long-standing XY autoclean clamp strategy, modulo vr/va vs vx/vy. Blend
    # with the data-derived er/ea using the same navg weight as drExpected: pure
    # reference floor where no accumulated data (wData=0), pure velStats sigma
    # where plentiful. er/ea are NaN-safed since they carry zero weight (wData=0)
    # exactly where they are NaN (out of the velStats grid), and 0*NaN != 0.
    vrRef = (velSim.rgOffSim - geom.rgGeom) / scaleR
    vaRef = (velSim.azOffSim - geom.azGeom) / scaleA
    speedRef = np.sqrt(vrRef ** 2 + vaRef ** 2)
    sigRefR = referenceSigmaFloor(vrRef, speedRef)
    sigRefA = referenceSigmaFloor(vaRef, speedRef)
    erSafe = np.nan_to_num(er, nan=0.)
    eaSafe = np.nan_to_num(ea, nan=0.)
    sigmaDr = (wRef * sigRefR + wData * erSafe) * scaleR
    sigmaDa = (wRef * sigRefA + wData * eaSafe) * scaleA

    # This frame's offsets have not yet had a baseline fit (rparams/azparams)
    # applied, so drExpected/daExpected -- built from velStats/velSim, both of
    # which are baseline-independent -- can differ from the raw measurement by
    # a systematic, roughly frame-constant bias (the uncorrected orbit
    # baseline) on top of any genuine per-pixel outliers. Remove that shared
    # bias (per component, via the median -- robust even when a large
    # fraction of the frame is currently reading as "bad", unlike a plain
    # mean, which that same population would drag off-center) before testing
    # against sigma, or the sigma test conflates "biased" with "bad".
    diffR = rangeMeasured - drExpected
    diffA = azOff.off - daExpected
    biasR = np.nanmedian(diffR)
    biasA = np.nanmedian(diffA)
    print(f'autocleanNISAR: median range/azimuth bias (pre-baseline-fit) = '
          f'{biasR:.3f} / {biasA:.3f} px')

    # NaN residuals (measurement gaps, or velStats/interp gaps) compare False
    # in both directions, so they correctly fall out as "not good" here
    # rather than needing explicit NaN handling.
    goodR = np.abs(diffR - biasR) <= sigThresh * sigmaDr
    goodA = np.abs(diffA - biasA) <= sigThresh * sigmaDa

    rangeMaskTif = os.path.join(frameDir, 'range.offsets.good.tif')
    azMaskTif = os.path.join(frameDir, 'azimuth.offsets.good.tif')
    writeGoodMask(rangeVrt, goodR, rangeMaskTif)
    writeGoodMask(azVrt, goodA, azMaskTif)
    embedMaskBand(rangeVrt, rangeMaskTif)
    embedMaskBand(azVrt, azMaskTif)
    return goodR, goodA, float(biasR), float(biasA)


def findFrameDirs(trackDir, framePattern='00??'):
    ''' Return merged/virtual frame dirs in trackDir matching *_<framePattern>
    -- the same glob setupNISARTracks.py uses for merged frames (framePattern's
    own '?' chars serve directly as glob wildcards, e.g. '00??' -> '*_00??'). '''
    return sorted(glob.glob(os.path.join(trackDir, f'*_{framePattern}')))


def frameHasMask(frameDir):
    ''' Return True if frameDir/range.offsets.vrt already carries an embedded
    autoclean mask (<MaskBand> block written by embedMaskBand). Mirrors
    setupNISARTracks.py:_track_has_autoclean_mask() at the single-frame
    granularity -- used by the -new filter to skip already-cleaned frames.
    False if the VRT is missing or has no <MaskBand>. '''
    vrt = os.path.join(frameDir, 'range.offsets.vrt')
    if not os.path.exists(vrt):
        return False
    with open(vrt) as f:
        return '<MaskBand' in f.read()


def filterNewFrames(frameDirs, new):
    ''' When new is True, drop frame dirs that already have an embedded
    autoclean mask, so only not-yet-cleaned frames are (re)processed. When
    False, return frameDirs unchanged. '''
    if not new:
        return frameDirs
    return [d for d in frameDirs if not frameHasMask(d)]


def _evaluateFrameThread(frameDir, trackDir, sigThresh, framePattern, results):
    ''' Thread target: run evaluateFrame for one frame, recording either
    (nBadR, nBadA, total, biasR, biasA) or the caught exception into
    results[frameDir]. Range and azimuth are evaluated (and reported)
    completely independently, same as evaluateFrame itself -- unlike
    autoclean.py's XY mode, a bad range pixel does not imply or get combined
    with a bad azimuth pixel. Catches SystemExit too, not just Exception --
    utilities.myerror() (used throughout u.offsets/u.geoimage for
    missing-file errors) calls sys.exit(), which is a BaseException, not an
    Exception, so it would otherwise silently kill this thread and leave
    results[frameDir] unset. '''
    try:
        goodR, goodA, biasR, biasA = evaluateFrame(frameDir, trackDir=trackDir,
                                                    sigThresh=sigThresh,
                                                    framePattern=framePattern)
        results[frameDir] = (goodR.size - int(np.sum(goodR)),
                             goodA.size - int(np.sum(goodA)), goodR.size,
                             biasR, biasA)
    except (Exception, SystemExit) as e:
        # str(SystemExit(1)) is just "1": u.myerror() prints the real message to stdout
        # and exits, so the reason is lost by the time it reaches the summary. stdout
        # cannot be captured per frame here -- redirect_stdout rebinds the process-wide
        # sys.stdout, which would corrupt every other thread in the pool -- but the
        # traceback is thread-local, so record where it came from instead. Enough to
        # tell a missing-input failure from a genuine computation error.
        stack = traceback.extract_tb(e.__traceback__)
        e._acDetail = f'{type(e).__name__}: {e}   at ' + ' <- '.join(
            f'{os.path.basename(fr.filename)}:{fr.lineno} {fr.line}'
            for fr in reversed(stack[-3:]) if fr.line)
        results[frameDir] = e


def writeTrackSummary(trackDir, frameDirs, results, sigThresh):
    '''
    Print the per-frame range/azimuth bad-pixel summary for one track and
    write <trackDir>/autocleanSummary.yaml (per-frame percent-good/bias/
    pixel-count stats plus a track-level aggregate). Factored out of
    evaluateTrack so evaluateAllTracks can share it.
    '''
    frameSummaries = {}
    nFailed = 0
    for frameDir in frameDirs:
        label = os.path.basename(frameDir)
        result = results.get(frameDir)
        if not isinstance(result, tuple):
            detail = getattr(result, '_acDetail', None) or str(result)
            print(f'{label}: FAILED ({detail})')
            frameSummaries[label] = {'failed': str(detail)}
            nFailed += 1
        else:
            nBadR, nBadA, total, biasR, biasA = result
            print(f'{label}: {nBadR} bad range px, {nBadA} bad '
                  f'azimuth px (of {total})')
            frameSummaries[label] = {
                'rangeGoodPercent': round(100. * (total - nBadR) / total, 2),
                'azimuthGoodPercent': round(100. * (total - nBadA) / total, 2),
                'rangeBadPixels': nBadR,
                'azimuthBadPixels': nBadA,
                'totalPixels': total,
                'rangeBiasPx': round(biasR, 4),
                'azimuthBiasPx': round(biasA, 4),
            }
    if nFailed > 0:
        u.mywarning(f'{nFailed} of {len(frameDirs)} frames failed')

    ok = [v for v in frameSummaries.values() if 'failed' not in v]
    summaryOut = {
        'sigThresh': sigThresh,
        'frames': frameSummaries,
        'summary': {
            'nFrames': len(frameDirs),
            'nFailed': nFailed,
            'meanRangeGoodPercent':
                round(float(np.mean([v['rangeGoodPercent'] for v in ok])), 2) if ok else None,
            'meanAzimuthGoodPercent':
                round(float(np.mean([v['azimuthGoodPercent'] for v in ok])), 2) if ok else None,
        },
    }
    summaryFile = os.path.join(trackDir, 'autocleanSummary.yaml')
    with open(summaryFile, 'w') as fp:
        yaml.safe_dump(summaryOut, fp, sort_keys=True, default_flow_style=False)
    print(f'Wrote {summaryFile}')


def evaluateTrack(trackDir='.', sigThresh=3.0, framePattern='00??', nThreads=8,
                  new=False):
    '''
    Run evaluateFrame over every merged frame directory in trackDir, in
    parallel (modeled on autoclean.py:main()'s threaded loop over
    getOffsetDirs()). Prints a per-frame range/azimuth bad-pixel summary
    (kept separate, never combined), writes <trackDir>/autocleanSummary.yaml
    (per-frame percent-good/bias/pixel-count stats plus a track-level
    aggregate), and returns {frameDir: (nBadR, nBadA, total, biasR, biasA)
    or Exception}. When new is True, frames that already carry an embedded
    autoclean mask are skipped -- only not-yet-cleaned frames are processed.
    '''
    frameDirs = filterNewFrames(findFrameDirs(trackDir, framePattern), new)
    if len(frameDirs) < 1:
        u.mywarning(f'No *_{framePattern} frame dirs found in {trackDir}'
                    + (' needing autoclean (all already masked)' if new else ''))
        return {}

    results = {}
    threads = [threading.Thread(target=_evaluateFrameThread,
                                args=[frameDir, trackDir, sigThresh,
                                      framePattern, results])
              for frameDir in frameDirs]
    u.runMyThreads(threads, nThreads, 'autocleanNISAR')

    writeTrackSummary(trackDir, frameDirs, results, sigThresh)
    return results


def evaluateAllTracks(trackDirs, sigThresh=3.0, framePattern='00??', nThreads=8,
                      new=False):
    '''
    Run evaluateFrame over every merged frame directory of EVERY track in
    trackDirs as one flat multithreaded pool (no per-track barrier -- unlike
    looping evaluateTrack per track, a straggler frame in one track does not
    hold up work in the next). Run from the project root with track-relative
    (or absolute) track dirs. Writes each track's autocleanSummary.yaml, and
    returns the combined {frameDir: result-or-Exception} dict. A single shared
    results dict is safe: keys are frameDir paths, unique across tracks
    because each is prefixed by its own trackDir. When new is True, frames
    that already carry an embedded autoclean mask are skipped -- only
    not-yet-cleaned frames are processed.
    '''
    results = {}
    allThreads = []
    trackFrames = {}
    summarised = set()
    summaryLock = threading.Lock()

    def frameThenMaybeSummary(frameDir, trackDir, sigThresh, framePattern, results):
        '''Run one frame, then write its track's summary the moment that frame is the
        last of the track to finish. Writing per track as it completes -- rather than
        only after the whole pool joins -- means a run that is killed part-way still
        leaves the diagnostics for every track it got through. Without this, an
        interrupted run loses the reason every frame failed, which is exactly what
        happened nightly on Antarctica40. The pool itself keeps its no-barrier
        behaviour: nothing waits for a track, the summary is just written by whichever
        thread happens to close it out.'''
        try:
            _evaluateFrameThread(frameDir, trackDir, sigThresh, framePattern, results)
        finally:
            with summaryLock:
                frames = trackFrames.get(trackDir, [])
                if trackDir in summarised or any(f not in results for f in frames):
                    return
                summarised.add(trackDir)
            print(f'--- {trackDir} ---')
            writeTrackSummary(trackDir, frames, results, sigThresh)

    for trackDir in trackDirs:
        allFrameDirs = findFrameDirs(trackDir, framePattern)
        if len(allFrameDirs) < 1:
            u.mywarning(f'No *_{framePattern} frame dirs found in {trackDir}')
            continue
        frameDirs = filterNewFrames(allFrameDirs, new)
        if len(frameDirs) < 1:
            continue
        trackFrames[trackDir] = frameDirs
        allThreads += [threading.Thread(target=frameThenMaybeSummary,
                                        args=[frameDir, trackDir, sigThresh,
                                              framePattern, results])
                       for frameDir in frameDirs]
    if not allThreads:
        u.mywarning('evaluateAllTracks: no frames found in any track')
        return {}
    print(f'autocleanNISAR: {len(allThreads)} frame jobs across '
          f'{len(trackFrames)} tracks, {nThreads} threads')
    u.runMyThreads(allThreads, nThreads, 'autocleanNISAR')

    # Anything not already summarised above -- a track whose last frame died in a way
    # that skipped the finally block, so its summary was never triggered.
    for trackDir, frameDirs in trackFrames.items():
        if trackDir in summarised:
            continue
        print(f'--- {trackDir} ---')
        writeTrackSummary(trackDir, frameDirs, results, sigThresh)
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Flag outlier range/azimuth offsets for one real NISAR '
                    'frame by comparison with offsets.velocity + velocityStats. '
                    'Run from the track-N/ directory.',
        epilog='Part of the mosaicworkflow package.')
    parser.add_argument('frame', metavar='orbit_frame', nargs='?', default=None,
                        help='Frame directory to evaluate, e.g. 1830_0000 '
                            '(omit when using -allFrames)')
    parser.add_argument('-sigThresh', '--sigThresh', type=float, default=3.0,
                        help='Discard threshold in units of sigma [3.0]')
    parser.add_argument('-framePattern', '--framePattern', default='00??',
                        help='Virtual-frame suffix pattern [00??]')
    parser.add_argument('-allFrames', '--allFrames', action='store_true',
                        help='Process every merged frame dir in the current '
                            'track-N/ directory instead of a single frame')
    parser.add_argument('-tracks', '--tracks', nargs='*', default=None,
                        metavar='TRACKDIR',
                        help='Run from the project root: process every merged '
                            'frame of the listed track dirs (default with no '
                            'values: all track-* dirs in the cwd) as ONE '
                            'multithreaded pool -- no per-track barrier, '
                            'unlike running -allFrames track by track')
    parser.add_argument('-threads', '--threads', type=int, default=8,
                        help='Number of parallel threads for '
                            '-allFrames/-tracks [8]')
    parser.add_argument('-new', '--new', action='store_true',
                        help='With -allFrames/-tracks, skip frames whose '
                            'range.offsets.vrt already carries an embedded '
                            'autoclean mask; only (re)clean not-yet-cleaned '
                            'frames. Ignored for a single explicit frame.')
    args = parser.parse_args()

    if args.tracks is not None:
        if args.frame is not None or args.allFrames:
            parser.error('-tracks cannot be combined with a frame '
                         'argument or -allFrames')
        trackDirs = args.tracks or sorted(
            d for d in glob.glob('track-*') if os.path.isdir(d))
        if not trackDirs:
            parser.error('-tracks: no track-* dirs found in the current '
                         'directory (run from the project root)')
        evaluateAllTracks(trackDirs, sigThresh=args.sigThresh,
                          framePattern=args.framePattern,
                          nThreads=args.threads, new=args.new)
        return

    if args.allFrames:
        evaluateTrack(trackDir='.', sigThresh=args.sigThresh,
                     framePattern=args.framePattern, nThreads=args.threads,
                     new=args.new)
        return

    if args.frame is None:
        parser.error('frame is required unless -allFrames is given')

    goodR, goodA, biasR, biasA = evaluateFrame(args.frame, trackDir='.',
                                               sigThresh=args.sigThresh,
                                               framePattern=args.framePattern)
    nBadR = goodR.size - np.sum(goodR)
    nBadA = goodA.size - np.sum(goodA)
    print(f'{args.frame}: {nBadR} bad range px, {nBadA} bad '
          f'azimuth px (of {goodR.size})')


if __name__ == '__main__':
    main()
