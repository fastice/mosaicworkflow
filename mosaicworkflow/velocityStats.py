#!/usr/bin/env python3
import argparse
import utilities as u
import numpy as np
import os
import shapefile
import sarfunc as s

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
    return args.nocull, args.velmap, args.sumVz, args.region, doRA, args.initializeReference


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
    velMap.readData(velMapFile)
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
    return orbit, frame


def setupFirstVel(vel, velFile, frameDir, velMap, doRA=False, tiffMode=False):
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
        setattr(velBase, cfg1, interp_result[0])
        setattr(velBase, cfg2, interp_result[1])
        velBase.v = interp_result[2]
    else:
        setattr(velBase, cfg1, np.full(c1.shape, np.nan))
        setattr(velBase, cfg2, np.full(c1.shape, np.nan))
        velBase.v = np.full(c1.shape, -1.0)

    velBase.v[np.isnan(velBase.v)] = -1.0

    baseOut = frameDir + ('/mosaicBaseRA' if doRA else '/mosaicBase')
    velBase.writeData(baseOut)

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


def writeStats(mux, muy, sigx, sigy, numAvg, froot, geo, doRA=False):
    ''' write stats results '''
    s1, s2, e1, e2 = ('.vr', '.va', '.er', '.ea') if doRA \
        else ('.vx', '.vy', '.ex', '.ey')
    u.writeImage(froot + s1, mux,   '>f4')
    u.writeImage(froot + s2, muy,   '>f4')
    u.writeImage(froot + e1, sigx,  '>f4')
    u.writeImage(froot + e2, sigy,  '>f4')
    u.writeImage(froot + '.navg', numAvg, '>f4')
    # write geodat sidecars using the geo object (works for both tiff and binary)
    for sfx in (s1, s2, e1, e2, '.navg'):
        geo.writeGeodat(froot + sfx + '.geodat')


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
        sigma = u.geoimage(geoType='error', verbose=False)
        mean  = u.geoimage(geoType=geoType, verbose=False)
        # load sigma components — stored under ex/ey or er/ea depending on mode
        sigma.readData(froot, geoType='error')
        # override attribute names if RA so findKeepers can access .er/.ea
        if doRA:
            sigma.er = sigma.ex
            sigma.ea = sigma.ey
        mean.readData(froot)
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
                                  doRA=doRA, tiffMode=tiffMode)
                shapeSave = getattr(vel, 'vr' if doRA else 'vx').shape
                sigShape  = myRegion.sigmaShape()
                sigMask   = makeSigMask(velBase, 1.0, sigShape)
                first = False
                if doVz:
                    muz = np.zeros(mux.shape)

            c1 = getattr(vel, 'vr' if doRA else 'vx')
            if shapeSave != c1.shape:
                u.myerror(f'velocityStats.py: shape {c1.shape} for '
                          f'{velFile} different from prior: {shapeSave}')

            iGood = findKeepers(vel, velBase, mean, sigma, sigMask, doRA=doRA)
            sumStats(vel, mux, muy, varx, vary, numAvg, iGood, doRA=doRA)
            if doVz:
                muz[iGood] += vZ.x[iGood]

    sigx, sigy = normStats(velBase, numAvg, mux, muy, varx, vary, froot,
                           velFileSave, sigMask, doRA=doRA)
    writeStats(mux, muy, sigx, sigy, numAvg, froot, outGeo, doRA=doRA)

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
    noCull, velMapFile, doVz, region, doRA, initRef = velocityStatsProcessArgs()

    frameDirs = u.dols("ls -d *-*")
    velFiles = [u.dols('ls -d ../*_*/velocity'),
                u.dols('ls -d ../*_*/velocity_nocull')][noCull]
    if len(velFiles) < 1:
        u.myerror('No files found; in velocityStats directory?')

    if region is None:
        region = getRegion(velFiles)
    print(f'region {region}')
    myRegion = s.defaultRegionDefs(region)

    if doRA:
        print('Mode: RA (vr/va)')
        velMap = None  # basemap loaded per-frameDir from sims/
    else:
        print('Mode: XY (vx/vy)')
        if velMapFile is None:
            velMapFile = myRegion.velMap()
        print(f'velMap {velMapFile}')
        velMap = getFullMap(velMapFile)

    for frameDir in frameDirs:
        frame1, frame2 = (int(x) for x in frameDir.split('-'))
        print(frame1, frame2)

        if doRA:
            vrPath = os.path.join(frameDir, 'sims', 'mosaicOffsets.vr')
            if initRef or not os.path.exists(vrPath):
                from mosaicworkflow.initRAReference import initRAReference
                myDem = myRegion.dem()
                initRAReference(frameDir, trackDir='../..', region=region,
                                dem=myDem)
            if not os.path.exists(vrPath):
                u.myerror(f'RA basemap missing: {vrPath}  '
                          f'(run with --initializeReference)')
            velMap = getRABasemap(frameDir)

        velFileLoop(velFiles, frame1, frame2, frameDir, velMap, noCull,
                    myRegion, doRA=doRA)
        velFileLoop(velFiles, frame1, frame2, frameDir, velMap, noCull,
                    myRegion, useSig=True, doVz=doVz, doRA=doRA)


if __name__ == '__main__':
    main()
