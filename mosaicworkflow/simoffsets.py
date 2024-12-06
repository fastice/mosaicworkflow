#!/usr/bin/env python3
import utilities as u
import numpy as np
import sys
import os
from subprocess import call
from utilities import geodatrxa
import sarfunc as s
import argparse
import glob


def simoffsetsUsage():
    region1 = s.defaultRegionDefs('greenland')
    region2 = s.defaultRegionDefs('antarctica')
    print('\033[1m\n Use a velocity map to simulate offsets for initial guess'
          ' in feature tracking\n')
    print('\tsimoffsets -help ')
    print('\t\t-azOffsets = file[offsets.da] ')
    print('\t\t-offsetsDat = dat or vrt file[offsets.dat]')
    print('\t\t-syncDat 1) Force siminsar to be run/rerun 2) use offsetsDat '
          'to size result [False]')
    print('\t\t-noVel Compute offsets for assuming 0 vel everywhere [False]')
    print('\t\t-fastMask =  Use the fast mask for mask generation [False]')
    print('\t\t-velMap = filebase[default for region]')
    print('\t\t-dem = file[default for region] ')
    print('\t\t-region = region name [greenland/antarctica for NH/SH]')
    print('\t\t-secondDir = dir[../secondslcdir] ')

    print('\t\t-geodatFile = file = [geodat10x2.in] ')
    print('\t\t-secondGeodatFile = dir[secondDir/geodatFile] ')
    print(f'\t\t-maskInputFile  file = [{region1.mask()} or {region2.mask()}]')
    print('\tNote:\n\t\t1) Run in directory for the desired image pair ')
    print('\t\t2) Default region is greenland for NH geodat and antarctica '
          'for SH geodat')
    print('\t\t3) produces the following: offsets.da(.da.dat, .dr, .dr.dat,'
          ' .dat, .lat, .lon, .mask, .simdat) \033[0m\n')
    exit()


def resolveDefault(args, key, default):
    if getattr(args, key) is not None:
        return getattr(args, key)
    return default

def simOffsetsProcessArgs1(fp):
    """
    Process arguments and defaults
    """

    # Resolve slcs other wise default to current dir
    slcs = glob.glob('ls -l *.slc')
    second = [s for s in slcs if '../' in s]
    if len(second) == 1:
        secondDirDefault = '../'+second[0].split('/')[1]
    else:
        secondDirDefault = '.'

    epilog = 'Notes: 1) Run in directory for the desired image pair ' \
        '2) Default region is greenland for NH or SH as determined from ' \
        'geodat 3) produces the following: offsets.da(.da.dat, .dr, ' \
        '.dr.dat, .dat, .lat, .lon, .mask, .simdat) '
    parser = argparse.ArgumentParser(description='\033[1m\n Use a velocity map'
                                     ' to simulate offsets for initial guess'
                                     ' in feature tracking\033[0m\n',
                                     epilog=epilog)

    parser.add_argument('-azOffsets', '--azOffsets', type=str,
                        default='offsets.da',
                        help='azOffsets file [offsets.da]')
    parser.add_argument('-offsetsDat', '--offsetsDat', type=str,
                        default='offsets.dat',
                        help='dat or vrt file to define output [offsets.dat]')
    parser.add_argument('-syncDat', '--syncDat',  action='store_true',
                        default=False, help='Force siminsar to be run '
                        '2) use offsetsDat to size result')
    parser.add_argument('-noVel', '--noVel',  action='store_true',
                        default=False,
                        help='Compute offsets for assuming 0 vel everywhere')
    parser.add_argument('-fastMask', '--fastMask',  action='store_true',
                        default=False,
                        help='Use the fast mask for mask generation')
    parser.add_argument('-LSB ', '--LSB',  action='store_true',
                        default=False,
                        help='Use the fast mask for mask generation')
    parser.add_argument('-velMap', '--velMap', type=str,
                        default=None, help='Root name for velocity map '
                        '[from Region file or default regionDefs]')
    parser.add_argument('-dem', '--dem', type=str, nargs=1, default=None,
                        help='DEM [from Region file or default regionDefs]')

    parser.add_argument('-region', '--region', type=str, default=None,
                        help='Predefined region if on GrIMP server'
                        ' [greenland/antarctica for NH/SH]]')
    parser.add_argument('-regionFile', '--regionFile', type=str,
                        default=None, help='Yaml file with locations of '
                        'velMap, DEM etc [None]')
    parser.add_argument('-secondDir', '--secondDir', type=str,
                        default=secondDirDefault,
                        help='Location of second directory for '
                        'the pair [../secondslcdir]')
    parser.add_argument('-geodatFile', '--geodatFile', type=str,
                        default='geodat10x2.in',
                        help='GrIMP geodat.in or geodat.geojson'
                        'file which defines the image geometry '
                        '[geoadat2x10.in]')
    parser.add_argument('-secondGeodatFile', '--secondGeodatFile', type=str,
                        default=None, help='GrIMP geodat.in or '
                        'geodat.geojson file for the second image in the pair'
                        '[secondir/geodat2x10.in]')
    parser.add_argument('-maskInputFile', '--maskInputFile', type=str,
                        default=None, help='Flat byte file with mask'
                        '[None]')

    #
    # Parse Args.
    args = parser.parse_args()
    #
    byteOrder = {True: 'LSB', False: 'MSB'}[args.LSB]
    #
    secondGeodatFile = resolveDefault(args, 'secondGeodatFile',
                                      f'{args.secondDir}/{args.geodatFile}')
    #
    myRegion = resolveRegion(args)
    #
    checkFiles(myRegion, secondGeodatFile, args, fp)
    #
    return args.azOffsets, args.offsetsDat, myRegion, args.geodatFile, \
        secondGeodatFile, args.syncDat, args.fastMask, not args.noVel, \
        byteOrder


def resolveRegion(args):
    '''
    Resolve region to either use default (GrIMP server only) or read a region
    file.
    '''
    #
    geo = geodatrxa(file=args.geodatFile, echo=False)
    region = args.region
    if region is None:
        if geo.corners[0][0] > 0:
            region = 'greenland'
        else:
            region = 'antarctica'
    print(region)
    #
    myRegion = s.defaultRegionDefs(region, regionFile=args.regionFile)
    if not hasattr(myRegion, 'region'):
        u.myerror('Region not defined for region and regionFile: '
                  f' {args.region} {args.regionFile}')
    # override region defaults
    if args.velMap is not None:
        myRegion.setRegionField('velMap', args.velMap)
    if args.dem is not None:
        myRegion.setRegionField('dem', args.dem)
    if args.maskInputFile is not None:
        if args.fastMask:
            myRegion.setRegionField('fastmask', args.maskInputFile)
        else:
            myRegion.setRegionField('mask', args.maskInputFile)
    return myRegion


def checkFiles(myRegion, secondDir, args, fp):
    if myRegion.velMap() is not None:
        print('velMap= ', myRegion.velMap())
        if not os.path.exists(myRegion.velMap()+'.vx'):
            printError('velMap not found: ', myRegion.velMap(), fp)
    # check dem
    print('dem = ', myRegion.dem())
    if not os.path.exists(myRegion.dem()):
        printError('dem not found: ', myRegion.dem(), fp)
    # check secondDir
    print('secondDir = ', secondDir)
    if not os.path.exists(secondDir):
        printError('secondDir not found: ', secondDir, fp)
    # check geodat
    print('geodatFile= ', args.geodatFile)
    if not os.path.exists(args.geodatFile):
        printError('geodatFile  not found: ', args.geodatFile, fp)
    # check geodat
    maskFile = {False: myRegion.mask(),
                True: myRegion.fastmask()}[args.fastMask]
    print('maskInputFile= ',  maskFile)
    if not os.path.exists(maskFile):
        printError('maskInputFile  not found: ', maskFile, fp)


def simOffsetsProcessArgs(fp):
    """
    Process arguments and defaults
    """
    #
    # find second dir default, relies on second slc to be a pointer
    #
    args = sys.argv[1:]
    for arg in args:
        if len(sys.argv) > 1:
            if '-help' in arg:
                simoffsetsUsage()
    # Find second image as slc that is a link with ../
    slcs = u.dols('ls -l *.slc')
    second = [s for s in slcs if '../' in s]
    if len(second) == 1:
        secondDir = '../'+second[0].split('/')[1]
    else:
        secondDir = '.'
    syncDat = False
    fastMask = False
    secondGeodatFile = None
    # frame = int(secondDir.split('_')[-1])
    #
    azOffsets = 'offsets.da'
    offsetsDat = 'offsets.dat'
    geodatFile = 'geodat10x2.in'
    #
    useVel = True
    region = None
    velMap, dem, maskInputFile = None, None, None
    byteOrder = 'MSB'
    for arg in args:
        if len(sys.argv) > 1:
            if '-help' in arg:
                simoffsetsUsage()
            elif '-azOffsets' in arg:
                tmp = arg.split('=')
                azOffsets = tmp[1]
            elif '-offsetsDat' in arg:
                tmp = arg.split('=')
                offsetsDat = tmp[1]
            elif 'syncDat' in arg:
                syncDat = True
            elif 'noVel' in arg:
                useVel = False
            elif 'LSB' in arg:
                byteOrder = 'LSB'
            elif 'fastMask' in arg:
                fastMask = True
            elif '-velMap' in arg:
                tmp = arg.split('=')
                velMap = tmp[1]
            elif '-dem' in arg:
                tmp = arg.split('=')
                dem = tmp[1]
            elif '-geodatFile' in arg:
                tmp = arg.split('=')
                geodatFile = tmp[1]
            elif '-secondDir' in arg:
                tmp = arg.split('=')
                secondDir = tmp[1]
            elif '-secondGeodatFile' in arg:
                tmp = arg.split('=')
                secondGeodatFile = tmp[1]
            elif '-region' in arg:
                tmp = arg.split('=')
                region = tmp[1]
            elif '-maskInputFile' in arg:
                tmp = arg.split('=')
                maskInputFile = tmp[1]
            else:
                print('\nInvalid argument ', arg)
                simoffsetsUsage()
    #
    #
    # Echo inputs and check existence
    #
    #
    # select epsg based on lat
    geo = geodatrxa(file=geodatFile, echo=False)
    if region is None:
        if geo.corners[0][0] > 0:
            region = 'greenland'
        else:
            region = 'antarctica'
    print(region)
    #
    myRegion = s.defaultRegionDefs(region)
    # override region defaults
    if velMap is not None:
        myRegion.setRegionField('velMap', velMap)
    if dem is not None:
        myRegion.setRegionField('dem', dem)
    if maskInputFile is not None:
        if fastMask:
            myRegion.setRegionField('fastmask', maskInputFile)
        else:
            myRegion.setRegionField('mask', maskInputFile)
    #
    # check velMap
    if myRegion.velMap() is not None:
        print('velMap= ', myRegion.velMap())
        if not os.path.exists(myRegion.velMap()+'.vx'):
            printError('velMap not found: ', myRegion.velMap(), fp)
    else:
        useVel = False
    # check dem
    print('dem = ', myRegion.dem())
    if not os.path.exists(myRegion.dem()):
        printError('dem not found: ', myRegion.dem(), fp)
    # check secondDir
    print('secondDir = ', secondDir)
    if not os.path.exists(secondDir):
        printError('secondDir not found: ', secondDir, fp)
    # check geodat
    print('geodatFile= ', geodatFile)
    if not os.path.exists(geodatFile):
        printError('geodatFile  not found: ', geodatFile, fp)
    # check geodat
    print('maskInputFile= ', [myRegion.mask(), myRegion.fastmask()][fastMask])
    if not os.path.exists([myRegion.mask(), myRegion.fastmask()][fastMask]):
        printError('maskInputFile  not found: ',
                   [myRegion.mask(), myRegion.fastmask()][fastMask], fp)
    #
    if secondGeodatFile is None:
        secondGeodatFile = f'{secondDir}/{geodatFile}'
    #
    return azOffsets, offsetsDat, myRegion, geodatFile, secondGeodatFile, \
        syncDat, fastMask, useVel, byteOrder


def LLtoRA(lat, lon, geodatFile, dem=None):
    if dem is None:
        dem = '/Volumes/insar7/ian/gimp/gimp2/270m/dem.gimp2.270m'
    # bail if no dem file
    pid = os.getpid()
    if not os.path.exists(dem):
        print('LLtoRA: DEM file does not exist ')
    # write lat/lon to temp file
    u.writeLLtoRAformat(lat, lon, tempfile=f'temp.{pid}.ll')
    # execute lltora
    command = f'lltora {geodatFile} {dem} temp.{pid}.ll temp.{pid}.ra'
    #
    call(command, shell=True, executable='/bin/csh')
    # read and reformat data
    r, a = u.readLLtoRA(tempfile=f'temp.{pid}.ra')
    r = np.reshape(r, lat.shape)
    a = np.reshape(a, lat.shape)
    # clean up
    os.remove(f'temp.{pid}.ra')
    os.remove(f'temp.{pid}.ll')
    return r, a


def groundToSlantRangeResolution(offsets):
    """ compute conversion from slant to ground range resolution
    This is a fairly crude approximation - e.g, earth curvature not included
    it is close enough for scaling offsets    """
    # get coordinates
    r, a = offsets.getRACoords()
    slpR, slpA = offsets.geodatrxa.singleLookResolution()
    # get geometry
    H = offsets.geodatrxa.satelliteAltm()
    Rc = offsets.geodatrxa.centerRangem()
    Re = offsets.geodatrxa.earthRadm()
    # single look pixel dimension
    nrls = offsets.geodatrxa.nlr * offsets.geodatrxa.nr
    # slant range
    R = (r-nrls/2) * slpR + Rc
    # compute look angle
    thetaInc = R**2 + 2 * H * Re + H**2
    thetaInc = thetaInc/(2*R*(Re+H))
    thetaInc = np.arccos(thetaInc)
    # compute incidence angle
    sinPhic = (Re+H)/Re * np.sin(thetaInc)
    groundToSlant = sinPhic/slpR
    return groundToSlant


def fixBad(d, coeff, r1, a1):
    ''' Fix potential outliers caused by DEM or other issues '''
    # evaluate polynomial
    dp = coeff[0]+r1*coeff[1] + a1*coeff[2]
    # gets really bad points with this threshold
    iBad = np.absolute(d-dp) > 2
    d[iBad] = dp[iBad]


def computeStaticOffsets(offsets1, secondGeoDat, dem):
    """ compute offsets based on geometry and topography - i.e. no motion """
    lat1, lon1 = offsets1.getLatLon()
    # get r1 coords in first and second image
    r1, a1 = offsets1.getRACoords()
    print('DEM', dem)
    r12, a12 = LLtoRA(lat1, lon1, secondGeoDat, dem=dem)
    # determing missing points
    missing12 = r12 < 0
    print('missing', np.sum(missing12))
    # compute difference
    dr0 = r12-r1
    dr0[missing12] = -2.0e9
    #
    da0 = a12-a1
    da0[missing12] = -2.0e9
    # compute linear polynomials for offset fit
    coeffR, coeffA = computePoly(dr0, da0, r1, a1)
    # u.myerror('stop')
    fixBad(dr0, coeffR, r1, a1)
    fixBad(da0, coeffA, r1, a1)
    return dr0, da0, coeffR, coeffA


def computePlane(x, y, z):
    # cull bad
    good12 = np.logical_and(z > -1.99e9, x >= 0)
    # flatten XY
    X = x[good12].flatten()
    Y = y[good12].flatten()
    # poly values
    A = np.array([X*0+1, X, Y]).T
    # diffs
    Brg = z[good12].flatten()
    # compute polynomials
    coeffXY, rR, rankR, sR = np.linalg.lstsq(A, Brg, rcond=None)
    return coeffXY


def computePoly(dr0, da0, r1, a1):
    # cull bad
    coeffR = computePlane(r1, a1, dr0)
    coeffA = computePlane(r1, a1, da0)
    return coeffR, coeffA


def writeOffPoly(fileName, coeffR, coeffA):
    fp = open(fileName, 'w')
    print(f'# range poly \n{coeffR[0]:7e} {coeffR[1]:7e} {coeffR[2]:7e}',
          file=fp)
    print(f'# azimuth poly \n{coeffA[0]:7e} {coeffA[1]:7e} {coeffA[2]:7e}',
          file=fp)


def computeVelocityRA(velMap, offsets1, srsInfo):
    """
    Read and interpolate velMap for points lat1, lon1 and return in radar
    coordinates
    """
    # get velocity
    vel = u.geoimage(geoType='velocity')
    vel.readData(velMap, epsg=srsInfo['epsg'], wktFile=srsInfo['wktFile'])

    # compute xy angle, heading, and angle for rotation
    lat1, lon1 = offsets1.getLatLon()
    xps1, yps1 = vel.geo.lltoxykm(lat1, lon1)
    xyAngle = np.arctan2(-yps1, -xps1)
    rotAngle = (offsets1.computeHeading() + np.pi/2.) - xyAngle
    #
    #  interpolate velocity
    vel.setupInterp()
    vxr, vyr, vr = vel.interpGeo(xps1, yps1)
    vr[np.isnan(vr)] = 0
    #
    # set fast regions to 8 to flag
    #
    mask = offsets1.getMask()
    if len(mask) > 0:
        fast = np.logical_and(vr > 200, mask > 0)
        mask[fast] = mask[fast] | 8
    #
    # Do rotation
    #
    cosRot = np.cos(rotAngle)
    sinRot = np.sin(rotAngle)
    vr = vxr * cosRot - vyr * sinRot
    va = vxr*sinRot + vyr * cosRot
    #
    return vr, va


def printError(msg, var, fp):
    '''
    Print error message
    '''
    print(msg, var)
    print(msg, var, file=fp)
    exit()


def runSim(geodatFile, offsetsDat, dem, maskInputFile, syncDat,
           byteOrder='MSB'):
    """
    runSim - run simulation with siminsar
    If no offsets.lat/lon/dat, create them
    """
    byteOrderFlag = {'MSB': '', 'LSB': '-LSB'}[byteOrder]
    geo = geodatrxa(file=geodatFile, echo=False)
    #
    # default values
    if not syncDat:
        r0, a0 = 180, 180
        dr, da = 24, 18
        nr, na = geo.singleLookSize()
        sr = int(((nr-2*r0)/dr)/6)*6
        sa = int(((na-2*a0)/da)/6)*6
        fpDat = open(offsetsDat, 'w')
        print(r0, a0, sr, sa, dr, da, file=fpDat)
        fpDat.close()
    else:
        # use values from specified offsets.dat file
        if not os.path.exists(offsetsDat):
            u.myerror(f'siminsar.py - runSim - missing {offsetsDat:s} ')
    #
    offsetsRoot = offsetsDat.replace('.dat', '').replace('.vrt', '')
    command = f'siminsar {byteOrderFlag} -center -toLL {offsetsDat} -mask ' \
        f'-xyDEM {dem}  {maskInputFile} {geodatFile} {offsetsRoot}'
    call(command, shell=True,  executable='/bin/csh')


def main():
    """
    Simulate offsets for input to strack programs
    """
    #
    # this will indicate a fail, unless program completes and removes
    #
    fp = open('fail.simoffsets', 'w')
    azOffsets, offsetsDat, myRegion, geodatFile, secondGeoDatFile, syncDat, \
        fastMask, useVel, byteOrder = simOffsetsProcessArgs1(fp)
    for key in myRegion.region:
        print(myRegion.region[key])
    #
    if not os.path.exists(secondGeoDatFile):
        printError('Second geodatFile  not found: ', secondGeoDatFile, fp)
        exit()

    #
    # run simoff
    #
    latFile = offsetsDat.replace('.dat', '.lat').replace('.vrt', '.lat')
    lonFile = offsetsDat.replace('.dat', '.lon').replace('.vrt', '.lon')
    maskFile = offsetsDat.replace('.dat', '.mask').replace('.vrt', '.mask')
    maskInputFile = [myRegion.mask(), myRegion.fastmask()][fastMask]
    #
    if not os.path.exists(latFile) or not os.path.exists(lonFile) \
            or not os.path.exists(maskFile) \
            or not os.path.exists(offsetsDat) or syncDat:
        runSim(geodatFile, offsetsDat, myRegion.dem(), maskInputFile, syncDat,
               byteOrder=byteOrder)
    #
    # load offsets
    #
    latLonRoot = azOffsets.replace('.da', '')
    #
    vrtFile = None
    if 'vrt' in offsetsDat:
        vrtFile = offsetsDat
    print('++++++++', latLonRoot, offsetsDat, azOffsets, vrtFile)
    offsets1 = u.offsets(fileRoot=azOffsets, latlon=latLonRoot,
                         datFile=offsetsDat, geodatrxaFile=geodatFile,
                         maskFile=maskFile, vrtFile=vrtFile)
    #
    secondGeodatRxA = geodatrxa(file=secondGeoDatFile, echo=False)
    # compute deltaT
    deltaT = secondGeodatRxA.datetime.date() - \
        offsets1.geodatrxa.datetime.date()
    deltaT = deltaT.days
    #
    # Compute offsets with no velocity
    dr0, da0, coeffR, coeffA = computeStaticOffsets(offsets1, secondGeoDatFile,
                                                    myRegion.dem())
    #
    writeOffPoly(offsetsDat.replace('.dat', '.poly'), coeffR, coeffA)
    #
    # compute motion offsets
    #
    if useVel:
        vr, va = computeVelocityRA(myRegion.velMap(), offsets1,
                                   myRegion.srsInfo())
    else:
        vr, va = np.zeros(dr0.shape), np.zeros(dr0.shape)
    #
    # Combine results
    #
    groundToSlant = groundToSlantRangeResolution(offsets1)
    drv = vr * groundToSlant * deltaT / 365.
    dr = dr0
    imissing = dr0 < -2.e8
    dr[np.isfinite(vr)] += drv[np.isfinite(vr)]
    dr[imissing] = -2.e9
    #
    # merge offsets
    #
    slpR, slpA = offsets1.geodatrxa.singleLookResolution()
    #
    dav = (va * deltaT / 365.) / slpA
    da = da0
    da[np.isfinite(va)] += dav[np.isfinite(va)]
    da[imissing] = -2.0e9
    offsets1.azOff = da
    offsets1.rgOff = dr
    # output
    offsets1.writeOffsets(fileRoot=azOffsets, noDatFiles=True,
                          byteOrder=byteOrder)
    #
    offsets1.writeOffsetVrt(offsetsDat.replace('.dat', '.vrt'), [
                            os.path.basename(
                                offsetsDat.replace('.dat', '.dr')),
                            os.path.basename(
                                offsetsDat.replace('.dat', '.da'))],
                            ['RangeOffsets', 'AzimuthOffsets'],
                            byteOrder=byteOrder,
                            additionalMetaData={'deltaT': deltaT})
    #
    # remove flag file
    #
    fp.close()
    if os.path.isfile('fail.simoffsets'):
        os.remove('fail.simoffsets')


if __name__ == "__main__":
    main()
