#!/usr/bin/env python3
import utilities as u
import numpy as np
import sys
import os
import shutil
import tempfile
from subprocess import call
from utilities import geodatrxa
import sarfunc as s
import argparse
import glob
from osgeo import gdal


def _writeArrayAsTiff(filename, data, noDataValue=-2.e9, compress=True):
    """Write float32 numpy array to GeoTIFF with pixel-coord geotransform."""
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(filename, data.shape[1], data.shape[0], 1, gdal.GDT_Float32,
                       options=['COMPRESS=DEFLATE'] if compress else [])
    ds.SetGeoTransform([-0.5, 1., 0., -0.5, 0., 1.])
    band = ds.GetRasterBand(1)
    band.WriteArray(data.astype(np.float32))
    band.SetNoDataValue(noDataValue)
    ds = None


def _writeByteArrayAsTiff(filename, data):
    """Write byte (uint8) numpy array to GeoTIFF with pixel-coord geotransform."""
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(filename, data.shape[1], data.shape[0], 1, gdal.GDT_Byte,
                       options=['COMPRESS=DEFLATE'])
    ds.SetGeoTransform([-0.5, 1., 0., -0.5, 0., 1.])
    band = ds.GetRasterBand(1)
    band.WriteArray(data.astype(np.uint8))
    ds = None


def _maskedBoxMean(arr, valid, size):
    """
    Masked box mean with window (size, size), via scipy.ndimage.uniform_filter (an O(1)-
    per-pixel running-sum mean, the Python-side equivalent of computeSmoothRadius.c's
    boxAverage). Filtering value*mask and mask separately and dividing cancels the window-area
    normalization uniform_filter applies, leaving a true masked average.
    """
    from scipy.ndimage import uniform_filter
    validF = valid.astype(np.float64)
    num = uniform_filter(np.where(valid, arr, 0.0), size=size, mode='constant', cval=0.0)
    den = uniform_filter(validF, size=size, mode='constant', cval=0.0)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den > 0, num / den, arr)


def computeSmoothRadiusMapC(dr, da, toleranceDr, toleranceDa, maxRadiusR, maxRadiusA,
                            nIter, threads, outputTif):
    """
    Same map as computeSmoothRadiusMap(), computed by the smoothradius binary
    (mosaicSource/smoothRadius) instead of in Python. The sweep is ~7x faster serially and
    the independent radii run in parallel, which matters because this step dominated the
    ROFF stage (~873 s of a ~884 s frame at maxRadius 50).

    The binary reads and writes rasters, so dr/da/tolerances go out to a scratch directory
    (TMPDIR, i.e. local disk -- not the frame directory, which is usually NFS) and the
    resulting <base>.smr.tif is moved into place as outputTif.

    Returns True if the binary produced the map; False if it is missing or failed, so the
    caller can fall back to the Python and still finish the product.
    """
    tmpDir = tempfile.mkdtemp(prefix='smr.')
    try:
        for arr, name in ((dr, 'dr'), (da, 'da'),
                          (toleranceDr, 'tolDr'), (toleranceDa, 'tolDa')):
            _writeArrayAsTiff(f'{tmpDir}/{name}.tif', np.asarray(arr, dtype=np.float32),
                              compress=False)
        command = f'smoothradius -maxRadiusR {maxRadiusR} -maxRadiusA {maxRadiusA} ' \
            f'-nIter {nIter} -ompThreads {threads} ' \
            f'{tmpDir}/dr.tif {tmpDir}/da.tif {tmpDir}/tolDr.tif {tmpDir}/tolDa.tif ' \
            f'{tmpDir}/out'
        print(command)
        status = call(command, shell=True)
        if status != 0 or not os.path.exists(f'{tmpDir}/out.smr.tif'):
            print(f'WARNING: smoothradius failed (status {status}); falling back to the '
                  'Python smoothing-radius sweep, which is far slower')
            return False
        shutil.move(f'{tmpDir}/out.smr.tif', outputTif)
        return True
    finally:
        shutil.rmtree(tmpDir, ignore_errors=True)


def computeSmoothRadiusMap(dr, da, toleranceDr, toleranceDa, maxRadius, nIter):
    """
    Per-pixel single-look azimuth-pixel half-width: the largest contiguous radius r
    (1..maxRadius) for which repeatedly (nIter times, box->triangular->Gaussian-ish) box-
    filtering dr/da with half-width (r, r) changes the pixel by no more than
    toleranceDr/toleranceDa. First-violation cutoff: once a pixel fails at some r, its answer
    is locked at r-1 and never revisited (same convention as computeSmoothRadiusMap in
    mosaicSource/simInSAR/computeSmoothRadius.c). A pixel is locked the moment *either* dr or
    da exceeds its tolerance, giving one combined (more conservative) radius map.
    """
    valid = (dr > -1.e6) & (da > -1.e6) & np.isfinite(dr) & np.isfinite(da)
    radius = np.zeros(dr.shape, dtype=np.uint8)
    locked = ~valid
    for r in range(1, maxRadius + 1):
        size = 2 * r + 1
        sDr, sDa = dr, da
        for _ in range(nIter):
            sDr = _maskedBoxMean(sDr, valid, size)
            sDa = _maskedBoxMean(sDa, valid, size)
        ok = (~locked) & (np.abs(sDr - dr) <= toleranceDr) & (np.abs(sDa - da) <= toleranceDa)
        radius[ok] = r
        locked = locked | ~ok
        if locked.all():
            break
    return radius


def _writeTiffOffsetVrt(vrtFile, drTif, daTif, meta, additionalMetaData=None):
    """Write a two-band VRT referencing .dr.tif and .da.tif offset files.

    Built by hand (not gdal.BuildVRT, which rejects the "positive NS
    resolution" geotransform used by _writeArrayAsTiff and the other
    GrIMP tiff outputs).
    """
    if additionalMetaData:
        meta = {**meta, **{str(k): str(v) for k, v in additionalMetaData.items()}}
    srcDs = gdal.Open(drTif)
    xSize, ySize = srcDs.RasterXSize, srcDs.RasterYSize
    geoTransform = srcDs.GetGeoTransform()
    dataType = srcDs.GetRasterBand(1).DataType
    srcDs = None

    driver = gdal.GetDriverByName('VRT')
    vrt = driver.Create(vrtFile, xSize, ySize, 0)
    vrt.SetGeoTransform(geoTransform)
    vrt.SetMetadata({str(k): str(v) for k, v in meta.items()})

    vrtDir = os.path.dirname(os.path.abspath(vrtFile))
    for tif, description in ((drTif, 'RangeOffsets'), (daTif, 'AzimuthOffsets')):
        vrt.AddBand(dataType)
        band = vrt.GetRasterBand(vrt.RasterCount)
        band.SetMetadataItem('Description', description)
        sourceXml = (
            '<SimpleSource>'
            f'<SourceFilename relativeToVRT="1">{os.path.relpath(tif, vrtDir)}</SourceFilename>'
            '<SourceBand>1</SourceBand>'
            f'<SrcRect xOff="0" yOff="0" xSize="{xSize}" ySize="{ySize}"/>'
            f'<DstRect xOff="0" yOff="0" xSize="{xSize}" ySize="{ySize}"/>'
            '</SimpleSource>'
        )
        band.SetMetadataItem('source_0', sourceXml, 'new_vrt_sources')
    vrt = None


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

def simOffsetsProcessArgs1():
    """
    Process arguments and defaults
    """
    # Resolve slcs other wise default to current dir
    slcs = [os.readlink(x) if os.path.islink(x) else x
            for x in glob.glob('*.slc')]
    second = [s for s in slcs if '../' in s]
    if len(second) == 1:
        secondDirDefault = '../'+second[0].split('/')[1]
    else:
        secondDirDefault = '.'

    epilog = 'Notes: 1) Run in directory for the desired image pair ' \
        '2) Default region is greenland for NH or SH as determined from ' \
        'geodat 3) produces the following: offsets.da(.da.dat, .dr, ' \
        '.dr.dat, .dat, .lat, .lon, .mask, .simdat) ' \
        '\nPart of the mosaicworkflow package.'
    parser = argparse.ArgumentParser(
        description='\033[1m\n Use a velocity map'
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
    parser.add_argument('--iceMask', '-iceMask', action='store_true',
                        default=False,
                        help='Use binary ice-extent mask (icemask field in region YAML, '
                        'values: 0=not ice, 1=ice) instead of tracking mask')
    parser.add_argument('--iceRockWaterMask', '-iceRockWaterMask', action='store_true',
                        default=False,
                        help='Use ice-rock-water mask (icerockwatermask field in region '
                        'YAML, values: water=0, rock=1, ice=2) instead of tracking mask')
    parser.add_argument('-LSB', '--LSB', action='store_true',
                        default=False,
                        help='Write output in LSB byte order [MSB]')
    parser.add_argument('-velMap', '--velMap', type=str,
                        default=None, help='Root name for velocity map '
                        '[from Region file or default regionDefs]')
    parser.add_argument('-dem', '--dem', type=str, default=None,
                        help='DEM [from Region file or default regionDefs]')
    parser.add_argument('-verticalCorrection', '--verticalCorrection', type=str,
                        default=None,
                        help='Vertical correction (submergence/emergence) grid, '
                        'm/yr ice-equivalent, scalar GeoTIFF, same grid as passed '
                        'to siminsar -verticalCorrection [None - no correction]')

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
    parser.add_argument('--ompThreads', type=int, default=4,
                        help='Number of OpenMP threads for siminsar [4]')
    parser.add_argument('--tiff', action='store_true', default=False,
                        help='Write lat/lon and offset outputs as GeoTIFF '
                        'instead of GrIMP binary flat files. '
                        'lat/lon → .lat.tif/.lon.tif + .ll.vrt; '
                        'offsets → .dr.tif/.da.tif + .vrt [False]')
    parser.add_argument('--minTol', type=float, default=None,
                        help='Variable smoothing-radius map: m/yr floor for the adaptive '
                        'tolerance clip(percentSpeed/100*speed, minTol, maxTol). Required '
                        'together with --percentSpeed/--maxTol to enable the map [None]')
    parser.add_argument('--percentSpeed', type=float, default=None,
                        help='Variable smoothing-radius map: percent of local speed (e.g. '
                        '1 = 1%%) used in the adaptive tolerance. Required together with '
                        '--minTol/--maxTol [None]')
    parser.add_argument('--maxTol', type=float, default=None,
                        help='Variable smoothing-radius map: m/yr ceiling for the adaptive '
                        'tolerance. Required together with --minTol/--percentSpeed [None]')
    parser.add_argument('--maxSmoothRadius', type=int, default=50,
                        help='Variable smoothing-radius map: sweep cap in single-look '
                        'pixels, clamped to <= 255 (byte output) [50]')
    parser.add_argument('--maxSmoothRadiusA', type=int, default=None,
                        help='Variable smoothing-radius map: azimuth sweep cap, if it is '
                        'to differ from the range cap --maxSmoothRadius. Offsets grids are '
                        'resampled to square pixels, so the default (same as '
                        '--maxSmoothRadius, an isotropic sweep) is normally right [None]')
    parser.add_argument('--smoothThreads', type=int, default=2,
                        help='Variable smoothing-radius map: threads for the smoothradius '
                        'binary. Memory is ~0.9 GB + 0.3 GB per thread on a 13 Mpixel '
                        'grid [2]')
    parser.add_argument('--smoothNIter', type=int, default=3,
                        help='Variable smoothing-radius map: repeated box-filter passes '
                        'per sweep step (Gaussian-ish) [3]')

    #
    # Parse Args.
    args = parser.parse_args()
    smoothFlags = [args.minTol is not None, args.percentSpeed is not None,
                  args.maxTol is not None]
    if any(smoothFlags) and not all(smoothFlags):
        u.myerror('simoffsets: --minTol/--percentSpeed/--maxTol must be given together')
    if args.maxSmoothRadiusA is None:
        args.maxSmoothRadiusA = args.maxSmoothRadius
    for radiusArg in ['maxSmoothRadius', 'maxSmoothRadiusA']:
        if getattr(args, radiusArg) > 255:
            print(f'WARNING: --{radiusArg} {getattr(args, radiusArg)} exceeds byte range, '
                  'clamping to 255')
            setattr(args, radiusArg, 255)
    #
    byteOrder = {True: 'LSB', False: 'MSB'}[args.LSB]
    #
    # Fail marker lives alongside the offsets product it describes, so
    # failures can be located per-product rather than in whatever directory
    # the caller happened to be in.
    offsetsRoot = args.offsetsDat.replace('.dat', '').replace('.vrt', '')
    failFile = os.path.join(os.path.dirname(os.path.abspath(offsetsRoot)),
                            f'fail.simoffsets.{os.path.basename(offsetsRoot)}')
    fp = open(failFile, 'w')
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
        byteOrder, args.ompThreads, args.iceRockWaterMask, args.iceMask, args.tiff, \
        args.verticalCorrection, fp, failFile, args.minTol, args.percentSpeed, \
        args.maxTol, args.maxSmoothRadius, args.maxSmoothRadiusA, \
        args.smoothThreads, args.smoothNIter


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
    if myRegion.velMap() is not None and args.noVel == False:
        print('velMap = ', myRegion.velMap())
        if '.vrt' in myRegion.velMap() or '.tif' in myRegion.velMap():
            velFileVx = \
                myRegion.velMap().replace('*', 'vx').replace('vv', 'vx')
            velFileVy = velFileVx.replace('vx', 'vy')
        else:
            velFileVx = f'{myRegion.velMap()}.vx'
            velFileVy = f'{myRegion.velMap()}.vy'
        if not os.path.exists(velFileVx) or not os.path.exists(velFileVy):
            print(velFileVx, velFileVy)
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
    # check mask
    if args.iceRockWaterMask:
        maskFile = myRegion.iceRockWaterMask()
        if maskFile is None:
            printError('--iceRockWaterMask: icerockwatermask not defined in region YAML', '', fp)
    elif args.iceMask:
        maskFile = myRegion.icemask()
        if maskFile is None:
            printError('--iceMask: icemask not defined in region YAML', '', fp)
    else:
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
    call(command, shell=True)  # executable='/bin/csh')
    # read and reformat data
    r, a = u.readLLtoRA(tempfile=f'temp.{pid}.ra')
    r = np.reshape(r, lat.shape)
    a = np.reshape(a, lat.shape)
    # clean up
    os.remove(f'temp.{pid}.ra')
    os.remove(f'temp.{pid}.ll')
    return r, a


# nodata sentinel for vertical-correction grids, matches mosaicSource/common/common.h MINVCORRECT
MINVCORRECT = -100


def computeSinPsi(offsets):
    """ sin of the local incidence angle psi = (Re+H)/Re * sin(look angle)
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
    return (Re+H)/Re * np.sin(thetaInc)


def groundToSlantRangeResolution(offsets):
    """ compute conversion from slant to ground range resolution
    This is a fairly crude approximation - e.g, earth curvature not included
    it is close enough for scaling offsets    """
    slpR, slpA = offsets.geodatrxa.singleLookResolution()
    sinPhic = computeSinPsi(offsets)
    groundToSlant = sinPhic/slpR
    return groundToSlant


def computeVerticalCorrectionOffset(vcFile, offsets1, srsInfo, deltaT):
    """
    LOS contribution of the vertical-correction (submergence/emergence) grid to
    the simulated range offset, in slant-range pixels.

    Mirrors mosaic3d's make3DOffsets.c, which removes this contribution from a
    measured range offset via
        dDelta += dzdtSubmergence * cos(psi) * nDays/365.25
    (in metres) before solving for horizontal velocity. For the forward
    simulation here, the same term is SUBTRACTED (same convention chosen for
    siminsar's -verticalCorrection phase case) so mosaic3d's += round-trips
    correctly. Vertical motion has no azimuth-direction LOS component, so only
    the range offset is affected.
    """
    vc = u.geoimage(geoType='scalar')
    vc.readData(vcFile, tiff=True, epsg=srsInfo['epsg'], wktFile=srsInfo['wktFile'])
    lat1, lon1 = offsets1.getLatLon()
    xc1, yc1 = vc.geo.lltoxykm(lat1, lon1)
    vc.setupInterp()
    dzdtSubmergence = vc.interpGeo(xc1, yc1)
    dzdtSubmergence = np.nan_to_num(dzdtSubmergence, nan=0.0)
    dzdtSubmergence[dzdtSubmergence <= MINVCORRECT] = 0.0
    slpR, slpA = offsets1.geodatrxa.singleLookResolution()
    cosPsi = np.sqrt(1. - computeSinPsi(offsets1)**2)
    return dzdtSubmergence * cosPsi * deltaT / 365.25 / slpR


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
    # Modern velMaps are a single GDAL file (.tif/.vrt, often with vx/vy as
    # separate bands -- see geoimage._multibandComponents); legacy velMaps
    # are flat binary + .vx.geodat/.vy.geodat sidecars, with no extension
    # to key off of.
    isTiff = velMap.lower().endswith(('.tif', '.vrt'))
    vel.readData(velMap, epsg=srsInfo['epsg'], wktFile=srsInfo['wktFile'],
                 tiff=isTiff)

    # compute xy angle, heading, and angle for rotation
    lat1, lon1 = offsets1.getLatLon()
    xps1, yps1 = vel.geo.lltoxykm(lat1, lon1)
    xyAngle = np.arctan2(-yps1, -xps1)
    rotAngle = (offsets1.computeHeading() + np.pi/2.) - xyAngle
    print('rotAngle', np.nanmean(rotAngle))
    #
    #  interpolate velocity
    vel.setupInterp()
    vxr, vyr, vr = vel.interpGeo(xps1, yps1)
    print('vx, vy', np.nanmean(vxr), np.nanmean(vyr), np.nanmean(vr))
    vr[np.isnan(vr)] = 0
    speed = vr  # ground-speed magnitude, before vr is overwritten by the rotation below
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
    return vr, va, speed


def printError(msg, var, fp):
    '''
    Print error message
    '''
    print(msg, var)
    print(msg, var, file=fp)
    exit()


def runSim(geodatFile, offsetsDat, dem, maskInputFile, syncDat,
           byteOrder='MSB', ompThreads=4, tiff=False):
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
    tiffFlag = '-tiff' if tiff else ''
    offsetsRoot = offsetsDat.replace('.dat', '').replace('.vrt', '')
    command = f'siminsar {byteOrderFlag} {tiffFlag} -ompThreads {ompThreads} ' \
        f'-center -toLL {offsetsDat} -mask ' \
        f'-xyDEM {dem}  {maskInputFile} {geodatFile} {offsetsRoot}'
    print(command)
    call(command, shell=True)  # executable='/bin/csh')


def main():
    """
    Simulate offsets for input to strack programs
    """
    #
    # this will indicate a fail, unless program completes and removes
    #
    azOffsets, offsetsDat, myRegion, geodatFile, secondGeoDatFile, syncDat, \
        fastMask, useVel, byteOrder, ompThreads, iceRockWaterMask, iceMask, tiff, \
        verticalCorrection, fp, failFile, minTol, percentSpeed, maxTol, \
        maxSmoothRadius, maxSmoothRadiusA, smoothThreads, smoothNIter = \
        simOffsetsProcessArgs1()
    #
    for key in myRegion.region:
        print(myRegion.region[key])
    #
    if not os.path.exists(secondGeoDatFile):
        printError('Second geodatFile  not found: ', secondGeoDatFile, fp)
        exit()
    #
    # run simoff
    #
    maskFile = offsetsDat.replace('.dat', '.mask').replace('.vrt', '.mask')
    if tiff:
        llVrt = offsetsDat.replace('.dat', '.ll.vrt').replace('.vrt', '.ll.vrt')
        # siminsar -tiff writes the mask as <root>.mask.tif + <root>.mask.vrt,
        # so testing for the raw <root>.mask would re-run the sim every time.
        needSim = not os.path.exists(llVrt) \
                  or not os.path.exists(f'{maskFile}.vrt') \
                  or not os.path.exists(offsetsDat) or syncDat
    else:
        latFile = offsetsDat.replace('.dat', '.lat').replace('.vrt', '.lat')
        lonFile = offsetsDat.replace('.dat', '.lon').replace('.vrt', '.lon')
        needSim = not os.path.exists(latFile) or not os.path.exists(lonFile) \
                  or not os.path.exists(maskFile) \
                  or not os.path.exists(offsetsDat) or syncDat
    if iceRockWaterMask:
        maskInputFile = myRegion.iceRockWaterMask()
    elif iceMask:
        maskInputFile = myRegion.icemask()
    else:
        maskInputFile = [myRegion.mask(), myRegion.fastmask()][fastMask]
    #
    if needSim:
        runSim(geodatFile, offsetsDat, myRegion.dem(), maskInputFile, syncDat,
               byteOrder=byteOrder, ompThreads=ompThreads, tiff=tiff)
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
    if tiff:
        # getMask() below reads maskVrtFile when set, otherwise the raw maskFile
        # -- which siminsar -tiff no longer writes.
        offsets1.maskVrtFile = f'{maskFile}.vrt'
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
        print('Using velocity')
        vr, va, speed = computeVelocityRA(myRegion.velMap(), offsets1,
                                          myRegion.srsInfo())
        if offsets1.geodatrxa.lookdir == 'left':
            vr *= -1
    else:
        print('No velocity')
        vr, va = np.zeros(dr0.shape), np.zeros(dr0.shape)
        speed = np.zeros(dr0.shape)
    #
    # Combine results
    #
    groundToSlant = groundToSlantRangeResolution(offsets1)
    drv = vr * groundToSlant * deltaT / 365.
    dr = dr0
    imissing = dr0 < -2.e8
    dr[np.isfinite(vr)] += drv[np.isfinite(vr)]
    if useVel and verticalCorrection is not None:
        print('Using vertical correction', verticalCorrection)
        drVCorrect = computeVerticalCorrectionOffset(verticalCorrection, offsets1,
                                                      myRegion.srsInfo(), deltaT)
        dr[~imissing] -= drVCorrect[~imissing]
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
    #
    # Variable smoothing-radius map (--minTol/--percentSpeed/--maxTol)
    #
    if minTol is not None:
        Xv = np.clip(percentSpeed / 100. * speed, minTol, maxTol)  # m/yr
        toleranceDr = Xv * groundToSlant * deltaT / 365.  # slant-range pixels
        toleranceDa = Xv * deltaT / 365. / slpA  # azimuth pixels
        smrTif = offsetsDat.replace('.dat', '.smr.tif')
        if not computeSmoothRadiusMapC(dr, da, toleranceDr, toleranceDa,
                                       maxSmoothRadius, maxSmoothRadiusA, smoothNIter,
                                       smoothThreads, smrTif):
            radius = computeSmoothRadiusMap(dr, da, toleranceDr, toleranceDa,
                                            maxSmoothRadius, smoothNIter)
            _writeByteArrayAsTiff(smrTif, radius)
    # output
    extraMeta = {'deltaT': deltaT}
    if useVel and verticalCorrection is not None:
        extraMeta['verticalCorrection'] = verticalCorrection
    if tiff:
        drTif = os.path.abspath(offsetsDat.replace('.dat', '.dr.tif'))
        daTif = os.path.abspath(offsetsDat.replace('.dat', '.da.tif'))
        _writeArrayAsTiff(drTif, dr.astype(np.float32))
        _writeArrayAsTiff(daTif, da.astype(np.float32))
        if offsets1.meta is None:
            offsets1.genMeta()
        _writeTiffOffsetVrt(offsetsDat.replace('.dat', '.vrt'),
                            drTif, daTif,
                            offsets1.meta,
                            additionalMetaData=extraMeta)
    else:
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
                                additionalMetaData=extraMeta)
    #
    # remove flag file
    #
    fp.close()
    if os.path.isfile(failFile):
        os.remove(failFile)


if __name__ == "__main__":
    main()
