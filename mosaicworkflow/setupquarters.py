#!/usr/bin/env python3
import utilities as u
import os
import geopandas as gpd
import sarfunc as s
from subprocess import call
from datetime import datetime, timedelta, date
import numpy as np
import threading
import shapefile
import shutil
import glob
import argparse
import mosaicfunc as mosf
from matplotlib import colors
from osgeo import gdal, osr
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
import gc
# Currently this is hardwired here since this program only works for greenland.
# It could be used for future expansion to Antarctica I have tried to use it to
# flag ice sheet specific stuff

iceSheetRegion = 'greenland'
iceSheetRegionFile = None
myProjWin = None

regionDefs = None
outputMaskTemplate = None
regionID = None

currentVersions = {'s1cycle': 2, 's1-12day': 2, 'Monthly': 5, 'Quarterly': 5,
                   'Annual': 5, 'multiYear': 1}
subVersions = {'s1cycle': 0, 'Monthly': 0, 'Quarterly': 0, 'Annual': 0,
               'multiYear': 0}


def processTemplate(args):
    '''
    Merge template inputs with commandline overrides.
    '''
    # Read template
    template = mosf.readVelocityTemplate(args.template)
    template['template'] = args.template
    # Command line will override dem
    for key in ['dem', 'lsFile', 'baseFlags', 'outputMask', 'mosaicMask',
                'inputFile']:
        try:
            # Override template if passed at command line
            # print(key, args)
            if key not in template:
                template[key] = None  # Aovid print fail, reset below
            if getattr(args, key) is not None:
                print(f'Overriding template {key} value of {template[key]}'
                      f'command line value {getattr(args, key)}')
                template[key] = getattr(args, key)
            else:
                if key not in template:
                    template[key] = None
        except Exception:
            u.mywarning(f'processTemplate: invalid key {key}')
        #
    # Check all other essential keys found
    for key in ['xll', 'yll', 'sx', 'sy', 'dx', 'dy', 'regionID', 'dem']:
        if key not in template:
            u.myerror(f'Missing parameter {key} in both template and command')
    #
    # Compute corners
    for c in ['x', 'y']:
        template[f'{c}min'] = template[f'{c}ll'] - template[f'd{c}'] * 0.5
        template[f'{c}max'] = template[f'{c}min'] + \
            template[f's{c}'] * template[f'd{c}']
    #
    if 'wktFile' in template and 'epsg' not in template:
        template['epsg'] = None
    return template


def processFlags(args, template):
    '''
    Create dictionary with relevant flags
    '''
    flags = {}
    # Dates
    for flag in ['firstdate', 'lastdate', 'noTSX', 'noCull', 'check',
                 'timeOverlapOff', 'noReprocess']:
        if flag in args:
            flags[flag] = getattr(args, flag)
        else:
            u.mywarning(f'processFlags: missing value for value {flag}')
    #
    try:
        flags['firstdate'] = datetime.strptime(flags['firstdate'], '%Y-%m-%d')
        flags['lastdate'] = datetime.strptime(flags['lastdate'], '%Y-%m-%d')
    except Exception:
        u.myerror('Could not process one or both dates '
                  f'{flags["firstdate"]} {flags["lastdate"]}')
    if template['baseFlags'] is None:
        u.myerror('processFlags: baseFlags not specified.')
    elif not flags['timeOverlapOff']:
        template['baseFlags'] += ' -timeOverlap '
    return flags


def extractRegion(template):
    '''
    Extract info to define regions
    '''
    region = {}

    for key in ['epsg', 'wktFile', 'dem', 'velMap', 'sigmaShape']:
        if key in template:
            region[key] = template[key]
        else:
            region[key] = None
            # Can miss epsg is if wktfile
            if not (key == 'epsg' and 'wktFile' in template):
                u.mywarning(f'Region def missing {key}')
    return region


def processQArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\033[1;34m Setup individual velocity mosaics \033[0m',
        epilog='Notes:  ', allow_abbrev='False')
    # flag
    parser.add_argument('--firstdate', type=str, default='2015-01-01',
                        help='First date for series of mosaics. '
                        'Use central dates >= first date [2015-01-01]')
    defaultEnd = date.today().strftime('%Y-%m-%d')
    # flag
    parser.add_argument('--lastdate', type=str, default=defaultEnd,
                        help='Last date for series of mosaics. '
                        'Use central dates <= lastdate [defaultDend]')
    #
    parser.add_argument('--noReprocess', action='store_true', default=False,
                        help='Just reformat data')
    # flag
    parser.add_argument('--noTSX', action='store_true', default=False,
                        help='Exclude TSX/CSK data')
    # flag
    parser.add_argument('--noCull', action='store_true', default=False,
                        help='Use unculled offsets')
    # flag
    parser.add_argument('-timeOverlapOff', action='store_true', default=False,
                        help='Set to turn timeOverlap off')
    # template
    parser.add_argument('--baseFlags', type=None,
                        default=None,
                        help='Override default flags for the mosaicker')
    # template
    parser.add_argument('--outputMask', type=str, default=None,
                        help='Shape mask for final output (shp)')
    parser.add_argument('--mosaicMask', type=str, default=None,
                        help='Mask passed to mosaic3d (shelfMask)')
    # flag
    parser.add_argument('--check', action='store_true', default=False,
                        help='Setup command,  check masks, but do not run')
    parser.add_argument('--dem', type=str, default=None,
                        help='Override dem from template')
    # template
    parser.add_argument('--inputFile', type=str, default=None,
                        help='Input file with SAR inputs')
    # template
    parser.add_argument('--lsFile', type=str, default=None,
                        help='Landsat input file for mosaic period')
    parser.add_argument('--template', type=str, default='mosaic.template.yaml',
                        help='template that defines mosaic')
    #
    args = parser.parse_args()
    template = processTemplate(args)
    flags = processFlags(args, template)
    #
    if template['mosaicMask'] is not None:
        template['baseFlags'] += f' -shelfMask {template["mosaicMask"]}'
    #
    if 'regionFile' in template:
        regionDefs = s.defaultRegionDefs(None,
                                         regionFile=template['regionFile'])
    else:
        region = extractRegion(template)
        regionDefs = s.defaultRegionDefs(None, regionDef=region)
    #
    return template, flags, regionDefs


def readWKT(wktFile):
    ''' get wkt from a file '''
    with open(wktFile, 'r') as fp:
        return fp.readline()

# -----------------------------------------------------------------------------
# Read in input file, save as list with [date, line] entries.
# -----------------------------------------------------------------------------


def processInputFileQ(filename, noTSX):
    try:
        inputFile = open(filename, 'r')
    except Exception:
        u.myerror(f'could not open SAR input file: {filename}')
    dataTakes = []
    count = 0
    for line in inputFile:
        # use line length to skip over non data stuff
        if ';' not in line and len(line) > 80:
            if noTSX and 'TSX' in line:
                continue
            # setup cross flag, 1 if no TSX, 0 if TSX
            crossFlag = ['1', '0']['TSX' in line]
            line = line.strip('\n')
            parts = line.split()
            # get geodat
            geo = parts[1].strip()
            exclude = "/".join(geo.split('/')[0:-2])+'/Exclude'
            if not os.path.exists(exclude):
                with open(geo, 'r') as fgeo:
                    # read geodat to extract date
                    for gline in fgeo:
                        if 'Image date' in gline:
                            gline = gline.split(':')[1]
                            gline = gline.strip('\n').replace(' ', '')
                            mydate = datetime.strptime(gline, '%d%b%Y')
                            break
                count += 1
                if count % 500 == 0:
                    print('.', end='')
                line += ' '+crossFlag
                dataTakes.append([mydate, line])
    inputFile.close()
    return dataTakes


def computeTSXSeasonCount(dataTakes, date1, date2):
    ''' Count number of TSX summer scenes for each glacier in the year
    defined by date1 and date2'''
    seasonCount = {}
    for dataTake in dataTakes:
        seasonKey = {12: 'w', 1: 'w', 2: 'w',
                     3: 'sp', 4: 'sp', 5: 'sp',
                     6: 's', 7: 's', 8: 's',
                     9: 'f', 10: 'f', 11: 'f'
                     }
        if 'TSX' in dataTake[1]:
            if dataTake[0] >= date1 and dataTake[0] < date2:
                glacierID = dataTake[1].split()[1].split('/')[-4]
                if glacierID not in seasonCount:
                    seasonCount[glacierID] = {'w': 0, 'sp': 0, 's': 0, 'f': 0}
                seasonCount[glacierID][seasonKey[dataTake[0].month]] += 1.0
    return seasonCount


def setLineWeight(dataTake, seasonCount):
    ''' Set weight to line 1/count. '''
    seasonKey = {12: 'w', 1: 'w', 2: 'w',
                 3: 'sp', 4: 'sp', 5: 'sp',
                 6: 's', 7: 's', 8: 's',
                 9: 'f', 10: 'f', 11: 'f'
                 }
    # return orig if no summer count
    if seasonCount is None:
        return dataTake[1]
    #
    glacierID = dataTake[1].split()[1].split('/')[-4]
    # Return line if glacier not in summerCount
    if glacierID not in seasonCount:
        return dataTake[1]
    season = seasonKey[dataTake[0].month]
    count = min(max(1., seasonCount[glacierID][season]), 50.)
    pieces = dataTake[1].split()
    pieces[4] = f'{1/count:0.3f}'
    newLine = ' '.join(pieces)
    # print(newLine)
    return newLine


def writeInputFileQ(infile, xdim, ydim, dataTakes, flags):
    #
    fOut = open(infile, 'w')
    print('; automatically generated input file from setupannual mosaics-- ',
          file=fOut)
    # print limit
    # print('{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}'.format(*limit),
    #      file=fOut)
    mosf.writeInputHeader(xdim, ydim, fp=fOut)
    # first pass through to filter by date, and get count
    n = 0
    newList = []
    #
    # this is part of a hack to deal with isolated tsx frames in the winter,
    # which can dominate. If there is a winter tsx frame, then only keep it if
    # is more than half way in output window (i.e., the middle date falls in
    # the win) this is the first test to determine whether to apply screening
    screenTSX = False
    # don't screen summer months since there should be multiple TSX then
    if flags['firstdate'].month not in [6, 7, 8] or \
            flags['lastdate'].month not in [7, 8, 9]:
        # if more than a 15 day estimate, apply screen
        if (flags['lastdate'] - flags['firstdate']).days > 15:
            screenTSX = True
            # print(screenTSX)
    productType = prodType(flags['firstdate'], flags['lastdate'])
    #

    # now loop through data takes, keeping ones in range
    if productType == 'Annual':
        TSXSeasonCount = computeTSXSeasonCount(dataTakes,
                                               flags['firstdate'],
                                               flags['lastdate'])
    #

    for dataTake in dataTakes:
        firstDate = dataTake[0]
        tmp = dataTake[1].split()
        phaseFile = tmp[0]
        secondDate = firstDate + timedelta(float(tmp[3]))
        midDate = firstDate + timedelta(float(tmp[3])) * 0.5
        #
        # Reduce weight of TSX summer frame so large number of frames does not
        # skew results in annual average.
        if 'TSX' in dataTake[1] and productType == 'Annual':
            dataTake[1] = setLineWeight(dataTake, TSXSeasonCount)
        #
        # if this a TSX case and screenTSX on
        # Added only images with center for 6/12 day pairs.
        # Not sure what the data range is different, might have to do with
        # .5 dT non integer.
        if ('TSX' in dataTake[1] and screenTSX):
            if midDate < flags['firstdate'] or midDate >= flags['lastdate']:
                continue
        # then check midDate in the range, otherwise
        # Fixed 6/7/22 - replaced <= date1 with < date1
        # use continue to skip this iteration
        # Fixed 2/7/24 (changed >= to >)
        if productType == 's1cycle':
            if midDate < flags['firstdate'] or midDate > flags['lastdate']:
                continue
        # This is setup to pass phase data from Oct/year-1 to March/year+1
        # As such, it will provide roughly symmetrical annual average.
        # Offsets will be preserved within years since they tend have full year
        # at coasts.
        if productType == 'Annual' and 'nophase' not in phaseFile:
            dT1, dT2 = timedelta(-61), timedelta(121)
        else:
            dT1, dT2 = timedelta(0), timedelta(0)
        # pass dates with any amount of overlap
        if secondDate >= (flags['firstdate'] + dT1) and \
                firstDate <= (flags['lastdate'] + dT2):
            if flags['noCull']:
                new = changeToNoCull(dataTake[1])
            else:
                new = [dataTake[1]]
            n += len(new)
            newList += new
    #
    # now write items to list file
    # number of lines
    print(len(newList), file=fOut)
    # loop for actual lines
    for item in newList:
        print(';\n', item, file=fOut)
    fOut.close()
    return


def changeToNoCull(line):
    mainDir = os.path.dirname(line.split()[1])
    az = sorted(glob.glob(f'{mainDir}/*.interp.da'))[0]
    azname = os.path.basename(az)
    rg = sorted(glob.glob(f'{mainDir}/*.interp.dr'))[0]
    rgname = os.path.basename(rg)
    dataTakeSlow = line.replace('azimuth.offsets', azname)
    dataTakeSlow = dataTakeSlow.replace('range.offsets', rgname)
    new = [dataTakeSlow]
    if os.path.exists(f'{mainDir}/fast/azimuth.offsets.fast') and \
            os.path.exists(f'{mainDir}/fast/range.offsets.fast'):
        dataTakeFast = line.replace('azimuth.offsets',
                                    'fast/azimuth.offsets.fast')
        dataTakeFast = dataTakeFast.replace('range.offsets',
                                            'fast/range.offsets.fast')
        new.append(dataTakeFast)
    return new


def maskMosaic(pieceName, pieceGeodat, outDir, outputMaskShape,
               baseFlags):
    #
    # make mask for region quadName
    #print(':', pieceName, pieceGeodat, outputMaskShape)
    mask = u.makeMaskFromShape(pieceGeodat, outputMaskShape)
    # print(np.sum(mask), outputMaskShape)
    # print(outputRegions[quadName][1], outputMaskShape)
    # u.writeImage(outDir+'/masked/mask-'+quadName, mask, 'u1')
    # Process velocity
    velImage = u.geoimage(geoType='velocity', verbose=False)
    piecePath = f'{outDir}/pieces/{pieceName}'
    velImage.readData(piecePath)
    # mask velocity
    if len(mask) > 1:
        imask = mask < 1
        velImage.vx[imask] = np.nan
        velImage.vy[imask] = np.nan
        velImage.v[imask] = np.nan
    # write velocity
    velImage.writeData(f'{outDir}/masked/{pieceName}')
    velImage = []
    #
    # Process errors
    errImage = u.geoimage(geoType='error', verbose=False)
    # read errors
    errImage.readData(piecePath)
    # mask errors
    if len(mask) > 1:
        errImage.ex[imask] = np.nan
        errImage.ey[imask] = np.nan
    # write errors
    errImage.writeData(f'{outDir}/masked/{pieceName}')
    errImage = []
    # dT - mask dT data created with timeOverlap
    if '-timeOverlap' in baseFlags:
        dT = u.geoimage(geoType='scalar', verbose=False)
        dT.readData(f'{piecePath}.dT')
        if len(mask) > 1:
            dT.x[imask] = np.nan
        dT.writeData(f'{outDir}/masked/{pieceName}.dT')
        dT = []


def maskMosaics(piecesFiles, piecesGeodats, outDir, outputMaskShape,
                baseFlags):
    #
    # loop through each regional mosaic (north...)
    print('Mask mosaics...')
    threads = []
    for pieceName, pieceGeodat in zip(piecesFiles, piecesGeodats):
        threads.append(threading.Thread(target=maskMosaic,
                                        args=[pieceName,
                                              pieceGeodat,
                                              outDir,
                                              outputMaskShape,
                                              baseFlags]))
    u.runMyThreads(threads, 10, '', quiet=True)
    return

# ---------------------------------------------------------------------------
# call commands to cat mosaics together
# ---------------------------------------------------------------------------


def interpMosiacs(piecesFiles, piecesGeodats, outDir, interpSize, baseFlags,
                  stdouts, stderrs):
    '''
    Call intfloat to fill holes
    '''
    components = ['.vx', '.vy']
    if '-timeOverlap' in baseFlags:
        components.append('.dT')
    # loop through each regional mosaic (north...)
    print('Interpolating Mosaics...')
    threads = []
    for pieceName, pieceGeodat, stdout, stderr in zip(piecesFiles,
                                                      piecesGeodats,
                                                      stdouts, stderrs):
        # loop through each component (.vx....)
        sx, sy = pieceGeodat.sizeInPixels()
        #
        # loop over components to interp
        for component in components:
            # intfloat -thresh $intThresh1 -wdist -nr $NR -na $NA -ratThresh 1.
            command = f'cd {outDir}; intfloat  -thresh {interpSize} ' \
                f' -wdist -nr {sx}  -na {sy} -ratThresh 1. ' \
                f'masked/{pieceName}{component} > ' \
                f'interp/{pieceName}{component}'
            #
            threads.append(threading.Thread(target=call, args=[command],
                                            kwargs={'shell': True,
                                                    'executable': '/bin/csh',
                                                    'stdout': stdout,
                                                    'stderr': stderr}))

            # call(command, shell=True, executable='/bin/csh',
            #     stdout=stdout, stderr=stderr)
            pieceGeodat.writeGeodat(
                f'{outDir}/interp/{pieceName}{component}.geodat')
    u.runMyThreads(threads, 10, '', quiet=True)
    return

# ---------------------------------------------------------------------------
# call commands to cat mosaics together
# ---------------------------------------------------------------------------


def mkQDirs(outDir, outputMaskTemplate):
    if not os.path.isdir(outDir):
        # print('creating outDir '+ outDir)
        os.mkdir(outDir)
    # director to log output to
    if not os.path.isdir(outDir+'/io'):
        # print(outDir+'/io')
        os.mkdir(outDir+'/io')
    # directory for individual submosaics
    if not os.path.isdir(outDir+'/pieces'):
        # print(outDir+'/pieces')
        os.mkdir(outDir+'/pieces')
    # directory for mosaics assembled from pieces
    # if not os.path.isdir(outDir+'/intermediate'):
    #     # print(outDir+'/intermediate')
    #     os.mkdir(outDir+'/intermediate')
    # directory for masked versions of mosaics
    if not os.path.isdir(outDir+'/masked'):
        # print(outDir+'/masked')
        os.mkdir(outDir+'/masked')
    # directory for interpolated products
    if not os.path.isdir(outDir+'/interp'):
        # print(outDir+'/interp')
        os.mkdir(outDir+'/interp')
    # directory for tif products
    if not os.path.isdir(outDir+'/tiff'):
        # print(outDir+'/tiff')
        os.mkdir(outDir+'/tiff')
    # directory for shp used to mask
    if not os.path.isdir(outDir+'/shp'):
        print(outDir+'/shp')
        os.mkdir(outDir+'/shp')
    #
    # If no shape file specified, copy a template that takes out Ellesmere
    shpFileBase = f'{outDir}/shp/{outDir}'
    if outputMaskTemplate is not None:
        if not os.path.exists(outputMaskTemplate):
            u.myerror(
                f'mkQdirs: nonexistent mask template ({outputMaskTemplate})')
    #
    # Modify here for other regions.
    if outputMaskTemplate is not None:
        if not os.path.exists(f'{shpFileBase}.shp'):
            # get the files for that "shapefile"
            print(outputMaskTemplate)
            shpName = outputMaskTemplate.replace(".shp", ".*")
            print(shpName)
            shpFiles = sorted(glob.glob(shpName))
            # print(shpFiles)
            for shpFile in shpFiles:
                shutil.copyfile(shpFile,
                                f'{shpFileBase}.{shpFile.split(".")[-1]}')
    outputMaskShape = f'{shpFileBase}.shp'
    return outputMaskShape


def writeQTiff(pieceFile, outDir, baseFlags, epsg, wktFile):
    velImage = u.geoimage(geoType='velocity', verbose=False)
    velRootIn = f'{outDir}/interp/{pieceFile}'
    velRootOut = f'{outDir}/tiff/{pieceFile}'
    print(f'Reading {velRootIn}')
    velImage.readData(velRootIn)
    velImage.writeMyTiff(velRootOut, epsg=epsg, wktFile=wktFile,
                         computeStats=False)
    # Errors
    errImage = u.geoimage(geoType='error', verbose=False)
    errImage.readData(f'{outDir}/masked/{pieceFile}')
    errImage.writeMyTiff(velRootOut, epsg=epsg, wktFile=wktFile,
                         computeStats=False)
    # dT if applicable
    if '-timeOverlap' in baseFlags:
        dT = u.geoimage(geoType='scalar', verbose=False)
        dT.readData(f'{velRootIn}.dT')
        dT.writeMyTiff(f'{velRootOut}.dT', epsg=epsg,
                       wktFile=wktFile, noDataDefault=-2.0e9,
                       computeStats=False)


def writeQTiffs(piecesFiles, outDir, baseFlags):
    epsg = regionDefs.epsg()
    wktFile = regionDefs.wktFile()
    #
    # Loop over quadrants and write tiffs
    print('Writing tiffs....')
    #
    threads = []
    for pieceFile in piecesFiles:
        threads.append(threading.Thread(target=writeQTiff,
                                        args=[pieceFile,
                                              outDir,
                                              baseFlags,
                                              epsg,
                                              wktFile]))
    u.runMyThreads(threads, 10, '', quiet=True)
    return

# ---------------------------------------------------------------------------
# call to produce vrt and pyramids
# ---------------------------------------------------------------------------


def makeQVRTSs(piecesFiles, outDir, template):
    print('Making VRTS...')
    suffixes = ['.vx', '.vy', '.v', '.ex', '.ey']
    if '-timeOverlap' in template['baseFlags']:
        suffixes.append('.dT')
    noData = {'.vx': -2.e9, '.vy': -2.0e9, '.v': -1, '.ex': -1, '.ey': -1,
              '.dT': -2.0e9}
    print('Creating vrts....')
    for suffix in suffixes:
        command = f'cd  {outDir}/tiff ; ' \
            f'gdalbuildvrt -te {template["xmin"]} {template["ymin"]} ' \
            f'{template["xmax"]} {template["ymax"]} ' \
            f'-tr {template["dx"]} {template["dy"]} -vrtnodata' \
            f' {noData[suffix]} {outDir}{suffix}.vrt *{suffix}.tif'
        print(command)
        call(command, shell=True)  # executable='/bin/csh')

# ---------------------------------------------------------------------------
# read landsat list file and filter by date
# ---------------------------------------------------------------------------


def processLsList(lsFile, firstDate, lastDate, outDir, noReprocess=False):
    try:
        inputFile = open(lsFile, 'r')
    except Exception:
        u.myerror('could not open landsat input file: '+lsFile)
    windowWidth = lastDate - firstDate
    centerWindowDate = firstDate + windowWidth * 0.5
    lsCulledFile = f'lsData.{firstDate.strftime("%Y-%m-%d")}.' \
        f'{lastDate.strftime("%Y-%m-%d")}'
    if noReprocess:
        return lsCulledFile
    #
    outputFile = open(outDir+'/'+lsCulledFile, 'w')
    #
    for line in inputFile:
        # use line length to skip over non data stuff
        if '&' in line:
            break
        elif ';' in line:
            continue
        datFile = line.split()[0]+'.dat'
        lsdat = u.lsdat(noProjection=True)
        lsdat.readLSdat(datFile)
        #
        # compute wSkew - use to discard scenes where the center date is >
        # one windowWidth from the center of the window.
        obsWindow = lsdat.date2 - lsdat.date1
        centerObsDate = lsdat.date1 + obsWindow*.5
        percentCover = windowWidth/obsWindow
        wSkew = abs(centerWindowDate - centerObsDate)
        # modified - 10/9/2019 to add 0.5 before windowWidth to avoid
        # data to far out of bounds (see notes)
        if (lsdat.date1 < lastDate and lsdat.date2 > firstDate and
                percentCover >= 0.6 and wSkew <= windowWidth * 0.5):
            print(line, file=outputFile, end='')
    print('&', file=outputFile)
    inputFile.close()
    outputFile.close()
    return lsCulledFile

# ---------------------------------------------------------------------------
# call to produce shape fles (meta/LS8... and meta/SAR...) with frame info
# ---------------------------------------------------------------------------


def makeShapeOutputs(outDir, firstDate, lastDate):
    print('Making Shape Files...')
    #
    # setup radar shapes
    inputFile = sorted(glob.glob(f'{outDir}/pieces/inputFile.*'))[0]
    if not os.path.isdir(outDir+'/meta'):
        os.mkdir(outDir+'/meta')
    command = 'makeimageshapefile.py -velocity -firstdate ' \
        f'{firstDate.strftime("%Y:%m:%d")}' \
        f' -lastdate {lastDate.strftime("%Y:%m:%d")} {inputFile}' \
        f' {outDir}/meta/SAR.{firstDate.strftime("%Y-%m-%d")}.' \
        f'{lastDate.strftime("%Y-%m-%d")}'
    print(command)
    call(command, shell=True)  # , executable='/bin/csh')
    #
    # setup landsat shapes
    lsFile = f'{outDir}/lsData.{firstDate.strftime("%Y-%m-%d")}.' \
        f'{lastDate.strftime("%Y-%m-%d")}'
    if os.path.isfile(lsFile):
        command = f'pathrowshp.py ' \
            f'-firstdate={firstDate.strftime("%Y:%m:%d")}' \
            f' -lastdate={lastDate.strftime("%Y:%m:%d")}' \
            f' {lsFile} {outDir}/meta/LS.' \
            f'{firstDate.strftime("%Y-%m-%d")}.' \
            f'{lastDate.strftime("%Y-%m-%d")}'
        print(command)
        call(command, shell=True)  # , executable='/bin/csh')
    # done
    return


# ---------------------------------------------------------------------------
# make the preview
# ---------------------------------------------------------------------------


def getWKT(template):
    if template['wktFile'] is not None:
        wkt = readWKT(template['wktFile'])
    elif template['epsg'] is not None:
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(template['epsg'])
        wkt = sr.ExportToWkt()
    return wkt


def writeRGBTiff(tiffFile, rgb, template):
    '''
    Write an RGB tiff for preview
    '''
    options = ['BIGTIFF=NO', 'COMPRESS=LZW', 'GEOTIFF_VERSION=1.1',
               'RESAMPLING=AVERAGE', 'PREDICTOR=YES']
    # Create mem version
    driver = gdal.GetDriverByName('MEM')
    dst_ds = driver.Create('', template['sx'], template['sy'], 3,
                           gdal.GDT_Byte)
    dst_ds.SetGeoTransform((template['xmin'], template['dx'], 0,
                            template['ymax'], 0, -template['dy']))
    # set projection
    wkt = getWKT(template)
    dst_ds.SetProjection(wkt)
    # driver specific stuff
    #
    dst_ds.FlushCache()
    # Create copy for the COG.
    for index in range(0, 3):
        band = np.squeeze(rgb[:, :, index]*255).astype(np.uint8)
        dst_ds.GetRasterBand(index+1).WriteArray(band)
    # Now make a copy
    driver = gdal.GetDriverByName('COG')
    dst_ds2 = driver.CreateCopy(tiffFile, dst_ds, options=options)
    dst_ds2.FlushCache()
    dst_ds.FlushCache()
    del(dst_ds)
    del(dst_ds2)
    gc.collect()


def addLabel(rgb, template):
    '''
    Label image with ESA credit
    '''
    labels = ['NASA MEaSUREs GrIMP Ice Velocity Map',
              'Produced using one or more of the following:',
              'Copernicus Sentinel 1 data processed by ESA',
              'TerraSAR-X/TanDEM-X data processed by DLR',
              'Landsat 8 data processed by USGS']
    # Remove GrIMP for non Greenland products
    if template['epsg'] != 3413:
        labels = labels[1:]
    if 'xLabel' not in template or 'yLabel' not in template:
        return
    #
    font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans.ttf", 75)
    dy = 80
    xl, yl = template['xLabel'], rgb.shape[0] - template['yLabel']
    #
    #  Label each band
    for bnum in range(0, 3):
        img = Image.fromarray(np.squeeze(rgb[:, :, bnum]))
        draw = ImageDraw.Draw(img)
        # Burn in each line
        for textLabel, i, in zip(labels, range(0, len(labels))):
            location = [xl, yl + i * dy]
            draw.text(location, textLabel, fill=0, font=font)
            rgb[:, :, bnum] = np.array(img)
    return rgb


def writeJPG(browseTiff):
    ''' Make jpg quicklook browse '''
    jpgFile = browseTiff.replace('.tif', '.jpg')
    jpgRes = 500
    command = f'gdal_translate -co "QUALITY=99" -scale -of JPEG -r ' \
        f'average -tr {jpgRes} {jpgRes} {browseTiff} {jpgFile}'
    call(command, shell=True)  # , executable='/bin/csh')


def processFinalPreview(outDir, template):
    '''
    Create the final preview
    '''
    #
    finalSpeed = glob.glob(f'{outDir}/release/*vv*.tif')[0]
    finalBrowse = finalSpeed.replace('vv', 'browse')
    # Create a long scale rgb image the same way idl code did
    vmin, vmax = 1, 3000
    v = u.geoimage(geoType='scalar')
    v.readData(finalSpeed, tiff=True, epsg=template['epsg'],
               wktFile=template['wktFile'])
    speed = v.x
    background = speed < 0
    print('Speed', speed.shape)
    # Use fixed value
    value = np.full(speed.shape, 1).astype('float32')
    # Use speed dependent saturation
    saturation = np.clip((speed/125 + .5)/1.5, 0, 1).astype('float32')
    saturation[background] = 0
    # Compute speed dependent
    hue = np.log10(np.clip(speed, vmin, vmax)) / \
        (np.log10(vmax) - np.log(vmin)).astype('float32')
    # order axes so rgb bands indexed last for imshow
    hsv = np.moveaxis(
        np.array([hue, saturation, value]), 0, 2).astype('float32')
    # Flip for tiff output
    rgb = np.flipud(colors.hsv_to_rgb(hsv))
    print('RGB', rgb.shape)
    del(hsv)
    gc.collect()
    # Label image
    rgb = addLabel(rgb, template)
    # Write result to release dir
    writeRGBTiff(finalBrowse, rgb, template)
    del(rgb)
    gc.collect
    #
    writeJPG(finalBrowse)


def prodType(d1, d2):
    '''
    Use duration to determine monthly, quarterly.
    Expand types as needed in future
    '''
    dT = (d2-d1).days
    if dT <= 12:
        return 's1cycle'
    elif dT >= 27 and dT <= 31:
        return 'Monthly'
    elif dT >= 86 and dT <= 93:
        return 'Quarterly'
    elif dT >= 363 and dT <= 367:
        return 'Annual'
    elif dT > 367:
        return 'multiYear'
    else:
        u.myerror(f'dates {d1.strftime("%Y-%m-%d")} and '
                  f'{d2.strftime("%Y-%m-%d")} are not monthly or quarterly '
                  f' products, fix or update code for new type')


def finalName(vrt, flags, pType):
    '''
    Create final name from intermediate vrt name
    '''
    # dates
    d1 = flags['firstdate']
    d2 = flags['lastdate']
    # compute prefix
    myPrefix = prodType(d1, d2)
    # vx, vy...
    if pType == 'v':
        pType = 'vv'
    print(currentVersions)
    print(myPrefix)
    currentVersion = currentVersions[myPrefix]
    subVersion = subVersions[myPrefix]
    name = f'{regionID}_vel_mosaic_{myPrefix}_{d1.strftime("%d%b%y")}_' \
        f'{d2.strftime("%d%b%y")}_{pType}_v{currentVersion:02d}.' \
        f'{subVersion}.tif'
    return name


def processFinalVRTs(outDir, releaseDir, flags):
    ''' Use the vrts for the individual pieces to build the finall tiffs'''
    u.pushd(outDir+'/tiff')
    vrts = sorted(glob.glob('*.vrt'))
    #
    threads = []
    for vrt in vrts:
        pType = vrt.split('.')[-2]
        pName = finalName(vrt, flags, pType)
        # cloud optimize - use rio because it works with the
        # vrt (if single, would be automatic)
        command = 'gdal_translate -of COG -co COMPRESS=DEFLATE ' \
            '-co RESAMPLING=AVERAGE -co OVERVIEWS=IGNORE_EXISTING ' \
            '-co GEOTIFF_VERSION=1.1 -co BIGTIFF=NO -stats ' \
            f'{vrt} ../release/{pName}'
        print(command)
        threads.append(threading.Thread(target=call, args=[command],
                                        kwargs={'shell': True,
                                                'executable': '/bin/csh'}))
        #call(command, shell=True, executable='/bin/csh')
    u.runMyThreads(threads, 8, '', quiet=True)
    u.popd()


def prepareReleaseDir(releaseDir):
    '''if directory doesn't exist, then create'''
    if not os.path.isdir(releaseDir):
        os.mkdir(releaseDir)
    # otherwise clean out any old tif files
    else:
        u.pushd(releaseDir)
        tiffs = []
        for pattern in ['*.tif', '*.ovr', '*.vrt']:
            tiffs += glob.glob(pattern)
        for tiff in tiffs:
            os.remove(tiff)
        u.popd()


def shapeName(shapeFile):
    ''' Create final name from intermediate vrt name '''
    pieces = shapeFile.split('.')
    # dates
    satType = pieces[0]
    d1 = datetime.strptime(pieces[1], '%Y-%m-%d')
    d2 = datetime.strptime(pieces[2], '%Y-%m-%d')
    myPrefix = prodType(d1, d2)
    mySuffix = pieces[3]
    # regionID = getRegionID()
    currentVersion = currentVersions[myPrefix]
    subVersion = subVersions[myPrefix]
    name = f'{regionID}_vel_mosaic_{myPrefix}_{d1.strftime("%d%b%y")}_' \
        f'{d2.strftime("%d%b%y")}_{satType}_v{currentVersion:02d}.' \
        f'{subVersion}.{mySuffix}'
    return name


def checkFullShape(shapeFile):
    ''' Check all parts of of shape file exist '''
    if shapeFile is None:
        return False
    baseName = shapeFile.replace('.shp', '')
    for suffix in ['shp', 'dbf']:
        if not os.path.exists(f'{baseName}.{suffix}'):
            return False
    return True


def containsLS(LSshp, satType):
    ''' check how many records'''
    # no file so return zero records
    if not checkFullShape(LSshp):
        return False
    # return true if any found
    shapeData = gpd.read_file(LSshp)
    for column in ['LSImage1', 'LSImage2']:
        for sat in {'LS8': ['LC8', 'LC08'], 'LS9': ['LC09']}[satType]:
            if len(shapeData.loc[shapeData[column].str.contains(sat)]) > 0:
                return True
    return False


def SatTypes(outDir, shapeFilesNames, currentDir=False):
    ''' check what sat types in LS8 (empty?) and SAR shape files '''
    if not currentDir:
        u.pushd(f'{outDir}/release')
    # any records, record LS8 present
    sTypes = {'S1A': False, 'S1B': False, 'TSX': False, 'TDX': False,
              'CSK': False, 'LS8': False, 'LS9': False, 'RADARSAT-1': False,
              'ALOS1-PALSAR': False, 'ERS-1/2': False}
    print(shapeFilesNames)
    if 'LS' in shapeFilesNames:
        for LStype in ['LS8', 'LS9']:
            if containsLS(shapeFilesNames['LS'], LStype):
                sTypes[LStype] = True
        # check file exists
        if not checkFullShape(shapeFilesNames['SAR']):
            u.mywarning(f'SatTypes: Missing shapefile {shapeFilesNames["SAR"]}'
                        f' in directory {os.getcwd()}')
            return sTypes
    # check sar shapes
    if 'SAR' in shapeFilesNames:
        current = shapefile.Reader(shapeFilesNames['SAR'])
        myDict = dict(zip([x[0] for x in current.fields],
                          range(0, len(current.fields))))
        for myRecord in current.records():
            sTypes[myRecord[myDict['SAT1']]] = True
        current = []
    if not currentDir:
        u.popd()
    return sTypes


def makeShapePremet(outDir, myShapeName,  shapeType, sTypes, firstDate,
                    lastDate):
    ''' write premet data for shapefile '''
    # change to release dir
    u.pushd(f'{outDir}/release')
    #  create premet file name
    preMetFile = myShapeName.replace('.shp', '.premet')
    preMets = []
    # Split sar/optical platforms based on shapeType
    mySats = {'LS': ['LS8', 'LS9'],
              'SAR': ['S1A', 'S1B', 'TSX', 'TDX', 'CSK']}
    validPlatforms = {}
    for mySat in mySats[shapeType]:
        validPlatforms[mySat] = sTypes[mySat]
    # print('validPlatforms:', validPlatforms)
    #
    preMets = populatePremet(preMets, firstDate, lastDate,
                             validPlatforms=validPlatforms)
    prodName = myShapeName
    # now write the file
    fp = open(preMetFile, 'w')
    # forcee file name to the beginning
    preMets = [preMets[0],
               s.premet({'ProducerGranuleId': prodName.split('.shp')[0]},
                        container='DataGranule')]+preMets[1:]
    for preMet in preMets:  # Output results
        preMet.printParams(fp=fp)
    fp.close()
    u.popd()


def makeFinalShapes(firstDate, lastDate):
    ''' Create the final shape files by copying from meta directory'''
    shapeFileNames = {'SAR': None, 'LS': None}
    for shapeType in ['SAR', 'LS']:  # Loop for each type
        # get meta file versions
        shapeFiles = sorted(glob.glob(f'../meta/{shapeType}.*.*.*'))
        # if more than 1, then copy (avoids .shp only)
        if len(shapeFiles) > 1:
            # loop on files
            for shapeFile in shapeFiles:
                # get suffix
                suffix = shapeFile.split('.')[-1]
                # valid suffix, then copy
                if suffix in ['shp', 'dbf', 'shx', 'prj']:
                    # remap name to final version and cpy
                    myShapeName = shapeName(os.path.basename(shapeFile))
                    shutil.copyfile(shapeFile, myShapeName)
                    if suffix == 'shp':  # save shape file name
                        shapeFileNames[shapeType] = myShapeName
    return shapeFileNames


def makeFinalOutput(outDir, template, flags):
    ''' Call routines to build final outputs (e.g., release) '''
    print('Making Final ouputs...')
    # setup release dir
    releaseDir = outDir+'/release'
    # u.mywarning('uncomment prepare release')
    prepareReleaseDir(releaseDir)
    #
    # switch to tiff dir and make single tiffs
    # u.mywarning('uncomment process final')
    processFinalVRTs(outDir, releaseDir, flags)
    #
    # this will assemble preview from pieces into single cog tif
    print('Making preview....')
    processFinalPreview(outDir, template)
    #
    # Make final shape
    u.pushd(releaseDir)
    # copy shape files
    shapeFileNames = makeFinalShapes(flags['firstdate'], flags['lastdate'])
    u.popd()
    return shapeFileNames


def populatePremet(preMetData, firstDate, lastDate, validPlatforms={}):
    ''' populate premet data fields '''
    currentVersion = currentVersions[prodType(firstDate, lastDate)]
    preMetData.append(s.premet({'VersionID_local': f'{currentVersion}'}))
    preMetData.append(s.premet(
        {'Begin_date': f'{firstDate.strftime("%Y-%m-%d")}'}))
    preMetData.append(s.premet(
        {'End_date': f'{lastDate.strftime("%Y-%m-%d")}'}))
    preMetData.append(s.premet({'Begin_time': '00:00:01.000'}))
    preMetData.append(s.premet({'End_time': '23:59:59.0000'}))
    CSKsens = {'CSK': 'X-SAR', 'CSKS1': 'X-SAR', 'CSKS2': 'X-SAR',
               'CSKS3': 'X-SAR', 'CSKS4': 'X-SAR'}
    sensorShort = {'S1A': 'C-SAR', 'S1B': 'C-SAR', 'TSX': 'X-SAR',
                   'TDX': 'X-SAR', 'LS8': 'OLI', 'LS9': 'OLI-2',
                   'RADARSAT-1': 'C-SAR',
                   'ERS-1/2': 'C-SAR', 'ALOS1-PALSAR': 'PALSAR'}
    sensorShort = {**CSKsens, **sensorShort}
    CSKinst = {'CSK': 'X-SAR', 'CSKS1': 'X-SAR', 'CSKS2': 'X-SAR',
               'CSKS3': 'X-SAR', 'CSKS4': 'X-SAR'}
    instShort = {'S1A': 'C-SAR', 'S1B': 'C-SAR', 'TSX': 'X-SAR',
                 'TDX': 'X-SAR', 'LS8': 'OLI', 'LS9': 'OLI-2',
                 'RADARSAT-1': 'C-SAR',
                 'ERS-1/2': 'C-SAR', 'ALOS1-PALSAR': 'PALSAR'}
    instShort = {**CSKinst, **instShort}
    cskIns = {'CSKS1': 'COSMO-SKYMED 1', 'CSKS2': 'X-SAR',
              'CSKS3': 'COSMO-SKYMED 3', 'CSKS4': 'COSMO-SKYMED 4'}
    instrument = {'S1A': 'Sentinel-1A', 'S1B': 'Sentinel-1B', 'TSX': 'TSX',
                  'TDX': 'TDX', 'CSK': 'COSMO-SKYMED',
                  'LS8': 'LANDSAT-8', 'LS9': 'LANDSAT-9',
                  'RADARSAT-1': 'RADARSAT-1', 'ERS-1/2': 'ERS-1-2',
                  'ALOS1-PALSAR': 'ALOS1-PALSAR'}
    instrument = {**instrument, **cskIns}
    for platform in validPlatforms.keys():
        # print('VP Key', platform)
        if validPlatforms[platform]:
            myPlatform = {'AssociatedPlatformShortname': instrument[platform],
                          'AssociatedInstrumentShortname': instShort[platform],
                          'AssociatedSensorShortname': sensorShort[platform]}
            preMetData.append(
                s.premet(myPlatform,
                         container='AssociatedPlatformInstrumentSensor'))
    #
    pType = prodType(firstDate, lastDate)
    idString = 'MEaSUREs Greenland {pType} Ice Sheet Velocity '
    if pType != 's1cycle':
        idString += 'Mosaics from SAR and Landsat'
    else:
        idString += 'Mosaics from SAR'
    #
    dataSetID = {'DataSetId': idString}
    # force this to be first
    preMetData = [s.premet(dataSetID, container='Collection')] + preMetData
    return preMetData


def makePremets(preMets, outDir, sTypes):
    ''' write the premet file '''
    u.pushd(f'{outDir}/release')
    suffixes = ['vx']  # [],'vy','vv','ex','ey','dT','browse']
    # generate file name using velocity as template
    preMetFiles = []
    # regionID = getRegionID()
    for suffix in suffixes:
        velX = sorted(glob.glob(f'{regionID}_vel*{suffix}*.tif'))
        if len(velX) == 0:
            u.myerror(
                f'makePremets: could not find release/{regionID}_vel*{suffix}'
                '*.tif as a template')
        velX = velX[0]
        # prodName=os.path.basename(velX)
        granuleName = f'{velX.replace(".tif","").replace("_vx","")}'
        GranuleID = s.premet({'ProducerGranuleId': granuleName},
                             container='DataGranule')
        preMetFile = f'{granuleName}.premet'
        #
        preMetsSuffix = [preMets[0], GranuleID]+preMets[1:]
        # now write the file
        print(f'Writing premet {preMetFile}')
        fp = open(preMetFile, 'w')
        for preMet in preMetsSuffix:
            preMet.printParams(fp=fp)
        fp.close()
        preMetFiles.append(preMetFile)
    u.popd()
    return preMetFiles


def makeSpatial(outDir, spatialFile, template):
    ''' write spatial file '''
    u.pushd(f'{outDir}/release')
    print(f'Making spatial file {spatialFile} in {outDir}/release')
    #
    llCorners = mosf.makeImageLLCorners(template, template)
    # assume makePremet error check
    with open(spatialFile, 'w') as fp:
        # corner points
        print(f'{llCorners["ll"]["lon"]:10.5f} {llCorners["ll"]["lat"]:10.5f}',
              file=fp)
        print(f'{llCorners["ul"]["lon"]:10.5f} {llCorners["ul"]["lat"]:10.5f}',
              file=fp)
        print(f'{llCorners["ur"]["lon"]:10.5f} {llCorners["ur"]["lat"]:10.5f}',
              file=fp)
        print(f'{llCorners["lr"]["lon"]:10.5f} {llCorners["lr"]["lat"]:10.5f}',
              file=fp)
    #    fp.close()
    u.popd()


def runSubMosaic(outPath, inputFileName, outFileName, baseFlags, firstDate,
                 lastDate, lsFile, dem, stdout, stderr):
    '''
    run sub mosaics (called as a thread).
    '''
    #
    # open outputs and check
    try:
        #
        productType = prodType(firstDate, lastDate)
        firstDateAdjusted, lastDateAdjusted = firstDate, lastDate
        # This expands limits so the phase data can included.
        # Offsets data should be filtered through the input.
        if productType == 'Annual':
            firstDateAdjusted = firstDateAdjusted + timedelta(-61)
            lastDateAdjusted = lastDateAdjusted + timedelta(121)
        # run cull command if doCull set (otherwise only apply coverage
        flags = baseFlags + ' -date1 ' + firstDateAdjusted.strftime('%m-%d-%Y')
        flags += ' -date2 ' + lastDateAdjusted.strftime('%m-%d-%Y')
        if lsFile is not None:
            flags += f' -landSat {lsFile}'
        command = f'cd {outPath}; mosaic3d -writeBlank -center {flags} ' \
            f'pieces/{inputFileName} {dem} pieces/{outFileName}'
        print(command, file=stdout)
        # , executable='/bin/csh'
        call(command, shell=True, stdout=stdout, stderr=stderr)
    except Exception:
        # if missing files, reject to the NoResult directory
        u.mywarning(f'warning: could not run {outPath}/runOff')
    return


def processSectors(template, flags, outDir, dataTakes, lsCulledFile):
    '''
    Breakup into sections
    '''
    threads, outPieces, geoPieces, stderrs, stdouts = [], [], [], [], []
    xdims, ydims, srsInfo, mTemp = mosf.sectionTemplate(template['template'])
    if regionDefs.epsg() is not None:
        wkt = regionDefs.epsg()
    else:
        try:
            with open(regionDefs.wktFile(), 'r') as fp:
                wkt = fp.readline()
        except Exception:
            u.myerror(f'Problem reading wkt file: {regionDefs.wktFile()}')
    #
    for xdim, xi in zip(xdims, range(len(xdims))):
        for ydim, yi in zip(ydims, range(len(ydims))):
            stdout = open(f'{outDir}/io/pieces.stdout.{xi:03d}.{yi:03d}', 'w')
            stdouts.append(stdout)
            stderr = open(f'{outDir}/io/pieces.stderr.{xi:03d}.{yi:03d}', 'w')
            stderrs.append(stderr)
            #
            inputSubFile = f'inputFile.{xi:03d}.{yi:03d}'
            outFileSub = f'mosaic-{xi:03d}.{yi:03d}'
            outPieces.append(outFileSub)
            if not flags['noReprocess']:
                writeInputFileQ(f'{outDir}/pieces/{inputSubFile}',
                                xdim, ydim, dataTakes, flags)
            # Setup thread for this sector
            threads.append(threading.Thread(target=runSubMosaic,
                                            args=[outDir,
                                                  inputSubFile,
                                                  outFileSub,
                                                  template['baseFlags'],
                                                  flags['firstdate'],
                                                  flags['lastdate'],
                                                  lsCulledFile,
                                                  regionDefs.dem(),
                                                  stdout,
                                                  stderr]
                                            )
                           )
            gd = u.geodat(x0=xdim['xll'] * 0.001, y0=ydim['yll'] * 0.001,
                          dx=xdim['dx'], dy=ydim['dy'],
                          xs=xdim['sx'], ys=ydim['sy'],
                          verbose=False, wkt=wkt)
            geoPieces.append(gd)

    stdout = open(f'{outDir}/io/All.stdout', 'w')
    stderr = open(f'{outDir}/io/All.stderr', 'w')

    return threads,  outPieces, geoPieces, stdouts, stderrs

# ----------------------------------------------------------------------------
#  main
# ----------------------------------------------------------------------------


def main():
    '''Set up mosaic files, makes directories if needed, produces all files,
    and can run with runAll '''
    global regionDefs
    global outputMaskTemplate
    global regionID
    #
    # get command line args
    maxThreads = 24
    interpSize = 50
    #
    print('\nTemplate')
    template, flags, regionDefs = processQArgs()
    for key in template:
        print(f'{key}: {template[key]}')
    print('\nFlags')
    for key in flags:
        print(f'{key}: {flags[key]}')
    #
    # setup global defs
    regionID = template['regionID']
    # Make output directory
    outDir = f'Vel-{flags["firstdate"].strftime("%Y-%m-%d")}.' \
        f'{flags["lastdate"].strftime("%Y-%m-%d")}'
    #
    # make subdirectories, set up output mask (generic or specific)
    outputMaskShape = mkQDirs(outDir, template['outputMask'])
    #
    # process landsat list (filter by date)
    lsCulledFile = None
    if template['lsFile'] is not None:
        lsCulledFile = processLsList(template['lsFile'], flags['firstdate'],
                                     flags['lastdate'], outDir,
                                     noReprocess=flags['noReprocess'])
    #
    # parse source list of data takes
    dataTakes = []
    if not flags['noReprocess']:
        dataTakes = processInputFileQ(template['inputFile'], flags['noTSX'])
        #
    threads, piecesFiles, piecesGeodats, stdouts, stderrs = \
        processSectors(template, flags, outDir, dataTakes, lsCulledFile)
    #
    # Run threads if noReprocess=False
    if not flags['noReprocess']:
        u.runMyThreads(threads, maxThreads, 'sector vel')
    #
    # Combine submosaics
    if True:  # Debug: set false to bypass
        #
        # Mask data
        maskMosaics(piecesFiles, piecesGeodats, outDir, outputMaskShape,
                    template['baseFlags'])
        #
        # Interpolate mosaics
        interpMosiacs(piecesFiles, piecesGeodats, outDir, interpSize,
                      template['baseFlags'], stdouts, stderrs)
        #
        # Interpolate tiffs
        writeQTiffs(piecesFiles, outDir, template['baseFlags'])
        #
        # makevrts
        makeQVRTSs(piecesFiles, outDir, template)
        #
    # makeshapes (meta files with frame locations)
    makeShapeOutputs(outDir, flags['firstdate'], flags['lastdate'])
    # make final output
    shapeFileNames = makeFinalOutput(outDir, template, flags)
    #
    sTypes = SatTypes(outDir, shapeFileNames)
    #
    preMetData = []
    preMetData = populatePremet(preMetData,
                                flags['firstdate'], flags['lastdate'],
                                validPlatforms=sTypes)
    preMetFiles = makePremets(preMetData, outDir, sTypes)
    #
    for shapeNameKey in shapeFileNames:
        # print(shapeNameKey, shapeFileNames)
        if shapeFileNames[shapeNameKey] is not None:
            makeSpatial(outDir,
                        shapeFileNames[shapeNameKey].replace('shp', 'spo'),
                        template)
            makeShapePremet(outDir, shapeFileNames[shapeNameKey], shapeNameKey,
                            sTypes, flags['firstdate'], flags['lastdate'])
    #
    for preMetFile in preMetFiles:
        if preMetFile is not None:
            print(preMetFile)
            makeSpatial(outDir, preMetFile.replace('premet', 'spo'), template)

if __name__ == "__main__":
    main()
