#!/usr/bin/env python3
import utilities as u
import sys
import os
import sarfunc as s
from subprocess import call
from datetime import datetime, timedelta
import numpy as np
import threading
import pyproj
import shutil
import shapefile

# Currently this is hardwired here since this program only works for greenland.
# It could be used for future expansion to Antarctica I have tried to use it to
# flag ice sheet specific stuff

iceSheetRegion = 'greenland'
currentVersions = {'s1cycle': 1, 'Monthly': 3, 'Quarterly': 3, 'Annual': 3,
                   'multiYear': 1}
subVersions = {'s1cycle': 0, 'Monthly': 0, 'Quarterly': 0, 'Annual': 0,
               'multiYear': 0}
defaultGreenlandDEM = s.defaultRegionDefs('greenland').dem()
defaultAntarcticDEM = '/Volumes/insar7/ian/TSXsouth/PIG/DEM/SentinelPigDEM'


def usageQ():
    #
    print('\n\t\t\t\t \033[1;34m SETUP QUADRANT MOSAICS \033[0m')
    print('\n'+'\033[1m'+'Use a template inputFile to create a mosaic made of'
          ' four quadrants for a specified date range')
    print("\n\trunquarters.py -timeOverlapOff  -noReprocess "
          "-baseFlags='-flag1...' -input=inputFile ")
    print("\t-shelfMask=shelfMask -ouputMask=outputMask.shp -lsFile=lsFile "
          "-firstdate=YYYY-MM-DD -lastdate=YYYY-MM-DD -noTSX\n ")
    print('where:\n')
    print('\t -dem\t\t\tis the full path to the DEM file '
          f'[{defaultGreenlandDEM}]')
    print('\t -baseFlags\t\tmosaic3d flags that are applied to all cases'
          ' [-offsets -rOffsets -fl 10 -noVh -3dOff]')
    print('\t -inputFile\t\tinputFile with all tracks to be considered '
          '[track-all/inputFile]')
    print('\t -lsFile\t\tfile with all the landsat  tracks to be considered '
          '[none]')
    print('\t -firstdate\t\tfirst date  of image in mosiac (required)')
    print('\t -lastdate\t\tlast date of image in mosaic (required)')
    print('\t -timeOverlapOff\tBy default timeOverlap called with mosaic3d '
          '- use flag to turn off [On]')
    print('\t -noReprocess\t\tRebuild tiffs with existing mosaic3d '
          'output [Off]')
    print('\t -noCull\t\tUse unculled offsets [Off]')
    print('\t -noTSX\t\t\texclude TSX files [include TSX]')
    print('\t -amundsen\t\tamundsen map [greenland]')
    print('\t -greenland\t\tgreenland map [greenland]')
    print('\t -shelfMask\t\tMask to pass to mosaic3d (byte mask)')
    print('\t -outputMask\t\tMask to apply to ouput to remove bad areas '
          '(shape file with regions to mask)')
    print('\nNotes:\n')
    print('\033[0m')
    sys.exit()

# -----------------------------------------------------------------------------
# Process command line args.
# -----------------------------------------------------------------------------


def processQArg(argv):
    # Setup defaults
    baseFlags = '-offsets -rOffsets -fl 10 -noVh -3dOff'
    inputFile = 'track-all/inputFile'
    lsFile = None
    noTSX = False
    firstDate = ' '
    lastDate = ' '
    outputMaskTemplate = None
    shelfMask = ''
    timeOverlapOff = False
    noReprocess = False
    noCull = False
    helpFlag = False
    # for now this is not an option, but allows for a change later
    cloudOptimize = True
    print('\n')
    global iceSheetRegion
    try:
        for myStr in argv[1:]:
            print(myStr)
            if '-help' in myStr:
                helpFlag = True
                usageQ()
            elif '-noTSX' in myStr:
                noTSX = True
            elif '-baseFlags' in myStr:
                temp = myStr.split('=')[1:]
                baseFlags = ' '.join(temp)
                print('base flags', baseFlags)
            elif '-input' in myStr:
                inputFile = myStr.split('=')[-1]
                print('inputFile = ', inputFile)
            elif '-lsFile' in myStr:
                lsFile = myStr.split('=')[-1]
                print('lsFile = ', lsFile)
            elif '-firstdate' in myStr:
                firstDate = datetime.strptime(myStr.split('=')[-1], '%Y-%m-%d')
                print('First Date = ', firstDate)
            elif '-lastdate' in myStr:
                lastDate = datetime.strptime(myStr.split('=')[-1], '%Y-%m-%d')
                lastDate += timedelta(seconds=86400-1)
                print('Last Date = ', lastDate)
            elif '-shelfMask' in myStr:
                shelfMask = myStr.split('=')[-1]
                print('shelfMask = ', shelfMask)
            elif '-outputMask' in myStr:
                outputMaskTemplate = myStr.split('=')[-1]
                print('outputMask = ',  outputMaskTemplate)
            elif '-timeOverlapOff' in myStr:
                timeOverlapOff = True
                print('timeOverlapOff = ', timeOverlapOff)
            elif '-noReprocess' in myStr:
                noReprocess = True
                print('noReprocess = ', noReprocess)
            elif '-noCull' in myStr:
                noCull = True
                print('noCull = ', noCull)
            elif '-amundsen' in myStr:
                iceSheetRegion = 'amundsen'
            elif '-greenland' in myStr:
                iceSheetRegion = 'greenland'
            elif '-taku' in myStr:
                iceSheetRegion = 'taku'
            else:
                print(f'problem with {myStr} ')
                usageQ()
        if not timeOverlapOff:
            baseFlags += ' -timeOverlap '
        print('baseFlags = ', baseFlags)
        if type(firstDate) is str or type(lastDate) is str:
            u.myerror('no first and/or last date specficied')
        #
        print(iceSheetRegion)
        if helpFlag:
            usageQ()
        return baseFlags, inputFile, lsFile, firstDate, lastDate, noTSX,\
            outputMaskTemplate, shelfMask, noReprocess, cloudOptimize, noCull
    except Exception:
        print(myStr)
        usageQ()

# ----------------------------------------------------------------------------
# define output windows for each projection -
# ----------------------------------------------------------------------------


def readWKT(wktFile):
    ''' get wkt from a file '''
    with open(wktFile, 'r') as fp:
        return fp.readline()


def projWin(region):
    # add other regions here
    wktFile = s.defaultRegionDefs(region).wktFile()
    if wktFile is not None:
        epsg = None
        projStr = readWKT(wktFile)
    else:
        epsg = s.defaultRegionDefs(region).epsg()
        projStr = f'EPSG:{epsg}'
    projDict = {'greenland': {'xmin': -659100.0, 'xmax': 857900.0,
                              'ymin': -3379100.0, 'ymax': -639100.0,
                              'dx': 200, 'dy': 200, 'epsg': epsg,
                              'proj': projStr},
                'amundsen': {'xmin': -2000125.0, 'xmax': -1020125.0,
                             'ymin': -1250125.0, 'ymax': 289875.0,
                             'dx': 250, 'dy': 250, 'epsg': epsg,
                             'proj': projStr},
                'taku': {'xmin': -147100.0, 'xmax': 163900.,
                         'ymin': -3431100.0, 'ymax': -3208100.0,
                         'dx': 200, 'dy': 200, 'epsg': epsg, 'proj': projStr},
                'amundsenOld': {'xmin': -1904125.0, 'xmax': -1084125.0,
                                'ymin': -1200125.0, 'ymax': 199875.0,
                                'dx': 250, 'dy': 250, 'epsg': epsg,
                                'proj': projStr}}
    # use try to avoid invalid regions
    try:
        myProj = projDict[region]
    except Exception:
        u.myerror(f'invalid region={region} in projWin ')
    # compute Corners
    # llproj = pyproj.Proj('+init=EPSG:4326')
    # xyproj = pyproj.Proj(f'+init=EPSG:{myProj["epsg"]}')
    xyllXform = pyproj.Transformer.from_crs(projStr, 'epsg:4326')
    latll, lonll = xyllXform.transform(myProj['xmin'], myProj['ymin'])
    latur, lonur = xyllXform.transform(myProj['xmax'], myProj['ymax'])
    latlr, lonlr = xyllXform.transform(myProj['xmax'], myProj['ymin'])
    latul, lonul = xyllXform.transform(myProj['xmin'], myProj['ymax'])
    llCorners = {'ll': {'lat': latll, 'lon': lonll},
                 'lr': {'lat': latlr, 'lon': lonlr},
                 'ur': {'lat': latur, 'lon': lonur},
                 'ul': {'lat': latul, 'lon': lonul}}
    print(llCorners)
    # 'amundsen': ' -1904125.0 199875.0 -1084125.0 -1000125.0  -tr 250 250 ' }
    return myProj, llCorners

# -----------------------------------------------------------------------------
# set up predefined quadrants
# -----------------------------------------------------------------------------


def setupQuadrants():
    ''' For now hardwire regions as a dictonary '''
    if iceSheetRegion == 'greenland':
        names = ['north', 'northeast', 'northwest', 'south']
        quadrants = {'north': [-626., -1320., 1235., 625.0, .2, .2],
                     'northeast': [110., -2466., 740., 1146., .2, .2],
                     'northwest': [-626., -2466., 736., 1146., .2, .2],
                     'south': [-410., -3356., 915., 890., .2, .2]}
    elif iceSheetRegion == 'amundsen':
        names = ['amundsen']
        quadrants = {'amundsen': [-2000, -1250, 980, 1540, 0.25, 0.25]}
    elif iceSheetRegion == 'taku':
        names = ['taku']
        myProj, _ = projWin(iceSheetRegion)
        x0 = (myProj['xmin'] + myProj['dx']*0.5)/1000.
        y0 = (myProj['ymin'] + myProj['dy']*0.5)/1000.
        sx = (myProj['xmax'] - myProj['xmin'])/1000.
        sy = (myProj['ymax'] - myProj['ymin'])/1000.
        quadrants = {'taku': [x0, y0, sx, sy, myProj['dx']/1000.,
                              myProj['dy']/1000.]}
        # print(quadrants)
        # u.myerror('stop')
    else:
        u.myerror('Other regions not yet implemented '+iceSheetRegion)
    return names, quadrants

# -----------------------------------------------------------------------------
# Read in input file, save as list with [date, line] entries.
# -----------------------------------------------------------------------------


def processInputFileQ(filename, noTSX):
    try:
        inputFile = open(filename, 'r')
    except Exception:
        u.myerror('could not open SAR input file: '+filename)
    dataTakes = []
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
                fgeo = open(geo, 'r')
                # read geodat to extract date
                for gline in fgeo:
                    if 'Image date' in gline:
                        gline = gline.split(':')[1]
                        gline = gline.strip('\n').replace(' ', '')
                        mydate = datetime.strptime(gline, '%d%b%Y')
                        break
                fgeo.close()
                line += ' '+crossFlag
                dataTakes.append([mydate, line])
    inputFile.close()
    return dataTakes

# -----------------------------------------------------------------------------
# For a given overalll mosaics limits, break into sub regions.
# -----------------------------------------------------------------------------


def findLimitsQ(limits):
    if limits[4] != limits[5]:
        u.myerror(' X and Y resolutions differ:', limits)
    resolution = limits[4]
    #
    sX = int(limits[2]/resolution)
    sY = int(limits[3]/resolution)
    dxkm = resolution
    dykm = resolution
    y0orig = limits[1]
    x0 = limits[0]
    gd = u.geodat(x0, y0orig, sX, sY, resolution*1000, resolution*1000,
                  verbose=False)
    # maximum number of points in a file
    maxPoints = 7000000
    if iceSheetRegion == 'amundsen':
        maxPoints = int(maxPoints/4)
    # this will force more boxes for 0.5 km and greater to run faster
    if resolution > 0.49:
        maxPoints = maxPoints/4
    # force it to give integer size
    maxLines = int(maxPoints/sX)
    # iteration to adjust size so last case is not too small (>.25)
    done = False
    while not done:
        frac = float(sY)/maxLines - int(sY/maxLines)
        if frac > 0.25:
            done = True
        else:
            maxLines -= 100
    line = 0
    lastline = 0
    limits = []
    while line < sY:
        deltaLines = int(min(sY-line, maxLines))
        line += deltaLines
        y0 = y0orig+lastline*dykm
        lastline = line
        limits.append([x0, y0, sX*dxkm, deltaLines*dykm, dxkm, dykm])
    #
    return limits, gd

#
# -----------------------------------------------------------------------------
# Write an input file filtered for date range.
# ----------------------------------------------------------------------------


def writeInputFileQ(infile, limit, dataTakes, date1, date2, noCull):
    #
    fOut = open(infile, 'w')
    print('; automatically generated input file from setupannual mosaics-- ',
          file=fOut)
    # print limit
    print('{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}'.format(*limit),
          file=fOut)
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
    if date1.month not in [6, 7, 8] or date2.month not in [7, 8, 9]:
        # if more than a 15 day estimate, apply screen
        if (date2-date1).days > 15:
            screenTSX = True
            print(screenTSX)
    productType = prodType(date1, date2)
    # now loop through data takes, keeping ones in range
    for dataTake in dataTakes:
        firstDate = dataTake[0]
        tmp = dataTake[1].split()
        secondDate = firstDate + timedelta(float(tmp[3]))
        midDate = firstDate + timedelta(float(tmp[3])) * 0.5
        #
        # if this a TSX case and screenTSX on
        # Added only images with center for 6/12 day pairs.
        if ('TSX' in dataTake[1] and screenTSX) or productType == 's1cycle':
            # then check midDate in the range, otherwise
            # use continue to skip this iteration
            if midDate <= date1 or midDate >= date2:
                continue
        # pass dates with any amount of overlap
        if secondDate >= date1 and firstDate <= date2:
            if noCull:
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
    return


def changeToNoCull(line):
    mainDir = os.path.dirname(line.split()[1])
    az = u.dols(f'ls {mainDir}/*.interp.da')[0]
    azname = os.path.basename(az)
    rg = u.dols(f'ls {mainDir}/*.interp.dr')[0]
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
# ----------------------------------------------------------------------------
# run sub mosaics (called as a thread)
# ----------------------------------------------------------------------------


def runSubMosaic(outPath, inputFileName, outFileName, baseFlags, firstDate,
                 lastDate, lsFile, dem, stdout, stderr):
    #
    # open outputs and check
    try:
        #
        # run cull command if doCull set (otherwise only apply coverage
        flags = baseFlags + ' -date1 ' + firstDate.strftime('%m-%d-%Y')
        flags += ' -date2 ' + lastDate.strftime('%m-%d-%Y')
        if lsFile is not None:
            flags += ' -landSat '+lsFile
        command = f'cd {outPath}; mosaic3d -writeBlank -center {flags} ' \
            f'pieces/{inputFileName} {dem} pieces/{outFileName}'
        print(command)
        print(command, file=stdout)
        call(command, shell=True, executable='/bin/csh', stdout=stdout,
             stderr=stderr)
    except Exception:
        # if missing files, reject to the NoResult directory
        u.mywarning('warning: could not run ' + outPath + '/runOff')
    stdout.close()
    stderr.close()
    return

# ----------------------------------------------------------------------------
# call commands to cat mosaics together
# ----------------------------------------------------------------------------


def assembleMosiacs(quadNames, outputRegions, outDir, baseFlags):
    #
    # Assemble final mosaic
    print('Assemble mosaics...')
    components = ['.vx', '.vy', '.ex', '.ey']
    if '-timeOverlap' in baseFlags:
        components.append('.dT')
    print(components)
    # loop through each regional mosaic (north...)
    for quadName in quadNames:
        # loop through each component (.vx....)
        for component in components:
            # build command to cat data
            command = 'cd '+outDir + ' ; cat '
            for region in outputRegions[quadName][0]:
                command += f' pieces/{region}{component}'
            command += f' > intermediate/mosaic-{quadName}{component}'
            call(command, shell=True, executable='/bin/csh',
                 stdout=outputRegions[quadName][2],
                 stderr=outputRegions[quadName][3])
            # write geodat
            outputRegions[quadName][1].writeGeodat(
                outDir+'/intermediate/' + 'mosaic-' + quadName+component +
                '.geodat')
    return

# ---------------------------------------------------------------------------
# call commands to cat mosaics together
# ---------------------------------------------------------------------------


def maskMosaics(quadNames, outputRegions, outDir, outputMaskShape, baseFlags):
    #
    # loop through each regional mosaic (north...)
    print('Mask mosaics...')
    for quadName in quadNames:
        #
        # make mask for region quadName
        print(outputRegions[quadName][1], outputMaskShape)
        mask = u.makeMaskFromShape(outputRegions[quadName][1], outputMaskShape)
        # print(np.sum(mask), outputMaskShape)
        # print(outputRegions[quadName][1], outputMaskShape)
        # u.writeImage(outDir+'/masked/mask-'+quadName, mask, 'u1')
        # Process velocity
        velImage = u.geoimage(geoType='velocity', verbose=False)
        print(outDir+'/intermediate/mosaic-'+quadName)
        velImage.readData(outDir+'/intermediate/mosaic-'+quadName)
        # mask velocity
        if len(mask) > 1:
            imask = mask < 1
            velImage.vx[imask] = np.nan
            velImage.vy[imask] = np.nan
            velImage.v[imask] = np.nan
        # write velocity
        velImage.writeData(outDir+'/masked/mosaic-'+quadName)
        velImage = []
        #
        # Process errors
        errImage = u.geoimage(geoType='error', verbose=False)
        # read errors
        errImage.readData(outDir+'/intermediate/mosaic-'+quadName)
        # mask errors
        if len(mask) > 1:
            errImage.ex[imask] = np.nan
            errImage.ey[imask] = np.nan
        # write errors
        errImage.writeData(outDir+'/masked/mosaic-'+quadName)
        errImage = []
        # dT - mask dT data created with timeOverlap
        if '-timeOverlap' in baseFlags:
            dT = u.geoimage(geoType='scalar', verbose=False)
            dT.readData(outDir+'/intermediate/mosaic-'+quadName+'.dT')
            if len(mask) > 1:
                dT.x[imask] = np.nan
            dT.writeData(outDir+'/masked/mosaic-'+quadName+'.dT')
            dT = []
    # u.myerror('mask stop')
    return

# ---------------------------------------------------------------------------
# call commands to cat mosaics together
# ---------------------------------------------------------------------------


def interpMosiacs(quadNames, outputRegions, outDir, interpSize, baseFlags):
    components = ['.vx', '.vy']
    if '-timeOverlap' in baseFlags:
        components.append('.dT')
    # loop through each regional mosaic (north...)
    print('Interpolating Mosaics...')
    for quadName in quadNames:
        # loop through each component (.vx....)
        sx, sy = outputRegions[quadName][1].sizeInPixels()
        #
        # loop over components to interp
        for component in components:
            # build command to cat data
            # intfloat -thresh $intThresh1 -wdist -nr $NR -na $NA -ratThresh 1.
            command = f'cd {outDir}; intfloat  -thresh {interpSize} ' \
                f' -wdist -nr {sx}  -na {sy} -ratThresh 1. ' \
                f'masked/mosaic-{quadName}{component} > ' \
                f'interp/mosaic-{quadName}{component}'
            #  print(command)
            #
            call(command, shell=True, executable='/bin/csh',
                 stdout=outputRegions[quadName][2],
                 stderr=outputRegions[quadName][3])
            outputRegions[quadName][1].writeGeodat(
                f'{outDir}/interp/mosaic-{quadName}{component}.geodat')
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
    if not os.path.isdir(outDir+'/intermediate'):
        # print(outDir+'/intermediate')
        os.mkdir(outDir+'/intermediate')
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
    if outputMaskTemplate is None:
        template = True
        outputMaskTemplate = \
            f'/Users/ian/maskTemplates/{iceSheetRegion}Output.shp'
        if not os.path.exists(outputMaskTemplate):
            u.myerror(
                f'mkQdirs: nonexistent mask template ({outputMaskTemplate})')
        shpFileBase = outDir+'/shp/'+outDir
    else:
        shpFileBase = outDir+'/shp/' +\
            outputMaskTemplate.split('/')[-1].replace('.shp', '')
        template = False
    #
    # Modify here for other regions.
    if iceSheetRegion in ['greenland', 'amundsen', 'taku']:
        # if this is a template file and it doesnt exits then copy or
        # if its a custome file then copy (overwrite existing if need be)
        # print(shpFileBase, template)
        if not os.path.exists(shpFileBase+'.shp') or not template:
            # get the files for that "shapefile"
            shpFiles = u.dols(f'ls {outputMaskTemplate.replace(".shp", ".*")}')
            # print(shpFiles)
            for shpFile in shpFiles:
                command = f'cp {shpFile} {shpFileBase}.' \
                    f'{shpFile.split(".")[-1]}'
                call(command, shell=True, executable='/bin/csh')
    outputMaskShape = shpFileBase+'.shp'
    return outputMaskShape

# ---------------------------------------------------------------------------
# call commands to cat mosaics together
# ---------------------------------------------------------------------------


def writeQTiffs(quadNames, outDir, baseFlags):
    #
    # Loop over quadrants and write tiffs
    print('Writing tiffs....')
    epsg = s.defaultRegionDefs(iceSheetRegion).wktFile()
    wktFile = s.defaultRegionDefs(iceSheetRegion).wktFile()
    #
    for quadName in quadNames:
        velImage = u.geoimage(geoType='velocity', verbose=False)
        print('Reading '+outDir+'/interp/mosaic-'+quadName)
        velImage.readData(f'{outDir}/interp/mosaic-{quadName}')
        velImage.writeMyTiff(f'{outDir}/tiff/{quadName}', epsg=epsg,
                             wktFile=wktFile)
        velImage = []
        # Errors
        errImage = u.geoimage(geoType='error', verbose=False)
        errImage.readData(outDir+'/masked/mosaic-'+quadName)
        errImage.writeMyTiff(outDir+'/tiff/'+quadName, epsg=epsg,
                             wktFile=wktFile)
        errImage = []
        # dT if applicable
        if '-timeOverlap' in baseFlags:
            dT = u.geoimage(geoType='scalar', verbose=False)
            dT.readData(outDir+'/interp/mosaic-'+quadName+'.dT')
            dT.writeMyTiff(outDir+'/tiff/'+quadName+'.dT', epsg=epsg,
                           wktFile=wktFile, noDataDefault=-2.0e9)
            dT = []
    # moved to makeShapeOutputs
    # stdout,stderr=open(outDir+'/io/stderr','w'),open(outDir+'/io/stdout','w')
    # command='cd '+outDir+'/tiff ; makeimageshapefile.py -velocity
    # ../pieces/inputFile-'+quadNames[0]+'.000 '+outDir+'.sarswathes'
    # call(command,shell=True,executable='/bin/csh',stderr=stderr,stdout=stdout)
    # stderr.close()
    # stdout.close()

    return

# ---------------------------------------------------------------------------
# call to produce vrt and pyramids
# ---------------------------------------------------------------------------


def makeQVRTSs(quadNames, outDir, baseFlags):
    print('Making VRTS...')
    suffixes = ['.vx', '.vy', '.v', '.ex', '.ey']
    if '-timeOverlap' in baseFlags:
        suffixes.append('.dT')
    noData = {'.vx': -2.e9, '.vy': -2.0e9, '.v': -1, '.ex': -1, '.ey': -1,
              '.dT': -2.0e9}
    print('Creating vrts....')
    p, _ = projWin(iceSheetRegion)
    for suffix in suffixes:
        command = f'cd  {outDir}/tiff ; ' \
            f'gdalbuildvrt  -te {p["xmin"]}  {p["ymin"]} {p["xmax"]}' \
            f' {p["ymax"]} -tr {p["dx"]} {p["dy"]} -vrtnodata' \
            f' {noData[suffix]} {outDir+suffix}.vrt *{suffix}.tif'
        print(command)
        call(command, shell=True, executable='/bin/csh')

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
        obsWindow = lsdat.date2-lsdat.date1
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
    inputFile = u.dols('ls ' + outDir + '/pieces/inputFile-*')[0]
    if not os.path.isdir(outDir+'/meta'):
        os.mkdir(outDir+'/meta')
    command = 'makeimageshapefile.py -velocity -firstdate ' \
        f'{firstDate.strftime("%Y:%m:%d")}' \
        f' -lastdate {lastDate.strftime("%Y:%m:%d")} {inputFile}' \
        f' {outDir}/meta/SAR.{firstDate.strftime("%Y-%m-%d")}.' \
        f'{lastDate.strftime("%Y-%m-%d")}'
    print(command)
    call(command, shell=True, executable='/bin/csh')
    #
    # setup landsat shapes
    lsFile = f'{outDir}/lsData.{firstDate.strftime("%Y-%m-%d")}.' \
        f'{lastDate.strftime("%Y-%m-%d")}'
    if os.path.isfile(lsFile):
        command = f'pathrowshp.py ' \
            f'-firstdate={firstDate.strftime("%Y:%m:%d")}' \
            f' -lastdate={lastDate.strftime("%Y:%m:%d")}' \
            f' {lsFile} {outDir}/meta/LS8.' \
            f'{firstDate.strftime("%Y-%m-%d")}.' \
            f'{lastDate.strftime("%Y-%m-%d")}'
        print(command)
        call(command, shell=True, executable='/bin/csh')
    # done
    return

# ---------------------------------------------------------------------------
# make typical rsfig map and save in preview
# ---------------------------------------------------------------------------


def makeLogMap(outDir, quadNames):
    print('Making Log Maps...')
    #
    jpgRes = '500'
    # setup radar shapes
    if not os.path.isdir(outDir+'/preview'):
        os.mkdir(outDir+'/preview')
    u.pushd(outDir+'/preview')
    for quadName in quadNames:
        command = f"idl -e \"rsfig,maxv=3000,/nobar,file='../interp/mosaic-" \
            f"{quadName}',geohue='{quadName}.browse.tif'"
        if iceSheetRegion == 'amundsen':
            command += ',/south'
        elif iceSheetRegion == 'taku':
            command += ',/taku'
        if 'south' in quadName:
            command += ",labels=[" \
                "'NASA MEaSUREs GIMP Ice Velocity Map'," \
                "'Produced using one or more of the following:'," \
                "'Copernicus Sentinel 1 data processed by ESA'," \
                "'TerraSAR-X/TanDEM-X data processed by DLR'," \
                "'Landsat 8 data processed by USGS']"
        command += '"'
        print(command)
        call(command, shell=True, executable='/bin/csh')
    # clean up
    if os.path.exists('tmp.tif'):
        os.remove('tmp.tif')
    p, _ = projWin(iceSheetRegion)
    epsg = s.defaultRegionDefs(iceSheetRegion).epsg()
    if epsg is not None:
        myProj = f'EPSG:{epsg}'
    else:
        myProj = s.defaultRegionDefs(iceSheetRegion).wktFile()
    command = f'gdalbuildvrt -a_srs {myProj} -te {p["xmin"]}  ' \
        f'{p["ymin"]} {p["xmax"]} {p["ymax"]} -tr {p["dx"]} {p["dy"]} ' \
        f'-hidenodata -vrtnodata "255 255 255" {outDir}.browse.vrt' \
        f' *.tif'
    call(command, shell=True, executable='/bin/csh')
    command = f'gdaladdo -r average {outDir}.browse.vrt 2 4 8 16'
    call(command, shell=True, executable='/bin/csh')
    command = f'gdal_translate -co "QUALITY=99" -scale -of JPEG -r '\
        f'average -tr {jpgRes} {jpgRes} ' \
        f'-of jpeg {outDir}.browse.vrt {outDir}.browse.jpg'
    call(command, shell=True, executable='/bin/csh')
    u.popd()

# ---------------------------------------------------------------------------
# make the preview
# ---------------------------------------------------------------------------


def processFinalPreview(outDir):
    ''' Create the final preview '''
    # merge preview images
    u.pushd(outDir+'/preview')
    # get the file name
    vrtpreview = u.dols('ls *browse*.vrt')[0]
    #
    finalBrowse = finalName(vrtpreview)
    # command = f'rio cogeo create {vrtpreview} ../release/{finalBrowse}' \
    #    f' --overview-resampling average'
    #
    command = 'gdal_translate -of COG -co COMPRESS=DEFLATE ' \
        '-co RESAMPLING=AVERAGE -co OVERVIEWS=IGNORE_EXISTING ' \
        '-co GEOTIFF_VERSION=1.1 -co BIGTIFF=NO -stats ' \
        f'{vrtpreview} ../release/{finalBrowse}'
    print(command)
    call(command, shell=True, executable='/bin/csh')
    for suffix in ['jpg', 'jpg.aux.xml']:
        shutil.copyfile(f'{outDir}.browse.{suffix}',
                        f'../release/{finalBrowse.replace("tif", suffix)}')
    u.popd()

# ---------------------------------------------------------------------------
# make typical rsfig map and save in preview
# ---------------------------------------------------------------------------


def prodType(d1, d2):
    ''' Use duration to determine monthly, quarterly.
    Expand types as needed in future'''
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


def getRegionID():
    regionIDs = {'greenland': 'GL', 'amundsen': 'ASE', 'taku': 'TAKU'}
    try:
        regionID = regionIDs[iceSheetRegion]
    except Exception:
        u.error(f'finalName: Invalid region name {regionID} not in dict '
                f'{regionIDs}')
    return regionID


def finalName(vrt):
    ''' Create final name from intermediate vrt name '''
    regionID = getRegionID()
    pieces = vrt.split('.')
    # dates
    d1 = datetime.strptime(pieces[0], 'Vel-%Y-%m-%d')
    d2 = datetime.strptime(pieces[1], '%Y-%m-%d')
    # compute prefix
    myPrefix = prodType(d1, d2)
    # vx, vy...
    pType = pieces[2]
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


def processFinalVRTs(outDir, releaseDir, cloudOptimize):
    ''' Use the vrts for the individual pieces to build the finall tiffs'''
    u.pushd(outDir+'/tiff')
    vrts = u.dols('ls *.vrt')
    #
    for vrt in vrts:
        print(vrt)
        pName = finalName(vrt)
        # cloud optimize - use rio because it works with the
        # vrt (if single, would be automatic)
        if cloudOptimize:
            # command = f'rio cogeo create {vrt} ../release/{pName} ' \
            #    f'--overview-resampling average' \
            #    f' ; gdalinfo -stats  ../release/{pName}'
            command = 'gdal_translate -of COG -co COMPRESS=DEFLATE ' \
                     '-co RESAMPLING=AVERAGE -co OVERVIEWS=IGNORE_EXISTING ' \
                     '-co GEOTIFF_VERSION=1.1 -co BIGTIFF=NO -stats ' \
                     f'{vrt} ../release/{pName}'
            # command += '; gdaladdo -r average ../release/{pName} ' \
            #    '2 4 8 16 32 64 128'
        # old format
        else:
            p, _ = projWin(iceSheetRegion)
            pw = f'{p["xmin"]} {p["ymax"]} {p["xmax"]} {p["ymin"]} ' \
                f'-te {p["dx"]} {p["dy"]}'
            command = f'gdal_translate  -co "COMPRESS=LZW" -co "PREDICTOR=1"'\
                f' -co "TILED=YES" -projwin {pw+vrt} ../release/{pName}'
        print(command)
        call(command, shell=True, executable='/bin/csh')
    u.popd()
    #
    # Old format, only happens if not cloud optimized.
    # now make pyramids and vrt
    if not cloudOptimize:
        u.pushd(releaseDir)
        tiffs = u.dols('ls *.tif')
        for tiff in tiffs:
            command = f'gdalbuildvrt {tiff.replace(".tif", ".vrt")} tiff ; ' \
                f'gdaladdo  -r average {tiff.replace(".tif", ".vrt") }' \
                ' 32 16 8 4 2'
            print(command)
            call(command, shell=True, executable='/bin/csh')
        u.popd()


def makePyramidsForOldFormat(releaseDir):
    ''' Legacy routine to create pyramids for vrt '''
    vrtpreview = u.dols('ls *browse*.vrt')[0]
    command = f'gdaladdo  -r average {vrtpreview.replace(".vrt", ".tif")}' \
        ' 32 16 8 4 2'
    print(command)
    call(command, shell=True, executable='/bin/csh')


def prepareReleaseDir(releaseDir):
    '''if directory doesn't exist, then create'''
    if not os.path.isdir(releaseDir):
        os.mkdir(releaseDir)
    # otherwise clean out any old tif files
    else:
        u.pushd(releaseDir)
        tiffs = u.dols('ls *.tif *.ovr *.vrt')
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
    regionID = getRegionID()
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


def LS8Count(LS8shp):
    ''' check how many records'''
    # no file so return zero records
    if not checkFullShape(LS8shp):
        return 0
    # count recoreds
    current = shapefile.Reader(LS8shp)
    nRecords = len(current.records())
    current = []
    return nRecords


def SatTypes(outDir, shapeFilesNames, currentDir=False):
    ''' check what sat types in LS8 (empty?) and SAR shape files '''
    if not currentDir:
        u.pushd(f'{outDir}/release')
    # any records, record LS8 present
    sTypes = {'S1A': False, 'S1B': False, 'TSX': False, 'TDX': False,
              'CSK': False, 'LS8': False, 'RADARSAT-1': False,
              'ALOS1-PALSAR': False, 'ERS-1/2': False}
    if LS8Count(shapeFilesNames['LS8']) > 0:
        sTypes['LS8'] = True
    # check file exists
    if not checkFullShape(shapeFilesNames['SAR']):
        u.mywarning(f'SatTypes: Missing shapefile {shapeFilesNames["SAR"]} '
                    f'in directory {os.getcwd()}')
        return sTypes
    # check sar shapes
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
    mySats = {'LS8': ['LS8'], 'SAR': ['S1A', 'S1B', 'TSX', 'TDX', 'CSK']}
    validPlatforms = {}
    for mySat in mySats[shapeType]:
        validPlatforms[mySat] = sTypes[mySat]
    print(validPlatforms)
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
    shapeFileNames = {'SAR': None, 'LS8': None}
    for shapeType in ['SAR', 'LS8']:  # Loop for each type
        # get meta file versions
        shapeFiles = u.dols(f'ls ../meta/{shapeType}.*.*.*')
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


def makeFinalOutput(outDir, quadNames, cloudOptimize, firstDate, lastDate):
    ''' Call routines to build final outputs (e.g., release) '''
    print('Making Final ouputs...')
    # setup release dir
    releaseDir = outDir+'/release'
    prepareReleaseDir(releaseDir)
    #
    # switch to tiff dir and make single tiffs
    processFinalVRTs(outDir, releaseDir, cloudOptimize)
    #
    # this will assemble preview from pieces into single cog tif
    processFinalPreview(outDir)
    #
    # make pyramids and copy shape
    u.pushd(releaseDir)
    if not cloudOptimize:
        makePyramidsForOldFormat(releaseDir)
    # copy shape files
    shapeFileNames = makeFinalShapes(firstDate, lastDate)
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
                   'TDX': 'X-SAR', 'LS8': 'PAN', 'RADARSAT-1': 'C-SAR',
                   'ERS-1/2': 'C-SAR', 'ALOS1-PALSAR': 'PALSAR'}
    sensorShort = {**CSKsens, **sensorShort}
    CSKinst = {'CSK': 'X-SAR', 'CSKS1': 'X-SAR', 'CSKS2': 'X-SAR',
               'CSKS3': 'X-SAR', 'CSKS4': 'X-SAR'}
    instShort = {'S1A': 'C-SAR', 'S1B': 'C-SAR', 'TSX': 'X-SAR',
                 'TDX': 'X-SAR', 'LS8': 'OLI', 'RADARSAT-1': 'C-SAR',
                 'ERS-1/2': 'C-SAR', 'ALOS1-PALSAR': 'PALSAR'}
    instShort = {**CSKinst, **instShort}
    cskIns = {'CSKS1': 'COSMO-SKYMED 1', 'CSKS2': 'X-SAR',
              'CSKS3': 'COSMO-SKYMED 3', 'CSKS4': 'COSMO-SKYMED 4'}
    instrument = {'S1A': 'SENTINEL-1A', 'S1B': 'SENTINEL-1B', 'TSX': 'TSX',
                  'TDX': 'TDX', 'CSK': 'COSMO-SKYMED', 'LS8': 'LANDSAT-8',
                  'RADARSAT-1': 'RADARSAT-1', 'ERS-1/2': 'ERS-1-2',
                  'ALOS1-PALSAR': 'ALOS1-PALSAR'}
    instrument = {**instrument, **cskIns}
    for platform in validPlatforms.keys():
        if validPlatforms[platform]:
            myPlatform = {'AssociatedPlatformShortname': instrument[platform],
                          'AssociatedInstrumentShortname': instShort[platform],
                          'AssociatedSensorShortname': sensorShort[platform]}
            preMetData.append(s.premet(myPlatform, container='Platforms'))
    #
    dataSetID = {'DataSetId': 'MEaSUREs Greenland '
                 f'{prodType(firstDate,lastDate)} Ice Sheet Velocity '
                 f'Mosaics from SAR and Landsat'}
    # force this to be first
    preMetData = [s.premet(dataSetID, container='Collection')] + preMetData
    return preMetData


def makePremets(preMets, outDir, sTypes):
    ''' write the premet file '''
    u.pushd(f'{outDir}/release')
    suffixes = ['vx']  # [],'vy','vv','ex','ey','dT','browse']
    # generate file name using velocity as template
    preMetFiles = []
    regionID = getRegionID()
    for suffix in suffixes:
        velX = u.dols(f'ls -d {regionID}_vel*{suffix}*.tif')
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


def makeSpatial(outDir, spatialFile):
    ''' write spatial file '''
    u.pushd(f'{outDir}/release')
    # assume makePremet error check
    fp = open(spatialFile, 'w')
    _, corners = projWin(iceSheetRegion)
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


def processQuadrants(quadNames, quadrants, outDir, dateTakes, firstDate,
                     lastDate, baseFlags, lsCulledFile, dem, noCull,
                     noReprocess=False):
    ''' loop through and process each quadrant'''
    threads, outputRegions = [], {}
    for quadName in quadNames:
        print(quadName)
        #
        # get lmits and geodat for submosaic
        limits, geodat = findLimitsQ(quadrants[quadName])
        #
        nl = 0
        outPieces = []
        #
        # Setup each submosaic
        for limit in limits:
            print(limit)
            stdout = open(f'{outDir}/io/{quadName}.stdout.{nl:03d}', 'w')
            stderr = open(f'{outDir}/io/{quadName}.stderr.{nl:03d}', 'w')
            # files names for submoaic
            inFileSub = f'inputFile-{quadName}.{nl:03d}'
            outFileSub = f'mosaic-{quadName}.{nl:03d}'
            # save names of ouput mosaics
            outPieces.append(outFileSub)
            # create the input file used to make the mosaic
            if not noReprocess:
                writeInputFileQ(outDir+'/pieces/'+inFileSub, limit, dateTakes,
                                firstDate, lastDate, noCull)
            # Add submosaic to list of threads
            threads.append(threading.Thread(target=runSubMosaic,
                                            args=[outDir, inFileSub,
                                                  outFileSub, baseFlags,
                                                  firstDate, lastDate,
                                                  lsCulledFile, dem, stdout,
                                                  stderr]))
            nl += 1
        # add the list of pieces and corresp geodat to dict index by quad name
        stdout = open(f'{outDir}/io/{quadName}.stdout', 'w')
        stderr = open(f'{outDir}/io/{quadName}.stderr', 'w')
        outputRegions.update({quadName: [outPieces, geodat, stdout, stderr]})
    return threads, outputRegions

# ----------------------------------------------------------------------------
#  main
# ----------------------------------------------------------------------------


def main():
    '''Set up mosaic files, makes directories if needed, produces all files,
    and can run with runAll '''
    #
    # get command line args
    maxThreads = 16
    interpSize = 50
    baseFlags, inputFile, lsFile, firstDate, lastDate, noTSX,\
        outputMaskTemplate, shelfMask, noReprocess, cloudOptimize, noCull \
        = processQArg(sys.argv)
    #
    if iceSheetRegion == 'amundsen':
        dem = defaultAntarcticDEM
    else:
        dem = s.defaultRegionDefs(iceSheetRegion).dem()
    x, y = projWin(iceSheetRegion)
    #
    if len(shelfMask) > 1:
        baseFlags += f' -shelfMask {shelfMask}'
    #
    # grab quadrant info (hardwired for now)
    quadNames, quadrants = setupQuadrants()
    print(quadrants, quadNames)
    # Make output directory
    outDir = f'Vel-{firstDate.strftime("%Y-%m-%d")}.' \
        f'{lastDate.strftime("%Y-%m-%d")}'
    #
    # make subdirectories, set up output mask (generic or specific)
    outputMaskShape = mkQDirs(outDir, outputMaskTemplate)
    #
    # process landsat list (filter by date)
    lsCulledFile = None
    if lsFile is not None:
        lsCulledFile = processLsList(lsFile, firstDate, lastDate, outDir,
                                     noReprocess=noReprocess)
        #
        # parse source list of data takes
    dateTakes = []
    if not noReprocess:
        dateTakes = processInputFileQ(inputFile, noTSX)
    #
    # Loop over each piece and do the relevant processing (threads)
    # and region info (outputRegions)
    threads, outputRegions = processQuadrants(quadNames, quadrants, outDir,
                                              dateTakes, firstDate, lastDate,
                                              baseFlags, lsCulledFile, dem,
                                              noCull,
                                              noReprocess=noReprocess)
    #
    # Run threads if noReprocess=False
    if not noReprocess:
        u.runMyThreads(threads, maxThreads, 'quadrant vel')
    #
    # Combine submosaics
    if True:  # Debug - set false to turn off these steps
        assembleMosiacs(quadNames, outputRegions, outDir, baseFlags)
        # Mask data
        maskMosaics(quadNames, outputRegions, outDir, outputMaskShape,
                    baseFlags)
        #
        # Interpolate mosaics
        interpMosiacs(quadNames, outputRegions, outDir, interpSize, baseFlags)
        #
        # Interpolate tiffs
        writeQTiffs(quadNames, outDir, baseFlags)
        #
        # makevrts
        makeQVRTSs(quadNames, outDir, baseFlags)
        #
        # makeshapes (meta files with frame locations)
        makeShapeOutputs(outDir, firstDate, lastDate)
        #
        # make preview map
    makeLogMap(outDir, quadNames)
    # make final output
    shapeFileNames = makeFinalOutput(outDir, quadNames, cloudOptimize,
                                     firstDate, lastDate)
    #
    sTypes = SatTypes(outDir, shapeFileNames)
    print(sTypes)
    #
    preMetData = []
    preMetData = populatePremet(preMetData, firstDate, lastDate,
                                validPlatforms=sTypes)
    preMetFiles = makePremets(preMetData, outDir, sTypes)
    #
    for shapeNameKey in shapeFileNames:
        print(shapeNameKey, shapeFileNames)
        if shapeFileNames[shapeNameKey] is not None:
            makeSpatial(outDir,
                        shapeFileNames[shapeNameKey].replace('shp', 'spo'))
            makeShapePremet(outDir, shapeFileNames[shapeNameKey], shapeNameKey,
                            sTypes, firstDate, lastDate)
    #
    for preMetFile in preMetFiles:
        if preMetFile is not None:
            print(preMetFile)
            makeSpatial(outDir, preMetFile.replace('premet', 'spo'))


main()
