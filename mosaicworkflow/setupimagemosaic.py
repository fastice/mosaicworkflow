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
from osgeo import gdal

currentVersion = 4
subVersion = 0

suffixDict = {'uncalibrated': ['image'], 'calibrated': ['sigma0', 'gamma0']}
noDataDict = {'uncalibrated': 0, 'calibrated': -30.0}


def setupImageMosaicArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\033[1mSetup and run an image mosaics product \033[0m',
        epilog='Notes:  ', allow_abbrev='False')
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
    parser.add_argument('--check', action='store_true', default=False,
                        help='Setup command, but do not run')
    parser.add_argument('--sensor', type=str, default='S1',
                        help='Sensor (S1)')
    parser.add_argument('--prefix', type=str, default='GL_S1bks_mosaic',
                        help='Prefix for product name [GL_S1bks_mosaic]')
    parser.add_argument('--template', type=str, default='mosaic.template.yaml',
                        help='template that defines mosaic')
    args = parser.parse_args()
    #
    dem, smoothL, DBPaths, tracks, prefix = mosf.mosaicTempArgs(args.template)
    #
    firstDate = datetime.strptime(args.firstdate, "%Y-%m-%d")
    if args.lastdate is None:
        lastDate = firstDate + timedelta(days=args.nDays - 1)
    else:
        lastDate = datetime.strptime(args.lastdate, "%Y-%m-%d")
    #
    calFlag = ['uncalibrated', 'calibrated'][args.calibrate]
    prodPrefix = args.prefix
    if prefix is not None:
        prodPrefix = prefix
    #
    # Some args are hardwired for now but can be updated later
    if tracks is None:
        tracks = ['26', '74', '90', '112', '141', '170', '83']
    # Make sure no reprocessing, including mosaic3d
    if args.noReprocess:
        args.noMosaic = True
    myArgs = {'prefix': prodPrefix,
              'firstDate': firstDate, 'lastDate': lastDate,
              'noReprocess': args.noReprocess, 'check': args.check,
              'calFlag': calFlag, 'sensor': args.sensor,
              'dateDBs': DBPaths, 'geodatFile': 'geodat10x2.in',
              'tracks': tracks, 'noMosaic': args.noMosaic,
              'templateFile': args.template,
              'dem': dem, 'smoothL': smoothL, 'cleanUp': not args.noClean,
              'runMosaics': not args.setupOnly}
    return myArgs


def runGeoMosaic(inputFile, date1, date2, dem, smoothL, calFlag, srsInfo):
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
    Returns
    -------
    None.
    '''
    flags = {'uncalibrated': f'-removePad 5 {smoothL} ',
             'calibrated': f'-removePad 5 {smoothL} -S1Cal '}[calFlag]
    suffixes = suffixDict[calFlag]
    outDir = os.path.dirname(inputFile)
    outFile = os.path.basename(inputFile).replace('inputFile', 'mosaic')
    stdoutFile = f'{outDir}/log/{outFile.replace("mosaic","stdout")}'
    stderrFile = f'{outDir}/log/{outFile.replace("mosaic","stderr")}'
    command = f'cd {outDir}; ' \
        f'geomosaic -center -descending {flags} -xyDEM '\
        f'-date1 {date1.strftime("%m-%d-%Y")} ' \
        f'-date2 {date2.strftime("%m-%d-%Y")} ' \
        f'{os.path.basename(inputFile)} {dem} {outFile}'
    with open(stdoutFile, 'w') as stdout:
        with open(stderrFile, 'w') as stderr:
            print(command, file=stdout)
            # executable='/bin/csh',
            call(command, shell=True, stdout=stdout, stderr=stderr)
            time.sleep(1)  # Add a delay to see if this fixes random seg faults
            for suffix in suffixes:
                fpImageToTiff(f'{outDir}/{outFile}.{suffix}',
                              f'{outDir}/{outFile}.{suffix}',
                              srsInfo, calFlag=calFlag)
                # make browse
                time.sleep(1)


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
                                              myArgs['dem'], myArgs['smoothL'],
                                              myArgs['calFlag'], srsInfo]))
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
    Read the sat types (for now S1A, S1B) from the summary shape file.
    Assumes being called from pieces directory.
    Parameters
    ----------
    shapeFilesName : shapefile name (xxx.shp).
    Returns
    -------
    sTypes : {'S1A': True/False, 'S1B': True/False} .

    '''
    u.pushd('../release')
    # any records, record LS8 present
    sTypes = {'S1A': False, 'S1B': False}
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
    sensorShort = {'S1A': 'C-SAR', 'S1B': 'C-SAR'}
    instShort = {'S1A': 'C-SAR', 'S1B': 'C-SAR'}
    instrument = {'S1A': 'Sentinel-1A', 'S1B': 'Sentinel-1B'}
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
    # print(sTypes)
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
                    cleanUp=True, noReprocess=False, calFlag='uncalibrated'):
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
    Returns
    -------
    satTypes : dict indicating which sats used.
    '''
    # prodTypes = {'uncalibrated': ['uncalibrated'],
    #              'calibrated': ['calibrated', 'calibrated']}[calFlag]
    suffixes = suffixDict[calFlag]
    noData = noDataDict[calFlag]
    jpgRes = '500'
    try:  # Use try so cleanUp won't occur if there is a problem
        u.pushd(directory=f'{prodDir}/pieces')
        #
        # make the shape file
        prodName = makeProdName(prefix, date1, date2)
        shapeFileName = prodName.replace('*', 'shape').replace('_RESm', '')
        command = f'makeimageshapefile.py inputFile.0.0 ' \
            f'../release/{shapeFileName}'
        print(command)
        call(command, shell=True)  # , executable='/bin/csh')
        satTypes = imageSatTypes(f'{shapeFileName}.shp')
        makeImageSpatial(corners, '..', f'{shapeFileName}.spo')
        # added Aug 20 2020
        makeShapePremet(f'../release/{shapeFileName}.shp', satTypes, date1,
                        date2)
        # build the vrt
        if not noReprocess:
            for suffix in suffixes:  # Loop on prod types
                # Remove prexisting
                if os.path.exists(f'{prodDir}.{suffix}.vrt'):
                    os.remove(f'{prodDir}.{suffix}.vrt')
                # generate vrt from tiffs
                command = f'gdalbuildvrt -vrtnodata {noData} '\
                    f'{prodDir}.{suffix}.vrt *{suffix}.tif'
                print(command)
                call(command, shell=True)  # , executable='/bin/csh')
                # create the cloud optimized geo
                tifName = \
                    prodName.replace("*", suffix).replace("RES", str(posting))
                jpgName = prodName.replace("*", suffix).replace("RES", jpgRes)
                # create image and stats
                # command = f'rio cogeo create {prodDir}.{suffix}.vrt ' \
                #    f'../release/{tifName}.tif --overview-resampling ' \
                #    f'average ; gdalinfo -stats ../release/{tifName}.tif'
                newTif = f'../release/{tifName}.tif'
                origVrt = f'{prodDir}.{suffix}.vrt'
                command = 'gdal_translate -of COG -co COMPRESS=DEFLATE ' \
                    '-co RESAMPLING=AVERAGE -co OVERVIEWS=IGNORE_EXISTING ' \
                    '-co GEOTIFF_VERSION=1.1 -co BIGTIFF=NO -stats ' \
                    f'{origVrt} {newTif}'
                # print(command)
                # u.myerror('debug')
                call(command, shell=True)  # , executable='/bin/csh')
                #
                # make quicklook
                call('gdal_translate -co "QUALITY=99" -scale -of JPEG -r '
                     f'average -tr {jpgRes} {jpgRes} '
                     f'../release/{tifName}.tif ../release/{jpgName}.jpg',
                     shell=True)  # , executable='/bin/csh')
    except Exception:
        u.myerror('mergeImageTiles: Problem creating final mosaic')
    #
    if cleanUp:
        # remove all image files - leave inputs and vrts
        toClean = u.dols('ls mosaic.*.*')
        for mosaicFile in toClean:
            os.remove(mosaicFile)
    u.popd()
    #
    return satTypes


def makeImagePremets(outDir, prefix, firstDate, lastDate, sTypes, posting,
                     calFlag='uncalibrated'):
    ''' write the premet file '''
    u.pushd(f'{outDir}/release')
    preMetFiles = []
    suffixes = suffixDict[calFlag]
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
    print(type(myArgs['firstDate']),  type(sarDB.minDate))
    print(myArgs['firstDate'],  sarDB.minDate)
    print(myArgs['lastDate'], sarDB.maxDate)
    print(myArgs['firstDate'] >= sarDB.minDate)
    print(myArgs['lastDate'] <= sarDB.maxDate)
    return myArgs['firstDate'] >= sarDB.minDate and \
        myArgs['lastDate'] <= sarDB.maxDate


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
    print(myArgs)
    #
    checkGdalVersion()
    #
    myDB = mosf.sarDB()
    myDB.readDB(myArgs['dateDBs'])
    # myDB.printByDate()
    #
    print(checkDateRange(myArgs, myDB))
    if not checkDateRange(myArgs, myDB):
        u.myerror('Requested range falls out side of DB date range:'
                  f'{myDB.minDate} {myDB.maxDate}\nUpdate DB or adjust dates')
    # create single dict by date of products across multiple dbs with paths
    dateRangeAll = myDB.extractByDateRange(myArgs['firstDate'],
                                           myArgs['lastDate'])

    dateRangeTracks = dateRangeAll.selectTrackDateRangeDB(myArgs['tracks'])
    #
    prodDir = setupMosaicDir(myArgs['firstDate'], myArgs['lastDate'],
                             f'{myArgs["prefix"]}_{myArgs["calFlag"]}')
    #
    inputFiles, srsInfo, corners, posting = \
        mosf.makeSectionedPowerInputFiles(f'{prodDir}/pieces',
                                          myArgs['templateFile'],
                                          dateRangeTracks,
                                          myArgs['geodatFile'])
    # Determine max threads
    maxThreads = min(len(inputFiles), int(multiprocessing.cpu_count()/2))
    # setup tiles
    threads = setupGeoMosaic(inputFiles, myArgs, srsInfo)
    # produce tiles
    if myArgs['runMosaics'] and not myArgs['noMosaic']:
        u.runMyThreads(threads, maxThreads, 'geomosaic')
    # merge tiles and produce mosaic in release dir
    satTypes = mergeImageTiles(prodDir, myArgs['prefix'], myArgs['firstDate'],
                               myArgs['lastDate'], posting, corners,
                               calFlag=myArgs['calFlag'],
                               noReprocess=myArgs['noReprocess'],
                               cleanUp=myArgs['cleanUp'])
    # Create and save premet file(s)
    preMetFiles = makeImagePremets(prodDir, myArgs['prefix'],
                                   myArgs['firstDate'],
                                   myArgs['lastDate'], satTypes, posting,
                                   calFlag=myArgs['calFlag'])

    # write spatial file for each premet
    for preMetFile in preMetFiles:
        makeImageSpatial(corners, prodDir, preMetFile.replace('premet', 'spo'))


if __name__ == "__main__":
    main()
