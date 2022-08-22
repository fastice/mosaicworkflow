#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 13:08:29 2019
@author: ian
"""
import argparse
import utilities as u
from datetime import datetime, timedelta
import os
from subprocess import call
import sarfunc as s
import threading
import shutil
import random
from operator import itemgetter
currentVersion = '04'
subVersion = '0'


def setupCleanTSXArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(description='\n\n\033[1mRun steps to '
                                     'produce a TSX or CSK release - '
                                     'run in release dir \033[0m\n\n')
    parser.add_argument('--track', type=int, default=-1, help='specify integer'
                        ' track number [*]')
    parser.add_argument('--firstdate', type=str, default='1900-01-31',
                        help='Only use central dates>=first date [1900-01-31]')
    parser.add_argument('--lastdate', type=str, default='2100-01-01',
                        help='Only use central dates<= last date [2100-01-31]')
    parser.add_argument('--noPrompt', action='store_true', default=False,
                        help='Run with no prompt to start threads')
    parser.add_argument('--delivered', type=str, default=None,
                        help='Listing of files aleady delivered [None]')
    parser.add_argument('--keepFast', action='store_true', default=False,
                        help='Do not mask fast tracked areas to keep melange')
    parser.add_argument('--sensor', type=str, default='TSX',
                        choices=['TSX', 'CSK'], help='Sensor name')
    parser.add_argument('--region', type=str, default='greenland',
                        choices=['greenland', 'taku'], help='Sensor name')
    parser.add_argument('releaseDir', type=str, nargs=1,
                        help='Name of release dir for this release')
    #
    args = parser.parse_args()
    #
    track = args.track
    if track <= 0:
        track = getTrack()
    #
    firstdate = checkDate(args.firstdate)
    lastdate = checkDate(args.lastdate)
    flags = {'default': None}
    #
    usePrompt = not args.noPrompt
    #
    pwd = os.getcwd()
    if 'release' not in pwd.split('/'):
        u.myerror(f'Error for cwd {pwd}: Much run in release dir')
    if args.sensor not in pwd:
        u.myerror(f'Error for sensor {args.sensor} not path name: {pwd}')
    #
    delivered = readDelivered(args.delivered)
    return args.releaseDir[0], track, firstdate, lastdate, flags, usePrompt, \
        args.keepFast, delivered, args.sensor, args.region


def getMaskDict(glacierID):
    ''' this reads the special masks for Jak, Helheim, and Kang'''
    #
    if glacierID not in ['Jak', 'Helheim', 'Kang']:
        return None
    masks = u.dols(
        f'ls /Users/ian/greenlandmask/{glacierID}Output/mask*.????.???')
    currentYear, currentDay = 2009, 1
    maskDict = {}
    for mask in masks:
        pieces = mask.split('.')
        maskYear, maskDay = int(pieces[-2]), int(pieces[-1])
        # loop until date matches mask
        while True:
            maskDict[f'{currentYear}-{currentDay}'] = mask
            currentDay += 1
            if currentDay > datetime(currentYear, 12, 31).timetuple().tm_yday:
                currentDay = 1
                currentYear += 1
            if currentDay == maskDay and currentYear == maskYear:
                break
    return maskDict


def readDelivered(listingFile):
    '''
    Read a listing file that contains a listing of products already delivered

    Parameters
    ----------
    def listingFile : str
        File with list of prior releases "boxName productName".

    Returns
    -------
    delivered. Dictionary by box with products delivered
    '''
    if listingFile is None:
        return {}
    if not os.path.exists(listingFile):
        u.myerror('readListing: listing file {listingFile} does not exist')
    delivered = {}
    with open(listingFile) as fpListing:
        for line in fpListing:
            if 'Vel' in line:
                box, product = [x.strip() for x in line.split()]
                if box not in delivered:
                    delivered[box] = []
                delivered[box].append(product)
    return delivered


def checkDate(dateString):
    try:
        myDate = datetime.strptime(dateString, "%Y-%m-%d")
    except Exception:
        u.myerror("Invalid date string {dateString} - use YYYY-MM-DD")
    return myDate


def getTrack():
    tracks = u.dols('ls -d ../track-*')
    if len(tracks) == 0:
        u.myerror('No track dir present')
    if len(tracks) > 1:
        return '*'
    try:
        track = int(tracks[0].split('-')[1])
    except Exception:
        u.myerror(f'Cannot parse track from {tracks[0]}')
    return track


def mkMyDir(releaseDir, dirType, verbose=False):
    if not os.path.exists(releaseDir):
        try:
            if verbose:
                print(f'Creating {dirType} dir {releaseDir}')
            os.mkdir(releaseDir)
        except Exception:
            u.myerror(f'Error creating release dir: {releaseDir}')


def prepReleaseDir(releaseDir, boxName):
    ''' Make directory for release with subdirectories for intermediate'''
    if len(releaseDir) < 4:
        u.myerror(f'Invalid Release dir {releaseDir}: must be a least 4 char')
    #
    mkMyDir(releaseDir, 'release', verbose=True)
    mkMyDir(f'{releaseDir}/intermediate', 'release', verbose=True)
    mkMyDir(f'{releaseDir}/{boxName}', 'release', verbose=True)


def checkOffs(offDirsUse):
    parts = ['azimuth.offsets', 'range.offsets', 'geodat3x3.in',
             'motion/baseline.3x3', 'motion/az.est', 'motion/rBaseline',
             'velocity/inputFile']
    for offDir in offDirsUse:
        for part in parts:
            for myDir in offDir["dir"]:
                if not os.path.exists(f'{myDir}/{part}'):
                    u.myerror(f'{offDir["dir"]} missing {part} for date '
                              f'{offDir["date"]}')


def computeNominalTimes(offDirsUse):
    ''' compute times from meta data - average if more than one'''
    # loop through products
    for offDir in offDirsUse:
        seconds, count = 0., 0
        # loop through meta (typically 1 or 2)
        for meta in offDir['meta']:
            h, m, s = meta['Nominal Time for Pair (HH:MM:SS)'].split(':')
            seconds += int(h) * 3600. + int(m) * 60. + int(s)
            count += 1
        tmpDate = datetime(1, 1, 1) + timedelta(seconds=seconds/count)
        offDir['nomTime'] = tmpDate.strftime('%H-%M-%S')


def populatePremet(offDirsUse, sensor):
    ''' populate premet data fields '''
    platform = {'AssociatedPlatformShortName': sensor,
                'AssociatedInstrumentShortName': 'X-SAR',
                'AssociatedSensorShortName': 'X-SAR'}
    for offDir in offDirsUse:
        offDir['premet'].append(
            s.premet({'VersionID_local': f'0{currentVersion}'}))
        offDir['premet'].append(
            s.premet(
                {'Begin_date': f'{offDir["date1"].strftime("%Y-%m-%d")}'}))
        offDir['premet'].append(
            s.premet({'End_date': f'{offDir["date2"].strftime("%Y-%m-%d")}'}))
        offDir['premet'].append(s.premet({'Begin_time': '00:00:01.000'}))
        offDir['premet'].append(s.premet({'End_time': '23:59:59.000'}))
        # offDir['premet'].append(s.premet({f'DataSetId': 'MEaSUREs Greenland
        # Ice Velocity: Selected Site Velocity Maps from InSAR
        # V00{currentVrsion}' }, container='Collection'))
        offDir['premet'].append(
            s.premet(platform, container='AssociatedPlatformInstrumentSensor'))


def getOffsetsToUse(track, firstDate, lastDate, deliveredList, sensor):
    '''
    Go through velocity tracks and build a list of products, merging
    products form same day
    '''
    print(f'get offsets for {track}')
    velDirs = u.dols(f'ls -d ../track-{track}/*_*/velocity')
    # Gather pairs and dates
    myPairs = []
    for velDir in velDirs:
        myMeta = s.parseVelMeta(f'{velDir}/mosaicOffsets.meta')
        myDate1 = datetime.strptime(myMeta['First Image Date (MM:DD:YYYY)'],
                                    '%b:%d:%Y')
        myDate2 = datetime.strptime(myMeta['Second Image Date (MM:DD:YYYY)'],
                                    '%b:%d:%Y')
        Exclude = os.path.exists(velDir.replace('velocity', 'Exclude'))
        productDir = f'Vel.{myDate1.strftime("%Y-%m-%d")}.' \
            f'{myDate2.strftime("%Y-%m-%d")}'
        # Skip those marked Exclude or already delivered
        if not Exclude and productDir not in deliveredList:
            myPairs.append({'dir': [os.path.dirname(velDir)],
                            'productDir': productDir,
                            'date1': myDate1, 'date2': myDate2,
                            'meta': [myMeta], 'premet': []})
    # sort results
    myPairs = sorted(myPairs, key=itemgetter('date1'))
    #
    offDirsUse = []
    for myPair in myPairs:
        centralDate = myPair['date1'] + (myPair['date2']-myPair['date1'])*0.5
        if centralDate >= firstDate and centralDate <= lastDate:
            # check this is not a multipart case
            # if this is same as lastpair, append to prior-relies on sort above
            found = False
            if len(offDirsUse) > 0:
                if offDirsUse[-1]['date1'] == myPair['date1'] and \
                        offDirsUse[-1]['date2'] == myPair['date2']:
                    offDirsUse[-1]['dir'].append(myPair['dir'][0])
                    offDirsUse[-1]['meta'].append(myPair['meta'][0])
                    found = True
            # default case, single pair for the products
            if not found:
                offDirsUse.append(myPair)
    checkOffs(offDirsUse)
    computeNominalTimes(offDirsUse)
    populatePremet(offDirsUse, sensor)
    return offDirsUse


def getBoxInfo(region):
    if region == 'greenland':
        boxName = u.dols('ls -d *coast-*')[0]
    elif region == 'taku':
        boxName = u.dols('ls -d Alaska-*')[0]
    else:
        u.myerror(f'{region} not valid region')
    resFile = f'{boxName}/resolution'
    if os.path.exists(resFile):
        fp = open(resFile)
        for line in fp:
            if 'resolution' in line:
                resString = line.strip().split('=')[1].replace('"', '')
                fp.close()
                break
    else:
        u.myerror(f'Missing resolution file for {boxName}')

    resP = [float(x) for x in resString.split()]
    #
    myParts = os.getcwd().split('/')
    if myParts[-1] != 'release':
        u.myerror('must run in release directory')
    glacierID = myParts[-2]
    #
    if len(resP) != 6:
        u.myerror(f'Invalid resolution string {resString} for {boxName}')
    return boxName, resString, glacierID


def mkInputFile(offDir, outDir, resolution, myInputs):
    fp = open(f'{outDir}/inputFile', 'w')
    print(f'; Generated using inputs {outDir}/velocity/inputFile', file=fp)
    print(resolution, file=fp)
    print(';', file=fp)
    print(len(myInputs), file=fp)
    print(';', file=fp)
    for myInput in myInputs:
        print(myInput, file=fp)
    fp.close()


def getInputData(offDir):
    try:
        myLines = []
        for myDir in offDir["dir"]:
            fp = open(f'{myDir}/velocity/inputFile')
            for line in fp:
                if 'rBaseline' in line:
                    myLines.append(line)
                    break
            fp.close()
        return myLines
    except Exception:
        u.myerror(
            f'Error reading input file for {offDir["dir"]}/velocity/inputFile')
    return None


def makeRunOff(outDir, shelfMask, dem):
    try:
        fp = open(f'{outDir}/runOff', 'w')
        print(f'#\nset DEM={dem}\n#', file=fp)
        print(f'mosaic3d -fl 20 -offsets -rOffsets -center -no3d -noVh '
              f'-shelfMask {shelfMask} inputFile $DEM mosaicOffsets', file=fp)
    except Exception:
        u.myerror(f'Problem making runOff in {outDir}')


def getDem(offDir, boxName, demDict, region):
    ''' get dem from dict or return defaull'''
    # dem = '/Volumes/insar3/ian/gimp/gimpdem_90m.270'
    # updated 5/1/2020 to use gimp1 < 2015, ow gimp2
    if region == 'greenland':
        if offDir["date1"].year < 2015:
            dem = '/Volumes/insar7/ian/gimp/gimp1/270m/dem.gimp1.270m'
        else:
            # dem = '/Volumes/insar7/ian/gimp/gimp2/270m/dem.gimp2.270m'
            dem = s.defaultRegionDefs('greenland').dem()
    else:
        dem = s.defaultRegionDefs(region).dem()
    # no custom DEM return
    if demDict is None:
        return dem
    # if entry for year, use it.
    demKey = f'dem{offDir["date1"].year}'
    try:
        if os.path.exists(demDict[demKey]):
            dem = demDict[demKey]
    except Exception:
        u.mywarning(f'dem dict found but no entry for date {offDir["date1"]}')
    return dem


def getDemDict(boxName):
    ''' load the dem dictionary'''
    if not os.path.exists(f'{boxName}/dems'):
        return None
    # read DEM dict if it exists - the dem is specified in the box
    # (e.g. WcoastX.Y/dems)
    fp = open(f'{boxName}/dems')
    demDict = {}
    for line in fp:
        if len(line) > 10 and '=' in line:
            demDict[line.split('=')[0].strip()] = line.split('=')[1].strip()
    fp.close()
    return demDict


def makeMosaic3D(outDir):
    ''' run the  mosaicker to produce the product'''
    command = f'cd {outDir} ; csh runOff'
    n = random.randint(1, 10000)
    fout = open(f'{outDir}/stdout.{n}', 'w')
    ferr = open(f'{outDir}/stderr.{n}', 'w')
    call(command, shell=True, executable='/bin/csh', stdout=fout, stderr=ferr)
    fout.close()
    ferr.close()
    os.remove(f'{outDir}/stdout.{n}')
    os.remove(f'{outDir}/stderr.{n}')


def makeJPG(outDir, productDir, velRoot, sensor):
    ''' Make quicklook '''
    spaceAgency = {'TSX': 'DLR provided TerraSAR-X/TanDEM-X data',
                   'CSK': 'ASI provided COSMO-SkyMed data'}
    n = random.randint(1, 10000)
    with open(f'{outDir}/stdout.{n}', 'w') as fout, \
            open(f'{outDir}/stderr.{n}', 'w') as ferr:
        name = velRoot.replace('_TBD', '').replace('.tif', '')
        command = f'cd {outDir} ; showvel.py -fig={name} -title={name} ' \
            f'-subtitle="''NASA MEaSUREs GIMP product produced using ' \
            f'{spaceAgency[sensor]}''"'
        call(command, shell=True, executable='/bin/csh', stdout=fout,
             stderr=ferr)
        os.rename(f'{outDir}/{name}.jpg', f'{productDir}/{name}.jpg')
    #
    os.remove(f'{outDir}/stdout.{n}')
    os.remove(f'{outDir}/stderr.{n}')


def writeTiffProduct(productDir, outDir, velRoot, productType, epsg, wktFile):
    ''' write velocity or error product '''
    noData = {'velocity': -2.0e9, 'error': -2.0e9}
    suffixDictIn = {'velocity': ['vx', 'vy', 'v'], 'error': ['ex', 'ey']}
    suffixDictOut = {'velocity': ['vx', 'vy', 'vv'], 'error': ['ex', 'ey']}
    myVel = u.geoimage(geoType=productType, verbose=False)
    myVel.readData(f'{outDir}/mosaicOffsets', epsg=epsg)
    # write with default names
    myVel.writeMyTiff(f'{productDir}/mosaicOffsets', epsg=epsg,
                      wktFile=wktFile,
                      noDataDefault=noData[productType], predictor=1,
                      noV=False, overviews=[2, 4])
    #
    suffixesIn = suffixDictIn[productType]
    suffixesOut = suffixDictOut[productType]
    for suffixIn, suffixOut in zip(suffixesIn, suffixesOut):
        myTiff = f'{productDir}/mosaicOffsets.{suffixIn}.tif'
        if not os.path.exists(myTiff):
            u.myerror(f'Error processing {myTiff}, which does not exist')
        myNewTiff = f'{productDir}/{velRoot.replace("TBD", suffixOut)}'
        #
        os.rename(myTiff, myNewTiff)
    return myVel.computePixEdgeCornersLL()


def makeFinalProductDir(releaseDir, offDir, boxName, sensor):
    ''' make directory to write the final product too '''
    productDir = f'{releaseDir}/{boxName}/{offDir["productDir"]}'
    mkMyDir(productDir, 'productDir')
    boxShort = boxName.replace('coast-', '')
    #
    velRoot = f'{sensor}_{boxShort}_{offDir["date1"].strftime("%d%b%y")}_' \
        f'{offDir["date2"].strftime("%d%b%y")}_' \
        f'{offDir["nomTime"]}_TBD_v{currentVersion}.' \
        f'{subVersion}.tif'
    return productDir, velRoot


def makePremet(preMets, productDir, velRoot):
    ''' write the premet file '''
    premetFile = velRoot.replace(".tif", ".premet").replace("_TBD", "")
    fp = open(f'{productDir}/{premetFile}', 'w')
    # prepend filename
    DataSet = s.premet({'DataSetId': 'MEaSUREs Greenland Ice Velocity:'
                        f'Selected Site Velocity Maps from InSAR'
                        f' V{currentVersion}.{subVersion}'},
                       container='Collection')
    GranID = s.premet({'ProducerGranuleID':
                       velRoot.replace("_TBD", "").replace(".tif", "")},
                      container='DataGranule')
    preMets = [DataSet, GranID]+preMets
    for preMet in preMets:
        preMet.printParams(fp=fp)
    fp.close()


def makeSpatial(corners, productDir, velRoot):
    ''' write spatial file '''
    fp = open(
        f'{productDir}/{velRoot.replace(".tif", ".spo").replace("_TBD", "")}',
        'w')
    print(f'{corners["ll"]["lon"]:10.5f} {corners["ll"]["lat"]:10.5f}',
          file=fp)
    print(f'{corners["ul"]["lon"]:10.5f} {corners["ul"]["lat"]:10.5f}',
          file=fp)
    print(f'{corners["ur"]["lon"]:10.5f} {corners["ur"]["lat"]:10.5f}',
          file=fp)
    print(f'{corners["lr"]["lon"]:10.5f} {corners["lr"]["lat"]:10.5f}',
          file=fp)

    fp.close()


def runReleaseVel(releaseDir, outDir, offDir, boxName, epsg, wktFile, sensor):
    ''' this routine is threaded to make the products'''
    # make the products
    makeMosaic3D(outDir)
    #
    productDir, velRoot = makeFinalProductDir(releaseDir, offDir, boxName,
                                              sensor)
    # write the tiff files
    cornersV = writeTiffProduct(productDir, outDir, velRoot, 'velocity', epsg,
                                wktFile)
    _ = writeTiffProduct(productDir, outDir, velRoot, 'error', epsg, wktFile)
    #
    shutil.copy(f'{outDir}/mosaicOffsets.meta',
                f'{productDir}/'
                f'{velRoot.replace(".tif", ".meta").replace("_TBD", "")}')
    #
    makeJPG(outDir, productDir, velRoot, sensor)
    #
    makePremet(offDir['premet'], productDir, velRoot)
    #
    makeSpatial(cornersV, productDir, velRoot)


def shelfMaskName(myDate, maskDict, glacierID, keepFast, region):
    ''' get the mask name '''
    #
    if region != 'greenland':
        return s.defaultRegionDefs(region).mask()
    if keepFast:
        return '/Users/ian/greenlandmask/MelangeMask/melangemask'
    doy = myDate.timetuple().tm_yday
    myYear = myDate.year
    # handle special cases
    if glacierID in ['Jak', 'Helheim', 'Kang'] and maskDict is not None:
        myKey = f'{myYear}-{doy}'
        if myKey in maskDict.keys():
            return maskDict[myKey]
    # use default seasonal masks
    if doy <= 30:
        myYear = myYear - 1
        season = 'ew'
    elif doy > 30 and doy <= 100:
        season = 'lw'
    elif doy > 100 and doy <= 160:
        season = 'sp'
    elif doy > 160 and doy <= 208:
        season = 'ms'
    elif doy > 208 and doy <= 293:
        season = 'ls'
    elif doy > 293:
        season = 'ew'
    # nominal case return
    shelfMask = f'/Users/ian/greenlandmask/output{myYear}/' \
        f'greenlandmask{season}{myYear}'
    if os.path.exists(shelfMask):
        return shelfMask
    if myYear == 2008:
        shelfMask = f'/Users/ian/greenlandmask/output{myYear}/'\
            f'greenlandmaskew{myYear}'
        return shelfMask
    # try one back with warning
    backDict = {'ew': 'ls', 'ls': 'ms', 'ms': 'sp', 'sp': 'lw', 'lw': 'ew'}
    shelfMaskBackup = f'/Users/ian/greenlandmask/output{myYear}/'\
        f'greenlandmask{backDict[season]}{myYear}'
    # check it exists
    if os.path.exists(shelfMaskBackup):
        u.mywarning(f'For {myDate.strftime("%Y-%m-%d")} shelfmask, {shelfMask}'
                    f', does not exist. Proceeding with {shelfMaskBackup}')
        return shelfMaskBackup
    u.myerror(f'Shelfmask does not exist: {shelfMask} and back up '
              f'{shelfMaskBackup} does not either')


def setupProductGen(releaseDir, offDirsUse, boxName, resolution, firstDate,
                    lastDate, epsg, wktFile, maskDict, glacierID, keepFast,
                    sensor, region):
    ''' set up mosaicker '''
    threads = []
    demDict = getDemDict(boxName)
    for offDir in offDirsUse:
        outDir = f'{releaseDir}/intermediate/{offDir["productDir"]}'
        offDir['intermediate'] = outDir
        # mk the directory
        mkMyDir(outDir, 'intermediate')
        # grab the input data string from velocity
        myInputs = getInputData(offDir)
        # create the input file
        mkInputFile(offDir, outDir, resolution, myInputs)
        #
        shelfMask = shelfMaskName(offDir["date1"], maskDict, glacierID,
                                  keepFast, region)
        #
        dem = getDem(offDir, boxName, demDict, region)
        # make the runOff
        makeRunOff(outDir, shelfMask, dem)
        #
        threads.append(threading.Thread(target=runReleaseVel,
                                        args=[releaseDir, outDir, offDir,
                                              boxName, epsg, wktFile, sensor]))
    return threads


def main():
    # assume for now this greenland only
    # epsg = 3413
    # get args
    releaseDir, track, firstDate, lastDate, flags, usePrompt, keepFast, \
        delivered, sensor, region = setupCleanTSXArgs()
    print(f'track {track} firstdate {firstDate} lastdate {lastDate}')
    #
    boxName, resolution, glacierID = getBoxInfo(region)
    maskDict = getMaskDict(glacierID)
    prepReleaseDir(releaseDir, boxName)
    print(f'Box {boxName} with resolution {resolution}')
    deliveredList = []
    if boxName in delivered:
        deliveredList = delivered[boxName]
    # get offsets to use based on date and track info
    offDirsUse = getOffsetsToUse(track, firstDate, lastDate, deliveredList,
                                 sensor)
    # print(deliveredList)
    # setup the product threads
    epsg = s.defaultRegionDefs(region).epsg()
    wktFile = s.defaultRegionDefs(region).wktFile()
    velThreads = setupProductGen(releaseDir, offDirsUse, boxName, resolution,
                                 firstDate, lastDate, epsg, wktFile, maskDict,
                                 glacierID, keepFast, sensor, region)
    # exit()
    # run threads
    u.runMyThreads(velThreads, 12, 'Release Vel', delay=0, prompt=usePrompt)


main()
