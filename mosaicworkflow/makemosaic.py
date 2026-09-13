#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 08:45:12 2020

@author: ian
"""

import argparse
import sys
import utilities as u
import os
from datetime import datetime, timedelta, date
from subprocess import call
import mosaicfunc as mosf


# Keys understood in mosaic.template.yaml (and, for sepAscDesc/projectFile,
# optionally in project.yaml -- see setupquarters.resolveSepAscDesc). This
# list is informational only, kept by hand since the schema isn't read from
# a single place in the code.
TEMPLATE_KEYS = [
    ('Required', [
        ('dem', None, 'DEM file (overridable via --dem)'),
        ('xll', None, 'lower-left x coordinate (km)'),
        ('yll', None, 'lower-left y coordinate (km)'),
        ('sx', None, 'mosaic width (pixels)'),
        ('sy', None, 'mosaic height (pixels)'),
        ('dx', None, 'pixel spacing in x (km)'),
        ('dy', None, 'pixel spacing in y (km)'),
        ('nx', None, 'number of pieces to section the mosaic into (x)'),
        ('ny', None, 'number of pieces to section the mosaic into (y)'),
        ('regionID', None, 'region identifier'),
        ('baseFlags', None, 'flags passed to mosaic3d (overridable via --baseFlags)'),
    ]),
    ('Optional, overridable via matching command line flag', [
        ('lsFile', 'None', 'Landsat input file for mosaic period'),
        ('outputMask', 'None', 'shapefile for final output mask'),
        ('mosaicMask', 'None', 'mask passed to mosaic3d (shelfMask)'),
        ('inputFile', 'None', 'input file with SAR inputs'),
    ]),
    ('Optional, projection/region', [
        ('epsg', 'None', 'EPSG code; may be omitted if wktFile is given'),
        ('wktFile', 'None', 'WKT projection file (takes precedence over epsg)'),
        ('regionFile', None, 'path to a separate region-definition yaml; its keys '
                              '(epsg/wktFile/dem/velMap/sigmaShape) are merged into '
                              'the template'),
        ('velMap', None, 'velocity map for region (warned if missing and no regionFile)'),
        ('sigmaShape', None, 'sigma shapefile for region (warned if missing and no regionFile)'),
    ]),
    ('Optional, preview image labeling', [
        ('xLabel', None, 'x pixel offset for credit label on preview image'),
        ('yLabel', None, 'y pixel offset for credit label on preview image '
                          '(both xLabel and yLabel must be set, else labeling is skipped)'),
    ]),
    ('Optional, ascending/descending pairing (see setupquarters.resolveSepAscDesc)', [
        ('sepAscDesc', 'True (mosaic3d default)',
         'False adds -noSepAscDesc to baseFlags; None defers to project.yaml; '
         'errors if it contradicts -noSepAscDesc already present in baseFlags'),
        ('projectFile', 'None (search upward from the template for project.yaml)',
         'explicit path to project.yaml, resolved relative to the template\'s own '
         'directory if not absolute; must exist if given'),
    ]),
    ('Optional, squint correction (see setupquarters.resolveUseSquint)', [
        ('useSquint', 'False (program default)',
         'True appends -useSquint to baseFlags; False leaves it off; '
         'None defers to project.yaml useSquint key; '
         'template takes precedence over project.yaml'),
    ]),
    ('Optional, azimuth ionosphere correction (see setupquarters.resolveBooleanBaseFlag)', [
        ('applyAzimuthIonosphereCorrection', 'False (program default)',
         'True appends -useAzIonosphere to baseFlags, so mosaic3d applies the '
         'azimuth ionosphere correction for frames whose azparams fit chose it; '
         'False leaves it off; None defers to project.yaml'),
    ]),
    ('Optional, crossing-orbit time thresholds (see setupquarters.resolveNumericBaseFlag)', [
        ('timeThresh', 'mosaic3d default (12 days)',
         'Max days between crossing-orbit offset pairs; appends -timeThresh N to baseFlags; '
         'template takes precedence over project.yaml'),
        ('timePhaseThresh', 'mosaic3d default (548 days)',
         'Max days between crossing-orbit phase pairs; appends -timePhaseThresh N to baseFlags; '
         'template takes precedence over project.yaml'),
    ]),
    ('Optional, joint crossing-orbit solvers (mosaic3d default since 2026-08-29)', [
        ('jointMaxSigmaPhase', 'mosaic3d default (35 m/yr)',
         'Drop crossing-PHASE pixels whose sigmaWorst*sqrt(n) exceeds this -- i.e. whose '
         'effective per-measurement sigma is too large.  n-normalised so it does not '
         'penalise thin coverage.  0 disables; template takes precedence over project.yaml'),
        ('jointMaxSigmaRange', 'mosaic3d default (100 m/yr)',
         'Same for crossing-RANGE.  Looser because offset errors are proportionally smaller '
         'on fast ice (30 m/yr on 10 km/yr is 0.3%).  0 disables'),
        ('legacyPairPhase', 'False (joint solver is the default)',
         'True appends -legacyPairPhase, restoring the ORIGINAL pairwise crossing-phase solver'),
        ('legacyPairRange', 'False (joint solver is the default)',
         'True appends -legacyPairRange, restoring the ORIGINAL pairwise crossing-range solver'),
    ]),
]


def printTemplateKeys():
    ''' Print all keys understood in mosaic.template.yaml, with defaults '''
    print('\nKeys understood in mosaic.template.yaml:\n')
    for section, keys in TEMPLATE_KEYS:
        print(f'{section}:')
        for key, default, description in keys:
            defaultStr = '' if default is None else f' [default: {default}]'
            print(f'  {key}{defaultStr}\n      {description}')
        print()


def makemosaicArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\033[1mFront end to run velocity mosaics '
        'using setupquarters.py \033[0m',
        epilog='Notes:\nPart of the mosaicworkflow package.', allow_abbrev='False')
    parser.add_argument('--firstdate', type=str, default='2015-01-01',
                        help='Use central dates >= first date [2015-01-01]')
    defaultEnd = date.today().strftime('%Y-%m-%d')
    parser.add_argument('--lastdate', type=str, default=defaultEnd,
                        help=f'Use central dates <= lastdate [{defaultEnd}]')
    parser.add_argument('--inputFile', type=str, default=None,
                        help='Input file for radar data [from mosaicTemplate]')
    parser.add_argument('--interval', type=str,
                        choices=['s1cycle', 's1-12day', 'monthly', 'quarterly',
                                 'annual', 'multiYear', 'quarterlyJFM'],
                        default='monthly',
                        help='time interval spanned')
    # parser.add_argument('--region', type=str,
    #                     choices=['greenland', 'amundsen', 'taku', 'custom'],
    #                     default='greenland',
    #                     help='region code if not customRegion')
    parser.add_argument('--noReprocess', action='store_true', default=False,
                        help='Just reformat data')
    parser.add_argument('--noLandsat', action='store_true', default=False,
                        help='Do not include landsat')
    parser.add_argument('--landsatPath', type=str,
                        default=None,
                        help='Path to directory with lists of landsat scenes')
    parser.add_argument('--LSFitType', type=str,
                        choices=['custom.mask', 'custom', 'static',
                                 'static.mask'],
                        default='custom.mask',
                        help='Landsat fit type')
    parser.add_argument('--mosaicMask', type=str, default=None,
                        help='Mask applied by mosaic3d')
    parser.add_argument('--outputMask', type=str, default=None,
                        help='Shapefile for final output mask')
    parser.add_argument('--baseFlags', type=str, default=None,
                        help='Override default flags for the mosaicker')
    parser.add_argument('--keepFast', action='store_true', default=False,
                        help='Do not mask fast tracked areas to keep melange')
    parser.add_argument('--check', action='store_true', default=False,
                        help='Setup command,  check masks, but do not run')
    parser.add_argument('--template', type=str, default='mosaic.template.yaml',
                        help='template that defines mosaic')
    parser.add_argument('--noLabel', action='store_true', default=False,
                        help='No processing source label')
    parser.add_argument('--metaOnly', action='store_true', default=False,
                        help='Rebuild shapefiles and metadata only; '
                        'skip masking, interpolation, and tif/vrt generation')
    parser.add_argument('--useSquint', action='store_true', default=False,
                        help='Pass -useSquint to mosaic3d (overrides template/project.yaml)')
    parser.add_argument('--useAzIonosphere', action='store_true', default=False,
                        help='Pass -useAzIonosphere to mosaic3d (overrides '
                             'template/project.yaml applyAzimuthIonosphereCorrection)')
    parser.add_argument('--timeThresh', type=float, default=None,
                        help='Max days between crossing-orbit offset pairs (mosaic3d default 12)')
    parser.add_argument('--jointMaxSigmaPhase', type=float, default=None,
                        help='Reject crossing-PHASE pixels with sigmaWorst*sqrt(n) > this '
                             '(m/yr); overrides template and project.yaml [mosaic3d default 35]')
    parser.add_argument('--jointMaxSigmaRange', type=float, default=None,
                        help='Same for crossing-RANGE; overrides template and project.yaml '
                             '[mosaic3d default 100]')
    parser.add_argument('--legacyPairPhase', action='store_true', default=False,
                        help='Use the ORIGINAL pairwise crossing-phase solver')
    parser.add_argument('--legacyPairRange', action='store_true', default=False,
                        help='Use the ORIGINAL pairwise crossing-range solver')
    parser.add_argument('--timePhaseThresh', type=float, default=None,
                        help='Max days between crossing-orbit phase pairs (mosaic3d default 548)')
    parser.add_argument('--mosaicsSetupFile', type=str,
                        default='mosaicsSetup.yaml',
                        help='yaml with seasonal and mask information')
    parser.add_argument('--nThreads', type=int, default=24,
                        help='max number of parallel sector threads (default 24)')
    parser.add_argument('--listTemplateKeys', action='store_true', default=False,
                        help='Print all keys understood in mosaic.template.yaml, '
                        'with defaults, and exit')
    #
    args = parser.parse_args()
    #
    if args.listTemplateKeys:
        printTemplateKeys()
        sys.exit(0)
    #
    mosaicsSetup = mosf.readYaml(args.mosaicsSetupFile, returnEmpty=True)
    if args.landsatPath is not None or 'landsatPath' not in mosaicsSetup:
        mosaicsSetup['landsatPath'] = args.landsatPath
    if args.mosaicMask is not None or 'mosaicMask' not in mosaicsSetup:
        mosaicsSetup['mosaicMask'] = args.mosaicMask
    # OVerride with commandline
    inputFile = None
    if args.inputFile is not None:
        inputFile = args.inputFile
    # else default if not defined in file
    # multiYear dates not adjusted
    firstDate = adjustFirstDate(datetime.strptime(args.firstdate, "%Y-%m-%d"),
                                args.interval, mosaicsSetup['seasonData'])
    endDate = datetime.strptime(args.lastdate, "%Y-%m-%d")
    myArgs = {'firstDate': firstDate, 'lastDate': endDate,
              'mosaicsSetup': mosaicsSetup,
              'interval': args.interval, 'noReprocess': args.noReprocess,
              'inputFile': inputFile,
              'landsatPath': mosaicsSetup['landsatPath'],
              'check': args.check, 'template': args.template,
              'outputMask': args.outputMask,
              'keepFast': args.keepFast, 'baseFlags': args.baseFlags,
              'noLandsat': args.noLandsat, 'fitType': args.LSFitType,
              'noLabel': args.noLabel, 'metaOnly': args.metaOnly,
              'useSquint': args.useSquint,
              'useAzIonosphere': args.useAzIonosphere,
              'timeThresh': args.timeThresh,
              'timePhaseThresh': args.timePhaseThresh,
              'jointMaxSigmaPhase': args.jointMaxSigmaPhase,
              'jointMaxSigmaRange': args.jointMaxSigmaRange,
              'legacyPairPhase': args.legacyPairPhase,
              'legacyPairRange': args.legacyPairRange,
              'nThreads': args.nThreads}
    return myArgs


def adjustFirstDate(firstdate, interval, seasonData):
    ''' adjust first date to match predefined ranges for interval '''
    if 'multiYear' in interval:
        return firstdate
    #

    firstdateOrig = firstdate
    # if region == 'quarterlyJFM':
    #     customDates = [0, 1, 1, 1, 4, 4, 4, 7, 7, 7, 10, 10, 10]
    #     customQuarterlyMonth = 2
    #     interval = 'quarterly'
    # else:
    #     customDates = [0, 12, 12, 3, 3, 3, 6, 6, 6, 9, 9, 9, 12]
    #     customQuarterlyMonth = 0
    # quarterlyDates ={'greenland': [0, 12, 12, 3, 3, 3, 6, 6, 6, 9, 9, 9, 12],
    #                   'amundsen': [0, 1, 1, 1, 4, 4, 4, 7, 7, 7, 10, 10, 10],
    #                   'taku': [0, 12, 12, 3, 3, 3, 6, 6, 6, 9, 9, 9, 12],
    #                   'custom': customDates}
    # quarterlyMonths = {'greenland': 2, 'amundsen': 0, 'taku': 2,
    #                    'custom': customQuarterlyMonth}
    # annualDates = {'greenland': datetime(firstdate.year, 12, 1),
    #                'amundsen': datetime(firstdate.year, 1, 1),
    #                'custom': datetime(firstdate.year, 1, 1),
    #                'taku': datetime(firstdate.year, 1, 1)}
    #
    if interval == 's1cycle' or interval == 's1-12day':
        if interval == 's1-12day':
            myDates = mosf.standardDates(nDays=12)
        else:
            myDates = mosf.standardDates()
        myDates.reverse()
        for myDate in myDates:  # look through possible dates until first valid
            if myDate['date1'] <= firstdate:
                break
        firstdate = myDate['date1']
    elif interval == 'monthly':
        # force to start at beginning of month
        firstdate = firstdate.replace(day=1)
    elif interval == 'quarterly':
        # start at beginning of quarter that contains firstdate
        month = seasonData['quarterlyDates'][firstdate.month]
        year = firstdate.year
        if firstdate.month <= seasonData['quarterlyMonth']:
            year = year-1
            month = seasonData['quarterlyMonth'] + 10
        firstdate = datetime(year, month, 1)
    elif interval == 'annual':
        firstdate = datetime(firstdate.year, seasonData['annualFirstMonth'], 1)
    #
    # warning that date has been ajusted
    if firstdate != firstdateOrig:
        u.mywarning(
                f'Adjusting first date from '
                f'{firstdateOrig.strftime("%Y-%m-%d")} to '
                f'{firstdate.strftime("%Y-%m-%d")} to be consisent with '
                f'{interval} ranges')
    return firstdate


def getLists(landsatPath, fitType='custom.mask'):
    ''' get mask and list files, along with masks '''
    # myPath = '/Volumes/insar7/ian/LANDSAT/Greenland'
    listFiles = u.dols(f'ls {landsatPath}/listfiles/Listfile.*.{fitType}')
    return listFiles


def incrementDate(myDate, interval):
    '''  increment by a s1 cycle, month, quarter, or year'''
    if interval == 's1cycle':
        date6 = datetime(2016, 9, 20)  # Transition from 12 to 6 day
        if myDate > date6 and myDate <= datetime(2021, 12, 24):
            dT = timedelta(days=6)
        else:
            dT = timedelta(days=12)
        return myDate + dT
    if interval == 's1-12day':
        return myDate + timedelta(days=12)
    # All other intervales
    dT = {'monthly': 32, 'quarterly': 93, 'quarterlyJFM': 93,
          'annual': 367}[interval]
    myDate = myDate + timedelta(days=dT)
    # force to start of month
    myDate = myDate.replace(day=1)
    return myDate


def findLastDate(firstDate, myArgs):
    ''' Find last date given firstDate '''
    # Multiyear is a single time span so return last date
    if 'multiYear' in myArgs['interval']:
        return myArgs['lastDate'], myArgs['lastDate'].year
    #
    # go back one day to get end of month
    endDate = incrementDate(firstDate, myArgs['interval']) - timedelta(days=1)
    # back date to first of month,  then subtract 1 day to get to end of month
    tmpDate = firstDate + (endDate-firstDate)*.5
    year = tmpDate.year
    return endDate, year


def createMergedList(listFiles, firstDate, lastDate):
    ''' merge landsat lists  lists'''
    shortList = []
    for listFile in listFiles:
        try:
            listYear = int(listFile.split('.')[-3])
            if listYear >= firstDate.year and listYear <= lastDate.year:
                shortList.append(listFile)
        except Exception:
            u.myerror('createMergedList: Problem parsing list file')
    #
    # now build list
    mergedList = f'mergedList.{firstDate.year}-{lastDate.year}'
    fp = open(mergedList, 'w')
    for listFile in shortList:
        fpIn = open(listFile, 'r')
        for line in fpIn:
            if '&' not in line:
                print(line, file=fp, end='')
        fpIn.close()
    print('&', file=fp)
    fp.close()
    return mergedList


def getNewMask(interval, firstDate, lastDate, mosaicsSetup):
    '''
    Fill in mask templates with the season and year
    '''
    maskPath = mosaicsSetup['intervalMasks'][interval]
    whichMask = mosaicsSetup['whichMask'][interval]
    #
    month = (firstDate + (lastDate - firstDate) * 0.5).month
    year = (firstDate + (lastDate - firstDate) * 0.5).year
    # update mask template with relevant info
    print(maskPath)
    maskFile = maskPath.replace('YY', f'{year-2000}').replace(
        'SS', whichMask[month-1])
    #
    return maskFile


def getMosaicMask(mosaicsSetup, interval, firstDate, lastDate, keepFast):
    ''' Return appropriate mask file '''
    if keepFast:
        print('keepFast flag set so no mask applied')
        return None
    # Use explicit version if given
    if mosaicsSetup['mosaicMask'] is not None:
        return mosaicsSetup['mosaicMask']
    # Select seasonal mask
    if 'intervalMasks' not in mosaicsSetup:
        return None
    maskFile = getNewMask(interval, firstDate, lastDate, mosaicsSetup)
    if not os.path.exists(maskFile):
        u.myerror(f'getMosaicMask: maskfile {maskFile} does not exist')
    print(maskFile)
    return maskFile


def makeCommand(firstDate, lastDate, mergedList, mosaicMaskFile, myArgs):
    ''' setup and return command '''
    #
    outputMaskArg, templateArg, lsArg, baseFlagsArg, keepFastFlag, \
        mosaicMaskArg, noReprocessFlag, noTSXFlag, noLabelFlag, \
        metaOnlyFlag, useSquintFlag, useAzIonFlag, timeThreshArg, timePhaseThreshArg, \
        jointSigPhaseArg, jointSigRangeArg, legacyPhaseFlag, legacyRangeFlag, \
        nThreadsArg = [''] * 19
    #
    noReprocessFlag = {False: '', True: '--noReprocess'}[myArgs["noReprocess"]]
    #
    # output mask
    if myArgs['outputMask'] is not None:
        outputMaskArg = '--outputMask {myArgs["outputMask"]}'
    #
    # mosaicMask
    if mosaicMaskFile is not None:
        mosaicMaskArg = f'--mosaicMask {mosaicMaskFile} '
    #
    # template file
    if myArgs["template"] is not None:
        templateArg = f'--template {myArgs["template"]} '
    #
    # radar inputs
    inputFileArg = ' '
    if myArgs['outputMask'] is not None:
        inputFileArg = f' --inputFile {myArgs["inputFile"]}'
    #
    # landsat inputs
    if mergedList is not None:
        lsArg = f'--lsFile {mergedList} '
    #
    # baseFlags
    if myArgs["baseFlags"] is not None:
        baseFlagsArg = f'--baseFlags \"{myArgs["baseFlags"]}\" '
    # keepFast uses no cull for melanges
    if myArgs["keepFast"]:
        keepFastFlag = '--noCull '
    if myArgs["noLabel"]:
        noLabelFlag = '--noLabel '
    if myArgs["metaOnly"]:
        metaOnlyFlag = '--metaOnly '
    if myArgs.get('useSquint', False):
        useSquintFlag = '--useSquint '
    if myArgs.get('useAzIonosphere', False):
        useAzIonFlag = '--useAzIonosphere '
    if myArgs.get('timeThresh') is not None:
        timeThreshArg = f'--timeThresh {myArgs["timeThresh"]} '
    if myArgs.get('timePhaseThresh') is not None:
        timePhaseThreshArg = f'--timePhaseThresh {myArgs["timePhaseThresh"]} '
    if myArgs.get('jointMaxSigmaPhase') is not None:
        jointSigPhaseArg = f'--jointMaxSigmaPhase {myArgs["jointMaxSigmaPhase"]} '
    if myArgs.get('jointMaxSigmaRange') is not None:
        jointSigRangeArg = f'--jointMaxSigmaRange {myArgs["jointMaxSigmaRange"]} '
    if myArgs.get('legacyPairPhase', False):
        legacyPhaseFlag = '--legacyPairPhase '
    if myArgs.get('legacyPairRange', False):
        legacyRangeFlag = '--legacyPairRange '
    if myArgs.get('nThreads') is not None:
        nThreadsArg = f'--nThreads {myArgs["nThreads"]} '
    # TSX excluded from single sycle data
    if myArgs["interval"] == 's1cycle' or myArgs["interval"] == 's1-12day':
        noTSXFlag = ' --noTSX'

    command = f'setupquarters.py ' \
        f'{templateArg} ' \
        f'--firstdate {firstDate.strftime("%Y-%m-%d")} ' \
        f' --lastdate {lastDate.strftime("%Y-%m-%d")} ' \
        f'{noReprocessFlag} {keepFastFlag} {noTSXFlag} {noLabelFlag} ' \
        f'{metaOnlyFlag}' \
        f'{useSquintFlag}' \
        f'{useAzIonFlag}' \
        f'{timeThreshArg}' \
        f'{timePhaseThreshArg}' \
        f'{jointSigPhaseArg}{jointSigRangeArg}' \
        f'{legacyPhaseFlag}{legacyRangeFlag}' \
        f'{baseFlagsArg}' \
        f'{nThreadsArg}' \
        f'{outputMaskArg} '  \
        f'{mosaicMaskArg} ' \
        f'{inputFileArg} {lsArg} '

    print(command)
    return command


def main():
    ''' Analyze velocity series with respect to tiepoints '''
    # get args
    myArgs = makemosaicArgs()
    #
    if myArgs['mosaicsSetup']['landsatPath'] is not None and \
            not myArgs['noLandsat'] and not myArgs['metaOnly']:
        listFiles = getLists(myArgs['landsatPath'], fitType=myArgs['fitType'])
        mergedList = createMergedList(listFiles, myArgs['firstDate'],
                                      myArgs['lastDate'])
    else:
        mergedList = None
    print(f'mergedList = {mergedList}')
    #
    currentFirstDate = myArgs['firstDate']
    currentLastDate, year = findLastDate(currentFirstDate, myArgs)
    print(currentFirstDate, currentLastDate)
    # loop to produce products as defined by date range
    while currentLastDate <= myArgs['lastDate']:
        # get mask
        if not myArgs['metaOnly']:
            mosaicMaskFile = getMosaicMask(myArgs['mosaicsSetup'],
                                           myArgs['interval'],
                                           currentFirstDate,
                                           currentLastDate,
                                           myArgs['keepFast'])
        else:
            mosaicMaskFile = None
        # setup command
        command = makeCommand(currentFirstDate, currentLastDate, mergedList,
                              mosaicMaskFile, myArgs)
        #
        if not myArgs['check']:
            returnCode = call(command, shell=True)  # , executable='/bin/csh'
            # setupquarters.py exits nonzero if any sector's mosaic3d failed;
            # stop here rather than march on and build more mosaics that may
            # ship stale results.
            if returnCode != 0:
                print(f'\n\t\033[1;31m *** makemosaic: setupquarters.py '
                      f'failed (exit {returnCode}) for '
                      f'{currentFirstDate.strftime("%Y-%m-%d")} to '
                      f'{currentLastDate.strftime("%Y-%m-%d")}; stopping before '
                      f'building further mosaics. Check the Vel-* io/ logs. '
                      f'*** \033[0m\n')
                sys.exit(1)
        # update dates
        if 'multiYear' in myArgs['interval']:
            break  # multi Year one off product
        currentFirstDate = incrementDate(currentFirstDate, myArgs['interval'])
        currentLastDate, year = findLastDate(currentFirstDate, myArgs)


if __name__ == "__main__":
    main()
