#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 08:45:12 2020

@author: ian
"""

import argparse
import utilities as u
import os
from datetime import datetime, timedelta, date
from subprocess import call
import mosaicfunc as mosf


def makemosaicArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\033[1mFront end to run velocity mosaics '
        'using setupquarters.py \033[0m',
        epilog='Notes:  ', allow_abbrev='False')
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
    parser.add_argument('--mosaicsSetupFile', type=str,
                        default='mosaicsSetup.yaml',
                        help='yaml with seasonal and mask information')
    #
    args = parser.parse_args()
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
              'noLandsat': args.noLandsat, 'fitType': args.LSFitType}
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
        mosaicMaskArg, noReprocessFlag, noTSXFlag = [''] * 8
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
    # TSX excluded from single sycle data
    if myArgs["interval"] == 's1cycle' or myArgs["interval"] == 's1-12day':
        noTSXFlag = ' --noTSX'

    command = f'setupquarters.py ' \
        f'{templateArg} ' \
        f'--firstdate {firstDate.strftime("%Y-%m-%d")} ' \
        f' --lastdate {lastDate.strftime("%Y-%m-%d")} ' \
        f'{noReprocessFlag} {keepFastFlag} {noTSXFlag} ' \
        f'{baseFlagsArg}' \
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
            not myArgs['noLandsat']:
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
        mosaicMaskFile = getMosaicMask(myArgs['mosaicsSetup'],
                                       myArgs['interval'],
                                       currentFirstDate,
                                       currentLastDate,
                                       myArgs['keepFast'])
        # setup command
        command = makeCommand(currentFirstDate, currentLastDate, mergedList,
                              mosaicMaskFile, myArgs)
        #
        if not myArgs['check']:
            call(command, shell=True)  # , executable='/bin/csh')
        # update dates
        if 'multiYear' in myArgs['interval']:
            break  # multi Year one off product
        currentFirstDate = incrementDate(currentFirstDate, myArgs['interval'])
        currentLastDate, year = findLastDate(currentFirstDate, myArgs)


if __name__ == "__main__":
    main()
