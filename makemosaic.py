#!/usr/bin/env python3
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
    parser.add_argument('--interval', type=str,
                        choices=['s1cycle', 's1-12day', 'monthly', 'quarterly',
                                 'annual', 'multiYear'], default='monthly',
                        help='time interval spanned')
    parser.add_argument('--region', type=str,
                        choices=['greenland', 'amundsen', 'taku'],
                        default='greenland',
                        help='region code if not customRegion')
    parser.add_argument('--customRegion', type=str,
                        default=None,
                        help='yml file with custom region specification')
    parser.add_argument('--noReprocess', action='store_true', default=False,
                        help='Just reformat data')
    parser.add_argument('--noLandsat', action='store_true', default=False,
                        help='Do not include landsat')
    parser.add_argument('--LSFitType', type=str,
                        choices=['custom.mask', 'custom', 'static',
                                 'static.mask'],
                        default='custom.mask',
                        help='Landsat fit type')
    parser.add_argument('--mask', type=str, default=None,
                        help='Uses non-default mask')
    parser.add_argument('--baseFlags', type=str, default=None,
                        help='Override default flags for the mosaicker')
    parser.add_argument('--keepFast', action='store_true', default=False,
                        help='Do not mask fast tracked areas to keep melange')
    parser.add_argument('--check', action='store_true', default=False,
                        help='Setup command,  check masks, but do not run')
    args = parser.parse_args()
    # multiYear dates not adjusted
    firstDate = adjustFirstDate(datetime.strptime(args.firstdate, "%Y-%m-%d"),
                                args.interval, args.region)
    endDate = datetime.strptime(args.lastdate, "%Y-%m-%d")
    if args.customRegion is not None:
        args.region = 'custom'
    myArgs = {'firstDate': firstDate, 'lastDate': endDate,
              'interval': args.interval, 'noReprocess': args.noReprocess,
              'check': args.check, 'region': args.region,
              'customRegion': args.customRegion, 'mask': args.mask,
              'keepFast': args.keepFast, 'baseFlags': args.baseFlags,
              'noLandsat': args.noLandsat, 'fitType': args.LSFitType}
    return myArgs


def adjustFirstDate(firstdate, interval, region):
    ''' adjust first date to match predefined ranges for interval '''
    if 'multiYear' in interval:
        return firstdate
    firstdateOrig = firstdate
    quarterlyDates = {'greenland': [0, 12, 12, 3, 3, 3, 6, 6, 6, 9, 9, 9, 12],
                      'amundsen': [0, 1, 1, 1, 4, 4, 4, 7, 7, 7, 10, 10, 10],
                      'taku': [0, 1, 1, 1, 4, 4, 4, 7, 7, 7, 10, 10, 10],
                      'custom': [0, 1, 1, 1, 4, 4, 4, 7, 7, 7, 10, 10, 10]}
    quarterlyMonths = {'greenland': 2, 'amundsen': 0, 'taku': 0}
    annualDates = {'greenland': datetime(firstdate.year, 12, 1),
                   'amundsen': datetime(firstdate.year, 1, 1),
                   'custom': datetime(firstdate.year, 1, 1),
                   'taku': datetime(firstdate.year, 1, 1)}
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
        month = quarterlyDates[region][firstdate.month]
        year = firstdate.year
        if firstdate.month <= quarterlyMonths[region]:
            year = year-1
        firstdate = datetime(year, month, 1)
    elif interval == 'annual':
        firstdate = annualDates[region]
    #
    # warning that date has been ajusted
    if firstdate != firstdateOrig:
        u.mywarning(
                f'Adjusting first date from '
                f'{firstdateOrig.strftime("%Y-%m-%d")} to '
                f'{firstdate.strftime("%Y-%m-%d")} to be consisent with '
                f'{interval} ranges')
    return firstdate


def getLists(fitType='custom.mask'):
    ''' get mask and list files, along with masks '''
    myPath = '/Volumes/insar7/ian/LANDSAT/Greenland'
    listFiles = u.dols(f'ls {myPath}/listfiles/Listfile.*.{fitType}')
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
    dT = {'monthly': 32, 'quarterly': 93, 'annual': 367}[interval]
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


def getNewMask(interval, firstDate, lastDate):
    # template for type
    maskPath = '/Users/ian/greenlandmask'
    maskPath = {'quarterly': f'{maskPath}/quartermasks/greenlandmask20YYSS',
                'monthly': f'{maskPath}/output20YY/greenlandmaskSS20YY',
                's1cycle': f'{maskPath}/output20YY/greenlandmaskSS20YY',
                's1-12day': f'{maskPath}/output20YY/greenlandmaskSS20YY',
                'annual':
                    f'{maskPath}/annualmasks/greenlandmask20YYA',
                'multiYear':
                    f'{maskPath}/annualmasks/greenlandmask20YYA'}[interval]
    # ids for a given interval type
    whichMask = {'quarterly': ['Q1', 'Q1', 'Q2', 'Q2', 'Q2', 'Q3', 'Q3', 'Q3',
                               'Q4', 'Q4', 'Q4', 'Q1'],
                 'monthly': ['lw', 'lw', 'lw', 'sp', 'sp', 'ms', 'ms', 'ms',
                             'ls', 'ls', 'ew', 'ew'],
                 's1cycle': ['lw', 'lw', 'lw', 'sp', 'sp', 'ms', 'ms', 'ms',
                             'ls', 'ls', 'ew', 'ew'],
                 's1-12day': ['lw', 'lw', 'lw', 'sp', 'sp', 'ms', 'ms', 'ms',
                              'ls', 'ls', 'ew', 'ew'],
                 'annual': ['A']*12,
                 'multiYear': ['A']*12}[interval]
    #
    month = (firstDate+(lastDate-firstDate)*0.5).month
    year = (firstDate+(lastDate-firstDate)*0.5).year
    # update mask template with relevant info
    maskFile = maskPath.replace('YY', f'{year-2000}').replace(
        'SS', whichMask[month-1])
    #
    return maskFile


def getGreenlandMask(interval, firstDate, lastDate, keepFast):
    ''' Return appropriate mask file '''
    if keepFast:  # This uses the melange mask
        return '/Users/ian/greenlandmask/MelangeMask/melangemask'
    maskPath = '/Users/ian/greenlandmask'
    mqmasks = {2014: f'{maskPath}/Winter201415/greenlandmask201415',
               2015: f'{maskPath}/Winter201516/greenlandmask201516',
               2016: f'{maskPath}/Winter201617/greenlandmask201617'}
    annmasks = {2015: f'{maskPath}/annualmasks/greenlandmask2015A',
                2016: f'{maskPath}/annualmasks/greenlandmask2016A'}
    oldMasks = {'s1cycle': mqmasks, 's1-12day': mqmasks, 'monthly': mqmasks,
                'quarterly': mqmasks, 'annual': annmasks}
    #
    if firstDate < datetime(2016, 12, 1):
        if interval in ['s1cycle', 's1-12day', 'quarterly', 'monthly']:
            myYear = max((firstDate-timedelta(days=182)).year, 2014)
        else:
            myYear = (firstDate+timedelta(days=182)).year
        maskFile = oldMasks[interval][myYear]
    else:
        maskFile = getNewMask(interval, firstDate, lastDate)
    if not os.path.exists(maskFile):
        u.myerror(f'getMask: maskfile {maskFile} does not exist')
    print(maskFile)
    return maskFile


def makeCommand(firstDate, lastDate, noReprocess, listFile, maskFile, region,
                customRegion, interval, keepFast, baseFlags):
    ''' setup and return command '''
    #
    if customRegion is not None:
        regionUsed = f'customRegion={customRegion}'
    else:
        regionUsed = region
    command = f'setupquarters.py -{regionUsed} -firstdate='
    command += firstDate.strftime('%Y-%m-%d')
    command += ['', ' -noReprocess '][noReprocess]
    command += f' -lastdate={lastDate.strftime("%Y-%m-%d")} '
    if baseFlags is not None:
        command += f' -baseFlags="{baseFlags}" '
    if keepFast:
        command += ' -noCull '
    if region == 'greenland' and interval != 's1cycle' and \
            listFile is not None:
        command += f'-lsFile={listFile}'
    if interval == 's1cycle' or interval == 's1-12day':
        command += ' -noTSX'
    command += ' -inputFile=track-all/inputFile '
    if maskFile is not None:
        command += f'-shelfMask={maskFile}'
    print(command)
    return command


def main():
    ''' Analyze velocity series with respect to tiepoints '''
    # get args
    myArgs = makemosaicArgs()
    #
    print(myArgs['noLandsat'])
    if myArgs['region'] in ['greenland'] and not myArgs['noLandsat']:
        listFiles = getLists(fitType=myArgs['fitType'])
        mergedList = createMergedList(listFiles, myArgs['firstDate'],
                                      myArgs['lastDate'])
    else:
        mergedList = None
    print('Date range', myArgs['firstDate'], myArgs['lastDate'])
    print(myArgs)
    # exit()
    #
    currentFirstDate = myArgs['firstDate']
    currentLastDate, year = findLastDate(currentFirstDate, myArgs)
    print(currentFirstDate, currentLastDate)
    # loop to produce products as defined by date range
    while currentLastDate <= myArgs['lastDate']:
        # get mask
        if myArgs['region'] == 'greenland' and myArgs['mask'] is None:
            maskFile = getGreenlandMask(myArgs['interval'], currentFirstDate,
                                        currentLastDate, myArgs['keepFast'])
        else:
            maskFile = myArgs['mask']
        # print(currentFirstDate, currentLastDate, year, maskFile)
        # setup command
        command = makeCommand(currentFirstDate, currentLastDate,
                              myArgs['noReprocess'], mergedList, maskFile,
                              myArgs['region'], myArgs['customRegion'],
                              myArgs['interval'],
                              myArgs['keepFast'], myArgs['baseFlags'])
        if not myArgs['check']:
            call(command, shell=True, executable='/bin/csh')
        # update dates
        if 'multiYear' in myArgs['interval']:
            break  # multi Year one off product
        currentFirstDate = incrementDate(currentFirstDate, myArgs['interval'])
        currentLastDate, year = findLastDate(currentFirstDate, myArgs)


main()
