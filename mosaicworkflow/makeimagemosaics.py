#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 08:45:12 2020

@author: ian
"""

import argparse
import utilities as u
import os
from datetime import datetime, date, timedelta
from subprocess import call
import threading
import mosaicfunc as mosf

currentVersion = 4
subVersion = 0

# NISAR repeats every 12 days; cycle 1 started on this date.
NISAR_CYCLE_EPOCH = datetime(2025, 9, 23)


def nisarCycleDates(cycleNum):
    '''Return (date1, date2) for the given 1-based NISAR cycle number.'''
    date1 = NISAR_CYCLE_EPOCH + timedelta(days=(cycleNum - 1) * 12)
    date2 = date1 + timedelta(days=11)
    return date1, date2


def makemosaicArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\033[1mFront end to run image mosaics '
        'using setupimagemosaic.py \033[0m',
        epilog='Notes:\nPart of the mosaicworkflow package.', allow_abbrev='False')
    parser.add_argument('--firstdate', type=str, default='2015-01-01',
                        help='Use first dates >= first date [2015-01-01]')
    defaultEnd = date.today().strftime('%Y-%m-%d')
    parser.add_argument('--lastdate', type=str, default=defaultEnd,
                        help=f'Use last dates <= lastdate [{defaultEnd}]')
    parser.add_argument('--noReprocess', action='store_true', default=False,
                        help='Just reformat data in release dir')
    parser.add_argument('--noMosaic', action='store_true', default=False,
                        help='Do not rerun mosaics - reassemble from pieces')
    parser.add_argument('--noAdjustDate', action='store_true', default=False,
                        help='Force dates to start on first date, ow dates'
                        'will use standard invervals based on 12 then 6 days'
                        'from 2015-01-01')
    parser.add_argument('--nDays', type=int, default=None,
                        help='Use nDays interval starting at firstdate [6/12]')
    parser.add_argument('--noCheckEndDate', action='store_true', default=False,
                        help='Do not clip end date based on dates [False]')
    parser.add_argument('--reset', action='store_true', default=False,
                        help='Redo existing, ow only do new products')
    parser.add_argument('--noClean', action='store_true', default=False,
                        help='Save the tiles in the the pieces directory')
    parser.add_argument('--setupOnly', action='store_true', default=False,
                        help='Do all steps but run - can be used to clean')
    parser.add_argument('--calibrate', action='store_true', default=False,
                        help='Make an S1 calibrated image product')
    parser.add_argument('--corr', action='store_true', default=False,
                        help='Make a NISAR correlation (coherence) mosaic')
    parser.add_argument('--check', action='store_true', default=False,
                        help='Setup command,  check masks, but do not run')
    parser.add_argument('--ascendingOnly', action='store_true', default=False,
                        help='Only include ascending data')
    parser.add_argument('--descendingOnly', action='store_true', default=False,
                        help='Only include descending data')
    parser.add_argument('--prefix', type=str, default='GL_S1bks_mosaic',
                        help='Prefix for product name')
    parser.add_argument('template_pos', type=str, nargs='?', default=None,
                        help='template file (positional form)')
    parser.add_argument('--template', type=str, default='mosaic.template.yaml',
                        help='template that defines mosaic [mosaic.template.yaml]')
    parser.add_argument('--removePad', type=int, default=5,
                        help='Remove first and last n cols')
    parser.add_argument('--averageAll', action='store_true', default=False,
                        help='Create a single mosaic averaging all scenes '
                        'from --firstdate to --lastdate (no date splitting)')
    parser.add_argument('--min', action='store_true', default=False,
                        help='Mosaic = pixel-wise minimum across all inputs')
    parser.add_argument('--max', action='store_true', default=False,
                        help='Mosaic = pixel-wise maximum across all inputs')
    args = parser.parse_args()
    template = args.template_pos if args.template_pos is not None else args.template
    _, _, DBPaths, _, prefix, _, removePad, _, averageAllTemplate = mosf.mosaicTempArgs(template)
    if args.removePad != 5:
        removePad = args.removePad
    if prefix is None:
        prefix = args.prefix
    # geoMosaicMode: CLI flag wins over template key
    templateDict = mosf.readTemplate(template)
    geoMosaicMode = templateDict.get('geoMosaicMode', None)
    if args.min:
        geoMosaicMode = 'min'
    elif args.max:
        geoMosaicMode = 'max'
    # nisarCycle mode: generate products by 12-day cycle number
    nisarCycle = bool(templateDict.get('nisarCycle', False))
    #
    firstDate = datetime.strptime(args.firstdate, "%Y-%m-%d")
    endDate = datetime.strptime(args.lastdate, "%Y-%m-%d")
    # CLI flag wins; template value is the fallback default
    # min/max mode always produces a single mosaic over the full date range
    averageAll = args.averageAll or averageAllTemplate or geoMosaicMode in ('min', 'max')
    # Adjust dates to standard intervals (skip for averageAll or nisarCycle —
    # in cycle mode, getProductList generates exact NISAR cycle boundaries)
    if not args.noAdjustDate and not averageAll and not nisarCycle:
        firstDate = adjustFirstDate(firstDate, nisar=args.corr)
    if args.descendingOnly and args.ascendingOnly:
        u.myerror('descendingOnly and ascendingOnly flags mutually exclusive')
    myArgs = {'firstDate': firstDate, 'lastDate': endDate,
              'noReprocess': args.noReprocess, 'check': args.check,
              'noMosaic': args.noMosaic,
              'calibrate': args.calibrate, 'corr': args.corr, 'noClean': args.noClean,
              'setupOnly': args.setupOnly, 'reset': args.reset,
              'templateFile': template, 'nDays': args.nDays,
              'noCheckEndDate': args.noCheckEndDate,
              'noAdjustDate': args.noAdjustDate,
              'averageAll': averageAll,
              'ascendingOnly': args.ascendingOnly,
              'descendingOnly': args.descendingOnly,
              'removePad': removePad,
              'prefix': prefix, 'dateDBs': DBPaths,
              'geoMosaicMode': geoMosaicMode,
              'nisarCycle': nisarCycle}
    return myArgs


def adjustFirstDate(firstdate, nisar=False):
    '''
    Adjust the first date to the standard dates, so that it precedes
    the standard date  (e.g., s1 fd s2 s3) -> fd=s2
    Parameters
    ----------
    firstdate : datetime
        Nominal first date.
    nisar : bool, optional
        If True, never resume 6-day cycles (NISAR repeat is always 12 days).
    Returns
    -------
    myDate : date time
        Adjusted fir date..

    '''
    resume6Day = datetime(2100, 1, 1) if nisar else datetime(2025, 4, 1)
    myDates = mosf.standardDates(resume6DayDates=resume6Day)
    myDates.reverse()
    for myDate in myDates:
        if myDate['date1'] <= firstdate:
            return myDate['date1']


def getProductList(myArgs):
    '''
    Generate list of products (date1, interval) to process.

    In nisarCycle mode each entry also carries a 'cycle' key (int) so the
    product can be labelled with the NISAR cycle number.

    Parameters
    ----------
    myArgs : dict
        Arguments that control mosaics.

    Returns
    -------
    myProducts : list of dicts  {'date1', 'date2'[, 'cycle']}
    '''
    if myArgs.get('averageAll') or myArgs.get('geoMosaicMode') in ('min', 'max'):
        return [{'date1': myArgs['firstDate'], 'date2': myArgs['lastDate']}]

    firstDate = myArgs['firstDate']
    lastDate = myArgs['lastDate']
    print(f'Date Range = {firstDate.strftime("%Y-%m-%d")}'
          f' {lastDate.strftime("%Y-%m-%d")}')

    if myArgs.get('nisarCycle'):
        # Generate products aligned to NISAR 12-day repeat cycles.
        # Find the first cycle whose start is >= firstDate (backing up one
        # cycle to avoid missing a cycle that straddles the boundary).
        delta0 = max(0, (firstDate - NISAR_CYCLE_EPOCH).days)
        startCycle = max(1, delta0 // 12)  # one cycle before to be safe
        myDates = []
        cycleNum = startCycle
        while True:
            date1, date2 = nisarCycleDates(cycleNum)
            if date1 > lastDate:
                break
            if date1 >= firstDate and date2 <= lastDate:
                myDates.append({'date1': date1, 'date2': date2,
                                'cycle': cycleNum})
            cycleNum += 1
        return myDates

    # --- standard S1 date-range mode ---
    anchorDate = None
    if myArgs['noAdjustDate']:
        anchorDate = firstDate
    resume6Day = datetime(2100, 1, 1) if myArgs.get('corr') \
        else datetime(2025, 4, 1)
    sDateRanges = mosf.standardDates(nDays=myArgs['nDays'],
                                     firstDate=anchorDate,
                                     resume6DayDates=resume6Day)
    myDates = []
    for dateRange in sDateRanges:
        if dateRange['date1'] < firstDate or dateRange['date2'] > lastDate:
            continue
        myDates.append(dateRange)
    return myDates


def makeImageMosaic(date1, date2, myArgs, cycleNum=None):
    '''
    Spawn the mosaic command.
    Parameters
    ----------
    date1, date2 : datetime
        First and second date.
    myArgs : dict
        Arguments that control mosaics.
    cycleNum : int or None
        When set (nisarCycle mode), appends --cycleTag cycle<NNN> to the
        setupimagemosaic command so the product is labelled accordingly.
    Returns
    -------
    None.
    '''
    nDays = (date2-date1).days + 1
    command = f'setupimagemosaic.py --firstdate {date1.strftime("%Y-%m-%d")}' \
        f' --nDays {nDays} --template {myArgs["templateFile"]} ' \
        f'--removePad {myArgs["removePad"]}'
    for key in ['calibrate', 'corr', 'noReprocess', 'noMosaic', 'noClean',
                'setupOnly', 'descendingOnly', 'ascendingOnly', 'averageAll']:
        if myArgs[key]:
            command += f' --{key}'
    if myArgs.get('geoMosaicMode') == 'min':
        command += ' --min'
    elif myArgs.get('geoMosaicMode') == 'max':
        command += ' --max'
    if cycleNum is not None:
        command += f' --cycleTag cycle{cycleNum:03d}'
    command += f' --prefix {myArgs["prefix"]}'
    print(f'\n{command}')
    call(command, shell=True)  # , executable='/bin/csh')


def makeProdName(prefix, date1, date2):
    ''' Make product name  -- cross ref against same name prog in
    setupimagemosaic.py'''
    prodName = f'{prefix}_{date1.strftime("%d%b%y")}_' \
        f'{date2.strftime("%d%b%y")}_*_RESm_v{currentVersion:02d}.{subVersion}'
    return prodName


def checkProduct(date1, date2, myArgs, posting, cycleNum=None):
    '''
    Checks if a product has already been produced. Verifies that tifs exist.
    Parameters
    ----------
    date1 : datetime
        First date.
    date2 : datetime
        Second date.
    myArgs : dict
        parameters defining mosaic creation.
    cycleNum : int or None
        When set (nisarCycle mode), the cycle tag is folded into the prefix
        when checking for the product directory and files.
    Returns
    -------
    exists : bool — True if directory and image products both exist.
    '''
    if myArgs.get('corr'):
        calString = 'corr'
    elif myArgs['calibrate']:
        calString = 'calibrated'
    else:
        calString = 'uncalibrated'
    prefix = myArgs["prefix"]
    if cycleNum is not None:
        prefix = f'{prefix}_cycle{cycleNum:03d}'
    dirName = f'{prefix}_{calString}.' \
        f'{date1.strftime("%Y-%m-%d")}.{date2.strftime("%Y-%m-%d")}'
    prodName = makeProdName(prefix, date1, date2)
    if not os.path.exists(dirName):
        return False  # No product directory so return false
    # check image products exist
    if myArgs.get('corr'):
        suffixes = ['corr']
    elif myArgs['calibrate']:
        suffixes = ['gamma0', 'sigma0']
    else:
        suffixes = ['image']
    for suffix in suffixes:
        tifName = prodName.replace("*", suffix).replace("RES", str(posting))
        myTiff = f'{dirName}/release/{tifName}.tif'
        if not os.path.exists(myTiff):
            return False  # No image product so return false
    return True


def setupProducts(myProds, myArgs, posting):
    '''
    Setup individual mosaics to run via "setupimagemosaic.py"
    Parameters
    ----------
    myProds : list
        list of dicts {'date1': datetime, 'date2': datetime[, 'cycle': int]}.
    myArgs : dict
        Arguments that control mosaics.
    posting : int
        Pixel posting (metres) used to verify existing products.
    Returns
    -------
    threads : list of threads.
    '''
    threads = []
    for myProd in myProds:
        cycleNum = myProd.get('cycle')
        prodExists = checkProduct(myProd['date1'], myProd['date2'], myArgs,
                                  posting, cycleNum=cycleNum)
        # Put product in queue if doesn't exist or if reset mode
        if not prodExists or myArgs['reset'] or myArgs['noReprocess']:
            threads.append(threading.Thread(target=makeImageMosaic,
                                            args=[myProd['date1'],
                                                  myProd['date2'],
                                                  myArgs],
                                            kwargs={'cycleNum': cycleNum}))
    return threads


def checkLastDateAgainsDB(myArgs):
    '''
    Clip dates to date range. Issue warning if out of range.

    Parameters
    ----------
    myArgs : dict
        Parsed command line args.
    '''
    myDB = mosf.sarDB()
    myDB.readDB(myArgs['dateDBs'])
    if myArgs['lastDate'] > myDB.maxDate:
        u.mywarning(f'Last date in DB is {myDB.maxDate} but last date '
                    f'requested is {myArgs["lastDate"]}. \nClipping last'
                    'to last DB date: Update DB for later dates')
        myArgs['lastDate'] = myDB.maxDate


def main():
    ''' Produce image mosaics '''
    # get args
    myArgs = makemosaicArgs()
    posting = mosf.readTemplate(myArgs['templateFile'])['dx']
    if not myArgs['noCheckEndDate']:
        checkLastDateAgainsDB(myArgs)
    myProds = getProductList(myArgs)
    threads = setupProducts(myProds, myArgs, posting)
    maxThreads = 1
    u.runMyThreads(threads, maxThreads, 'geomosaic')


if __name__ == "__main__":
    main()
