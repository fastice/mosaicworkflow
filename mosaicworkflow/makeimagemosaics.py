#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 08:45:12 2020

@author: ian
"""

import argparse
import utilities as u
import os
from datetime import datetime, date
from subprocess import call
import threading
import mosaicfunc as mosf

currentVersion = 4
subVersion = 0


def makemosaicArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\033[1mFront end to run image mosaics '
        'using setupimagemosaic.py \033[0m',
        epilog='Notes:  ', allow_abbrev='False')
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
    parser.add_argument('--check', action='store_true', default=False,
                        help='Setup command,  check masks, but do not run')
    parser.add_argument('--prefix', type=str, default='GL_S1bks_mosaic',
                        help='Prefix for product name')
    parser.add_argument('--template', type=str, default='mosaic.template.yaml',
                        help='template that defines mosaic')
    args = parser.parse_args()
    _, _, DBPaths, _, prefix = mosf.mosaicTempArgs(args.template)
    if prefix is None:
        prefix = args.prefix
    #
    firstDate = datetime.strptime(args.firstdate, "%Y-%m-%d")
    endDate = datetime.strptime(args.lastdate, "%Y-%m-%d")
    # Adjust dates to standard intervals
    if not args.noAdjustDate:
        firstDate = adjustFirstDate(firstDate)

    myArgs = {'firstDate': firstDate, 'lastDate': endDate,
              'noReprocess': args.noReprocess, 'check': args.check,
              'noMosaic': args.noMosaic,
              'calibrate': args.calibrate, 'noClean': args.noClean,
              'setupOnly': args.setupOnly, 'reset': args.reset,
              'templateFile': args.template, 'nDays': args.nDays,
              'noCheckEndDate': args.noCheckEndDate,
              'noAdjustDate': args.noAdjustDate,
              'prefix': prefix, 'dateDBs': DBPaths}
    print(myArgs)
    return myArgs


def adjustFirstDate(firstdate):
    '''
    Adjust the first date to the standard dates, so that it precedes
    the standard date  (e.g., s1 fd s2 s3) -> fd=s2
    Parameters
    ----------
    firstdate : datetime
        Nominal first date.
    Returns
    -------
    myDate : date time
        Adjusted fir date..

    '''
    myDates = mosf.standardDates()
    myDates.reverse()
    for myDate in myDates:
        if myDate['date1'] <= firstdate:
            return myDate['date1']


def getProductList(myArgs):
    '''
    Generate list of products (date1, interval) to process
    Parameters
    ----------
    myArgs : dict
        Arguments that control mosaics.

    Returns
    -------
    myProducts: .
        list of products to process
    '''
    firstDate = None
    if myArgs['noAdjustDate']:
        firstDate = myArgs['firstDate']
    sDateRanges = mosf.standardDates(nDays=myArgs['nDays'],
                                     firstDate=firstDate)
    myDates = []
    print(f'Date Range = {myArgs["firstDate"].strftime("%Y-%m-%d")}'
          f' {myArgs["lastDate"].strftime("%Y-%m-%d")}')
    for dateRange in sDateRanges:
        # print(dateRange['date1'], myArgs['firstDate'],myArgs['lastDate'])
        if dateRange['date1'] < myArgs['firstDate'] or \
                dateRange['date2'] > myArgs['lastDate']:
            continue  # Out of range so skip
        myDates.append(dateRange)
    return myDates


def makeImageMosaic(date1, date2, myArgs):
    '''
    Spawn the mosaic command.
    Parameters
    ----------
    date1, date2 : datetime
        First and second date.
    calibrate : bool
        To calibrate or not calibrate.
    Returns
    -------
    None.

    '''
    nDays = (date2-date1).days + 1
    command = f'setupimagemosaic.py --firstdate {date1.strftime("%Y-%m-%d")}' \
        f' --nDays {nDays} --template {myArgs["templateFile"]} '
    for key in ['calibrate', 'noReprocess', 'noMosaic', 'noClean',
                'setupOnly']:
        if myArgs[key]:
            command += f' --{key}'
    command += f' --prefix {myArgs["prefix"]}'
    # print(f'\n{command}')
    # u.myerror('debug')
    # make add this mosaic to thread list
    call(command, shell=True)  # , executable='/bin/csh')


def makeProdName(prefix, date1, date2):
    ''' Make product name  -- cross ref against same name prog in
    setupimagemosaic.py'''
    prodName = f'{prefix}_{date1.strftime("%d%b%y")}_' \
        f'{date2.strftime("%d%b%y")}_*_RESm_v{currentVersion:02d}.{subVersion}'
    return prodName


def checkProduct(date1, date2, myArgs, posting):
    '''
    Checks if a product has already been produced. Verifies that tifs exist
    Parameters
    ----------
    date1 : datetime
        First date.
    date2 : datetime
        Second date.
    myArgs : dict
        parameters defining mosaic creation.

    Returns
    -------
    exists bool,
        return True if directory and image products exist.
    '''
    calString = 'uncalibrated'
    if myArgs['calibrate']:
        calString = 'calibrated'
        #
    dirName = f'{myArgs["prefix"]}_{calString}.' \
        f'{date1.strftime("%Y-%m-%d")}.{date2.strftime("%Y-%m-%d")}'
    prodName = makeProdName(myArgs["prefix"], date1, date2)
    if not os.path.exists(dirName):
        return False  # No product directory so return false
    # check image products exist
    suffixes = ['image']
    if myArgs['calibrate']:
        suffixes = ['gamma0', 'sigma0']
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
        list of start and end dates {'date1': datetime, 'date2': datetime}.
    Returns
    -------
    threads; list of threads.
    '''
    threads = []
    for myProd in myProds:
        prodExists = checkProduct(myProd['date1'], myProd['date2'], myArgs,
                                  posting)
        # Put product in queue if doesn't exist or if reset mode
        if not prodExists or myArgs['reset'] or myArgs['noReprocess']:
            threads.append(threading.Thread(target=makeImageMosaic,
                                            args=[myProd['date1'],
                                                  myProd['date2'],
                                                  myArgs]))
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
    print(myProds)
    # u.myerror('stop')
    #
    threads = setupProducts(myProds, myArgs, posting)
    maxThreads = 1
    u.runMyThreads(threads, maxThreads, 'geomosaic')
    print(myProds)


if __name__ == "__main__":
    main()
