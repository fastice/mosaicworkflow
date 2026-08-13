#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 09:03:25 2023

@author: ian
"""
import argparse
import utilities as u
import glob
import sarfunc as s
import os
from datetime import datetime, timedelta
from operator import itemgetter
import json
#geodatNames  = {'S1': 'geodat10x2.in', 'TSX': 'geodat9x9.in'}

#
def getSensor(sensor):
    ''' autodetect the sensor from the directory name and return '''
    # ------ modify this for other sensors
    if sensor == 'default':
        if 'TSX' in os.getcwd():
            sensor = 'TSX'
        elif 'CSK' in os.getcwd():
            sensor = 'CSK'
        elif 'Sentinel' in os.getcwd():
            sensor = 'S1'
        elif 'NISAR' in os.getcwd():
            sensor = 'NISARTest'
        else:
            print('\n**** Could not determine sensor ****\n')
            exit()
    #
    try:
        sensorDef = s.sensorDefinitions(sensor)
    except Exception:
        u.myerror(f'Invalid sensor {sensor}')
    #
    return sensorDef


def parseFrames(frame):
    if '*' in frame:
        # all frames
        return 0, 1000000
    try:
        frames = frame.split(',')
        if len(frames) == 1:
            return int(frames[0]), int(frames[0])
        else:
            return int(frames[0]), int(frames[1])
    except Exception:
        u.myerror(f'Could not parse frames: {frame}. Must be NNN or MMM,NNN')


def parseDate(date, dt):
    try:
        return datetime.strptime(date, "%Y-%m-%d") + dt
    except Exception:
        u.myerror('Could not parse date: {date}')


def grepdateArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\033[1mFind status of all pairs for a year\033[0m',
        epilog='Part of the mosaicworkflow package.')
    parser.add_argument('year', type=int, nargs=1, help='year')
    parser.add_argument('--sensor', type=str, default='default',
                        help='sensor name [S1,TSX,CSK] - only needed if '
                        'sensor name is not in path (Sentinel,TSX,CSK)')
    parser.add_argument('--firstdate', type=str, default=None,
                        help='first date in range to override year YYYY-MM-DD')
    parser.add_argument('--lastdate', type=str, default=None,
                        help='last date in range to override year YYYY-MM-DD')
    parser.add_argument('--frame', type=str, default='*',
                        help='frame (for single) frame1, frame2 (for range)'
                        ' [*, all frames]')
    args = parser.parse_args()
    frame1, frame2 = parseFrames(args.frame)
    year = args.year[0]
    sensor = getSensor(args.sensor)
    # process first date
    if args.firstdate is None:
        date1 = datetime(year, 1, 1)
    else:
        date1 = parseDate(args.firstdate, timedelta(seconds=0))
    # sprocess second date
    if args.lastdate is None:
        date2 = datetime(year, 12, 31, 23, 59, 59)
    else:
        date2 = parseDate(args.lastdate, timedelta(seconds=86399.9))
    #
    if year < 2008 or year > 2040:
        u.myerror('invalid year {0:d}'.format(year))
    return date1, date2, sensor, frame1, frame2


def getDate(product, sensor):
    '''
    Get Product date from geodat
    '''
    geoFile = f'{product}/{sensor.geodatName()}'
    if 'geojson' in geoFile and os.path.exists(geoFile):
        with open(geoFile, "r") as f:
            geojson = json.load(f)
            date = f'{geojson["properties"]["Date"]} ' \
                f'{geojson["properties"]["NominalTime"]}'
            return datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f')
    # Example: get "date" from the first feature
    elif os.path.exists(geoFile):
        with open(geoFile, 'r') as fp:
            _ = fp.readline()
            dateline = fp.readline().split(':')[-1].strip().replace(' ', '')
            timeline = fp.readline().split(':')[-1].strip().replace(' ', ':')
            date = f'{dateline} {timeline}'
            return datetime.strptime(date, '%d%b%Y %H:%M:%S.%f')
    else:
        # Return clearly out of range date
        return datetime(1900, 1, 1)


def getProducts(date1, date2, frame1, frame2, sensor):
    '''
    Get Product dirs
    Returns
    -------
    filtered list of products.

    '''
    productDirs = []
    for pattern in ['????_*', '?????_*', '???_*']:
        productDirs += glob.glob(pattern)
    products = []
    for product in productDirs:
        productDict = {}
        productDict['dir'] = product
        # Parse orbit/frame
        try:
            productDict['orbit'], productDict['frame'] = \
                        [int(x) for x in product.split('_')]
            # Filter by Frame
            if productDict['frame'] < frame1 or productDict['frame'] > frame2:
                # not in frame range so contiune
                continue
        except Exception:
            # It wasn't a directory NNNN_MMM so continue
            continue
        #
        # Filte by date
        date = getDate(product, sensor)
        if date >= date1 and \
                date <= (date2 + timedelta(days=sensor.SAR['maxDays'])):
            productDict['date'] = date
            products.append(productDict)
    products = sorted(products, key=itemgetter('date'))
    return products


def findPair(products, i, sensor):
    firstImage = products[i]
    for j in range(i+1, min(i+10, len(products))):
        if products[j]['frame'] == firstImage['frame']:
            return products[j]
    return None


def slcStatus(firstImage):
    #
    if os.path.exists(f'{firstImage["dir"]}/{firstImage["dir"]}.slc'):
        return '+'
    return '-'


def processingStatus(firstImage):
    if os.path.exists(f'{firstImage["dir"]}/azimuth.offsets'):
        thumbs = '-'
        if os.path.exists(f'{firstImage["dir"]}/velocity'):
            thumbs = 'v'
        fast = '.'
        if os.path.exists(f'{firstImage["dir"]}/fast/azimuth.offsets.fast'):
            fast = 'f'
        return f'x{fast}{thumbs}'
    elif os.path.exists(f'{firstImage["dir"]}/runboth'):
        return '.--'
    else:
        return 'o--'


def pairInterval(firstImage, secondImage):
    ''' Compute temporal separation and round to nearest day '''
    if secondImage is None:
        return 0
    dt = (secondImage['date'] - firstImage['date'] + timedelta(hours=12)).days
    return dt


def processPair(firstImage, secondImage):
    '''
    Give a pair of image directories compute summary line
    '''
    slc = slcStatus(firstImage)
    pStatus = processingStatus(firstImage)
    dt = pairInterval(firstImage, secondImage)
    orbFrame = f'{firstImage["orbit"]}_{firstImage["frame"]}'
    print(f'{pStatus} {dt:3}d {orbFrame:>10} '
          f' {firstImage["date"].strftime("%-d-%^b-%Y"):>11}  {slc}  '
          f'{firstImage["date"].strftime("%-j"):>3}')


def outputSummary(products, date1, date2, sensor):
    '''
    Loop through products to:
    1) find the next image to form a pair
    2) Compute and print the pair summary line
    '''
    for i in range(0, len(products)):
        firstImage = products[i]
        if firstImage['date'] >= date1 and firstImage['date'] <= date2:
            secondImage = findPair(products, i, sensor)
            processPair(firstImage, secondImage)


def printHeader(products, sensor):
    '''
    Print traditional grepdate/greptops header line
    Returns
    -------
    None.
    '''
    direction = None
    # First check for ascending/descending in track dir
    if os.path.exists('ascending'):
        direction = '****  Ascending\n'
    elif os.path.exists('descending'):
        direction = '***  Descending\n'
    # Fall and check products until result found
    else:
        # Loop through products until a direction is found (should be first)
        for product in products:
            if os.path.exists(f'{product["dir"]}/{sensor.geodatName()}'):
                with open(f'{product["dir"]}/{sensor.geodatName()}') as fp:
                    for line in fp:
                        if 'ing Pass' in line:
                            direction = f'***  {line.split()[1].strip()}\n'
                            break
            if direction is not None:
                break
    if direction is None:
        u.myerror('Could not finding ascending/descending direction.'
                  ' In track dir?')
    # print the header
    print(f'\n{direction}')


def main():
    #
    # Get command line args
    date1, date2, sensor, frame1, frame2 = grepdateArgs()
    #
    #
    products = getProducts(date1, date2, frame1, frame2, sensor)
    #
    # Print 4 line header
    printHeader(products, sensor)
    #
    # Print line by line results
    outputSummary(products, date1, date2, sensor)


if __name__ == '__main__':
    main()
