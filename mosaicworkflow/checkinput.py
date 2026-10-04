#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
checkinput - check every file a mosaic3d master input file references.

    checkinput.py track-all/inputFile [--phaseTieThresh N] [--phaseCoverage F]
                  [--sigmaTieThresh S]

See Documents/checkinput.md. Created on Mon Nov 22 10:26:23 2021, @author: ian
"""
import argparse
import utilities as u
import os
from datetime import datetime
import numpy as np
from osgeo import gdal

def checkinputArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(description='\n\n\033[1mCheck all files '
                                     'exist for mosaic 3d input file\033[0m\n')
    parser.add_argument('inputFile', type=str, nargs=1,
                        help='Inputfile to check')
    parser.add_argument('--phaseTieThresh', type=int, default=None,
                        help='Add exclude to phase directories with fewer '
                        'than phaseTieThresh tiepoints')
    parser.add_argument('--phaseCoverage', type=float, default=None,
                        help='Exclude if less than phaseCoverage and less '
                        'than 10 tiepoints')
    parser.add_argument('--sigmaTieThresh', type=float, default=None,
                        help='Issue warning for phase sigma > thresh')
    args = parser.parse_args()
    tieParams = {'phaseTieThresh': args.phaseTieThresh,
                 'sigmaTieThresh': args.sigmaTieThresh,
                 'phaseCoverage': args.phaseCoverage,
                 'azSigma': 1.0,
                 'rgSigma': 1.0,
                 'svConst': True}
    for tie in tieParams:
        print(f'{tie}:  {tieParams[tie]}')
    return args.inputFile[0], tieParams


def isYaml(paramFile):
    return paramFile.endswith('.yaml')


def readYamlParams(paramFile):
    ''' (sigma, nUsed) from a yaml parameter file (az.est.yaml,
    rBaseline.deltabp.yaml, baseline.*.yaml): the first `sigma:` and
    `nTiepointsUsed:`/`nTiepoints:` keys. rBaseline yaml keeps a legacy `;`
    header, so parse lines rather than load it as YAML. sigma <= 0 means no
    solution. '''
    sigma, nUsed = None, None
    with open(paramFile, 'r') as fp:
        for line in fp:
            key, _, value = line.strip().partition(':')
            if not value.strip():
                continue
            if key == 'sigma' and sigma is None:
                sigma = float(value.split()[0])
            elif key in ['nTiepointsUsed', 'nTiepoints'] and nUsed is None:
                nUsed = int(value.split()[0])
    return sigma, nUsed


def checkYamlParams(paramFile, count, sigmaThresh):
    ''' yaml version of the checkAzEst/checkRbaseline tests. Returns the
    number of tiepoints used (0 for no solution). '''
    if not os.path.exists(paramFile):
        print(f'\033[1;31m{count} Missing {paramFile}\033[0m')
        return 0
    sigma, nUsed = readYamlParams(paramFile)
    if sigma is None or sigma <= 0 or not nUsed:
        print(f'\033[1;31m{count} no solution for {paramFile}\033[0m')
        return 0
    if sigmaThresh is not None and sigma > sigmaThresh:
        print(f'{count} warning sigma = {sigma} for {paramFile}')
    if nUsed <= 10:
        print(f'{count} only {nUsed} tiepoints for {paramFile}')
    return nUsed


def readInput(inputFile, svConst):
    ''' Read input and store info as list of dictionaries '''
    keys = ['phaseFile', 'geodat', 'baseline', 'nDays', 'weight', 'azimuthOff',
            'azimuthParams', 'rangeOff', 'rBaseline', 'flag']
    types = [str, str, str, int, float, str, str, str, str, int]
    products = []
    with open(inputFile, 'r') as fp:
        for line in fp:
            pieces = line.split()
            if ';' not in line and len(pieces) >= 9:
                try:
                    pieces = [myType(x) for x, myType in zip(pieces, types)]
                    myDict = dict(zip(keys, pieces))
                    track = myDict['geodat'].split('/')[-4]
                    myDict['track'] = track
                    if 'TSX' in myDict['geodat']:
                        myDict['sensor'] = 'TSX'
                    elif 'Sentinel' in myDict['geodat']:
                        myDict['sensor'] = 'Sentinel1'
                    else:
                        myDict['sensor'] = 'unkown'
                    # yaml params already name the const/deltabp solution
                    if svConst:
                        if 'None' not in myDict['azimuthParams'] and \
                                not isYaml(myDict['azimuthParams']):
                            myDict['azimuthParams'] += '.const'
                        if 'None' not in myDict['rBaseline'] and \
                                not isYaml(myDict['rBaseline']):
                            myDict['rBaseline'] += '.deltabp'
                except Exception:
                    u.myerror(f'Cannot parse line: {line}')
                #
                products.append(myDict)
    return products


def parseDatFile(datFile):
    ''' parse data information'''
    if 'None' in datFile:
        return None, None
    if not os.path.exists(datFile):
        print(f'Missing dat file {datFile}')
    fp = open(datFile, 'r')
    firstLine = fp.readline().split()
    nr, na = [int(x) for x in firstLine[2:4]]
    return nr, na


def checkFile(myFile, nr, na, nGood=False):
    if 'None' in myFile:
        return
    # offsets moved to tif+vrt keep the bare name in the input file
    if not os.path.exists(myFile) and os.path.exists(f'{myFile}.vrt'):
        myFile = f'{myFile}.vrt'
    if not os.path.exists(myFile):
        print(f'Missing file {myFile}')
        return None
    if myFile.endswith('.vrt'):
        # VRT phase (S1 phase-only products): read through GDAL, no size test
        if nGood:
            x = gdal.Open(myFile).ReadAsArray()
            return float(np.sum(x > -2e8)) / float(x.size)
        return None
    size = os.path.getsize(myFile)
    if size < nr * na * 4:
        print(f'Offset file: {myFile} missing data')
    if nGood:
        x  = u.readImage(myFile, nr, na,'>f4')
        nValid = np.sum(x > -2e8)
        # print(nValid, nr*na)
        return float(nValid)/float(nr * na)
    else:
        return None

def checkAzEst(azEst, count, azSigma):
    ''' Check valid az.est'''
    noSolution = True
    sigmaFound = False
    if 'None' in azEst:
        return
    # yaml given, or the legacy name now stored only as <name>.yaml
    if not isYaml(azEst) and not os.path.exists(azEst) and \
            os.path.exists(f'{azEst}.yaml'):
        azEst = f'{azEst}.yaml'
    if isYaml(azEst):
        checkYamlParams(azEst, count, azSigma)
        return
    if not os.path.exists(azEst):
        if 'TSX' in azEst:
            azEst = azEst.replace('.const', '')
        else:
           print(f'\033Missing {azEst}\033[0m')
           return
    with open(azEst, 'r') as fp:
        for line in fp:
            if 'sigma' in line:
                sigmaFound = True
                sigma = float(line.split('=')[1].strip().split()[0])
                if sigma == 0.0 or sigma > azSigma:
                    print(f'{count} warning sigma = {sigma} for {azEst}')
            elif 'Ntiepoints' in line:
                nUsed, nGiven = [int(x) for x in
                                 line.split('=')[1].strip().split('/')]
                if nUsed <= 10:
                    print(f'{count} only {nUsed} tiepoints for {azEst}')
            elif ';' not in line:
                try:
                    params = [float(x) for x in line.split()]
                    if len(params) >= 4 and len(params) <= 6:
                        noSolution = False
                except Exception:
                    noSolution = True
    if noSolution or not sigmaFound:
        print(f'\033[1;31m{count} no solution for {azEst}\033[0m')


def checkRbaseline(rBaseline, count, sigmaThresh=20):
    ''' Check valid rBaseline '''
    noSolution = True
    sigmaFound = False
    if 'None' in rBaseline:
        return
    if not isYaml(rBaseline) and not os.path.exists(rBaseline) and \
            os.path.exists(f'{rBaseline}.yaml'):
        rBaseline = f'{rBaseline}.yaml'
    if isYaml(rBaseline):
        return checkYamlParams(rBaseline, count, sigmaThresh)
    if not os.path.exists(rBaseline):
        if 'TSX' in rBaseline:
            rBaseline = rBaseline.replace('.deltabp', '')
        else:
           print(f'\033Missing {rBaseline}\033[0m')
           return 0
    # print(rBaseline)
    with open(rBaseline, 'r') as fp:
        for line in fp:
            if 'sigma' in line:
                sigmaFound = True
                sigma = float(line.split('=')[1].strip().split()[0])
                if sigma == 0.0:
                    print(f'{count} warning sigma = {sigma} for {rBaseline}')
                if sigmaThresh is not None:
                    if sigma > sigmaThresh:
                        print(f'{count} warning sigma = {sigma} for '
                              f'{rBaseline}')
            elif 'Ntiepoints' in line:
                nUsed, nGiven = [int(x) for x in
                                 line.split('=')[1].strip().split('/')]
                if nUsed <= 10:
                    print(f'{count} only {nUsed} tiepoints for {rBaseline}')
            elif ';' not in line and '&' not in line:
                try:
                    params = [float(x) for x in line.split()]
                    # print(params)
                    if len(params) > 4 and len(params) <= 7:
                        noSolution = False
                except Exception:
                    # print(line)
                    noSolution = True
    if noSolution or not sigmaFound:
        print(f'\033[1;31m{count} no solution for {rBaseline}\033[0m')
        nUsed = 0
        # u.myerror('stop')
    return nUsed

def addExcludeIfNeeded(phaseThresh, nUsed, phaseFile, msg=None):
    ''' Add Excluded to phase directories for too few tiepoints'''
    if phaseThresh is None:
        return
    if nUsed <= phaseThresh:
        # legacy phase is <frame>/motion/phase.reflat.uw; phase-only
        # products keep it in the frame dir (<frame>/phase.uw.interp.vrt)
        if '/motion' in phaseFile:
            phaseDir = phaseFile.split('/motion')[0]
        else:
            phaseDir = os.path.dirname(phaseFile)
        with open(f'{phaseDir}/Exclude.pending', 'w') as fp:
            if msg is None:
                msg = \
                    f'{datetime.now()}: Excluding with {nUsed} phase tiepoints'
            print(msg, file=fp)
            print(msg)

def parseOffsetSize(azOff):
    if 'None' in azOff:
        return None, None
    vrtFile = f'{azOff}.vrt'
    datFile = f'{azOff}.dat'
    myOffsets = u.offsets(verbose=False, myPath=os.path.dirname(azOff))
    if os.path.exists(vrtFile):
        myOffsets.readOffsetsDat(vrtFile=vrtFile)
    elif os.path.exists(datFile):
        myOffsets.readOffsetsDat(datFile=os.path.basename(datFile))
    else:
        return None, None
    #
    return myOffsets.nr, myOffsets.na


def checkProduct(product, count, tieParams):
    ''' Validate files for a product '''
    nr, na = parseOffsetSize(product["azimuthOff"])
    checkFile(product['azimuthOff'], nr, na)
    checkFile(product['rangeOff'], nr, na)
    #
    geo = u.geodatrxa(file=product['geodat'])
    #
    checkAzEst(product['azimuthParams'], count, tieParams['azSigma'])
    nUsed = checkRbaseline(product['rBaseline'], count,
                           sigmaThresh=tieParams['rgSigma'])
    #
    if 'nophase' not in product['phaseFile']:
        nUsed = checkRbaseline(product['baseline'], count,
                               sigmaThresh=tieParams['sigmaTieThresh'])
        #
        # Check low coverage cases
        if nUsed <= 10 and tieParams['phaseCoverage'] is not None:
            # Compute fraction covered.
            phaseCoverage = checkFile(product['phaseFile'], geo.nr, geo.na,
                                      nGood=True)
            if phaseCoverage <  tieParams['phaseCoverage']:
                msg = f'Excluding {product["phaseFile"]} for no solution and' \
                    f' low coverage {phaseCoverage:0.3f}'
                addExcludeIfNeeded(10, nUsed, product['phaseFile'], msg=msg)
           # print(phaseCoverage)
        addExcludeIfNeeded(tieParams['phaseTieThresh'], nUsed,
                           product['phaseFile'])


    return geo.datetime


def updateDates(myDates, myDate, track):
    ''' Update dates to get the most recent date for the file'''
    if track in myDates:
        if myDates[track] > myDate:
            return
    # First or new date
    myDates[track] = myDate


def printDates(myDates):
    for myDate in myDates:
        print(f'{myDate} {myDates[myDate].strftime("%Y--%m-%d")}')


def main():
    '''
    main for checkinput
    '''
    # Parse Args
    inputFile, tieParams = checkinputArgs()
    print(f'Checking {inputFile}')
    # Get input data
    products = readInput(inputFile, tieParams['svConst'])
    print(products[-1])
    count = 1
    s1Dates = {}
    tsxDates = {}
    for product in products:
        myDate = checkProduct(product, count, tieParams)
        if product['sensor'] == 'TSX':
            updateDates(tsxDates, myDate, product['track'])
        elif product['sensor'] == 'Sentinel1':
            updateDates(s1Dates, myDate, product['track'])
        count += 1
    print('\n*** TSX MOST RECENT')
    printDates(tsxDates)
    print('\n*** Sentinel 1 MOST RECENT')
    printDates(s1Dates)

if __name__ == '__main__':
    main()
