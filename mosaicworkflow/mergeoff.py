#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 14:53:43 2018

@author: ian
"""

import utilities as u
import sarfunc as s
import argparse
import copy
import numpy as np
import scipy.io as io
import os
import rasterio

myLog = u.logger(fileRoot='mergeoff')


def mergeOffArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\033[1mMerge fast offsets with regular offsets \033[0m',
        epilog='Usually called as part of cleanoffmerge.py. '
               'Part of the mosaicworkflow package.')
    parser.add_argument('sensor', type=str, default='S1',
                        help='Sensor (CSK, S1, TSX)')
    parser.add_argument('--tiff', action='store_true', default=False,
                        help='Write the fast/merged offsets as GeoTIFF + '
                        'tiff-backed VRT instead of raw binary')
    #
    args = parser.parse_args()
    # log arguments
    #
    sarDef = s.sensorDefinitions(args.sensor)
    sensorInfo = sarDef.SAR
    #
    return sensorInfo, args.tiff


def readRegularOffsets(offsetFile, vrtFile=None, tiff=False):
    if vrtFile is not None:
        if not os.path.exists(vrtFile):
            vrtFile = None  # Older products have no vrt
    # Now read vrt or data files
    if vrtFile is not None:
        offR = u.offsets(vrtFile=vrtFile, tiff=tiff)
        offR.readOffsets(vrtFile=vrtFile)
    else:
        offR = u.offsets(offsetFile, datFile=f'{offsetFile}.dat', tiff=tiff)
        offR.readOffsets(datFile=f'{offsetFile}.dat')
    offR.readSigma()
    return offR


def readFastOffsets(fastFile, fastPath, vrtFile=None):
    # open offsets
    if vrtFile is not None:
        if not os.path.exists(vrtFile):
            vrtFile = None  # Older products have no vrt
    # now read
    if vrtFile is not None:
        offF = u.offsets(vrtFile=vrtFile, verbose=False)
        offF.readOffsets()
        offF.readSigma()
    else:
        offF = u.offsets(fastFile, datFile=f'{fastFile}.dat', myPath=fastPath)
        # read data
        offF.readOffsets()
        offF.readSigma(sigmaAFile=offF.sigmaAFile.replace('.noclean', ''),
                       sigmaRFile=offF.sigmaRFile.replace('.noclean', ''))
    #
    return offF


def cleanFastOffsets(offsetsFile, offReg, vrtFile=None, tiff=False):
    ''' This routine produces a file that applies the idl handcleaning, which
    is common for TSX but rare for S1. This input is used to produce
    velocity_nocull results. '''
    if vrtFile is not None:
        if not os.path.exists(vrtFile):
            vrtFile = None  # Older products have no vrt
    if vrtFile is not None:
        fastOffOrig = u.offsets(vrtFile=vrtFile, tiff=tiff)
    else:
        fastOffOrig = u.offsets(offsetsFile, datFile=offsetsFile+'.dat',
                                myPath='./fast', tiff=tiff)
    fastOffOrig.readOffsets(vrtFile=vrtFile)
    # readOffsets only loads the offsets bands; load sigma too so writeOffsets
    # emits the sigma outputs (required in tiff mode, where the sigma tiffs do
    # not already exist as they would as raw files).
    if vrtFile is not None:
        fastOffOrig.readSigma(vrtFile=vrtFile)
    if os.path.exists('fast/badoffsets.idl'):
        badoff = io.readsav('fast/badoffsets.idl')
        print(badoff)
        fastOffOrig.removeList(badoff.idel)
    fastOffOrig.offsetFileNames('azimuth.offsets.fast', myPath='./fast',
                                updateDatFileName=True)
    # Set sigma output names (the vrt constructor skips this and readSigma from
    # a vrt does not set them) so writeOffsets emits the fast sigma outputs.
    fastOffOrig.sigmaFileName()
    #
    # Copy the geo1/2 files so that velocity_nocull can use stateV solution
    fastOffOrig.geo1 = f'../{offReg.geo1}'
    fastOffOrig.geo2 = f'../{offReg.geo2}'
    if vrtFile is not None:
        fastOffOrig.writeOffsets(noDatFiles=True)
        fastOffOrig.writeOffsetVrt('./fast/offsets.fast.vrt',
                                   ['range.offsets.fast',
                                    'range.offsets.fast.sr',
                                    'azimuth.offsets.fast',
                                    'azimuth.offsets.fast.sa'],
                                   ['RangeOffsets', 'RangeSigma',
                                    'AzimuthOffsets', 'AzimuthSigma'],
                                   byteOrder=None)
        fastOffOrig.writeOffsetVrt('./fast/range.offsets.fast.vrt',
                                   ['range.offsets.fast',
                                    'range.offsets.fast.sr'],
                                   ['RangeOffsets', 'RangeSigma'],
                                   byteOrder=None)
        fastOffOrig.writeOffsetVrt('./fast/azimuth.offsets.fast.vrt',
                                   ['azimuth.offsets.fast',
                                    'azimuth.offsets.fast.sa'],
                                   ['AzimuthOffsets', 'AzimuthSigma'],
                                   byteOrder=None)
    else:
        fastOffOrig.writeOffsets()


def doTheMerge(offR, offF, sensorInfo):
    # copy speckle offsets
    offM = copy.deepcopy(offR)
    print('---', np.sum(offM.areValid()))
    # average where both exist
    iSame = np.logical_and(offR.areValid(), offF.areValid())
    wF = sensorInfo['fastW']
    wR = sensorInfo['regW']
    #
    offM.rgOff[iSame] = wR*offR.rgOff[iSame] + wF*offF.rgOff[iSame]
    offM.azOff
    offM.sigmaR[iSame] = np.sqrt(wR * offR.sigmaR[iSame]**2 +
                                 wF * offF.sigmaR[iSame]**2)
    offM.sigmaA[iSame] = np.sqrt(wR * offR.sigmaA[iSame]**2 +
                                 wF * offF.sigmaA[iSame]**2)
    print(np.min(offM.sigmaA[iSame]), np.max(offM.sigmaR[iSame]))
    print(np.min(offR.sigmaA[iSame]), np.max(offR.sigmaR[iSame]))
    print(np.min(offF.sigmaA[iSame]), np.max(offF.sigmaR[iSame]))
    # Make sure the fast offsets files have the geodat information so that
    # qa can use svfit for velocity_nocull

    # where there are only smooth, add those in
    iFast = np.logical_and(offR.notValid(), offF.areValid())
    offM.rgOff[iFast] = offF.rgOff[iFast]
    offM.azOff[iFast] = offF.azOff[iFast]
    offM.sigmaR[iFast] = offF.sigmaR[iFast]
    offM.sigmaA[iFast] = offF.sigmaA[iFast]
    print('+++', np.sum(offM.areValid()))
    return offM


def applySECorrection(off):
    SEFile = 'offsets.SECorrection'
    # raw <SEFile> in the legacy path, <SEFile>.tif with SETideOffsets --tiff;
    # the vrt is named the same either way and is what is actually read.
    if os.path.exists(f'{SEFile}.vrt') and \
            (os.path.exists(SEFile) or os.path.exists(f'{SEFile}.tif')):
        print('****** Applying SE Correction******')
        correction = np.squeeze(rasterio.open(f'{SEFile}.vrt').read())
        off.rgOff[off.areValid()] -= correction[off.areValid()]


def main():
    ''' Merges regular tracked offsets with fast tracked offsets. '''
    #
    # get command line args
    sensorInfo, tiff = mergeOffArgs()
    # filter offsets if needed to create azimuth.offsets.fast
    #
    # get offsets
    offReg = readRegularOffsets('azimuth.offsets.slow',
                                vrtFile='offsets.slow.vrt', tiff=tiff)
    #
    cleanFastOffsets('azimuth.offsets.noclean.fast', offReg,
                     vrtFile='./fast/offsets.noclean.fast.vrt', tiff=tiff)
    offFast = readFastOffsets('azimuth.offsets.fast', './fast',
                              vrtFile='./fast/offsets.fast.vrt')
    # Apply to fast, azimuth offsets already corrected in cleanoff.py
    applySECorrection(offFast)
    #
    offReg.sigmaFileName()
    #
    offMerge = doTheMerge(offReg, offFast, sensorInfo)
    #
    offMerge.offsetFileNames('azimuth.offsets', myPath='.',
                             updateDatFileName=True)
    offMerge.sigmaFileName()
    #
    offMerge.writeOffsets(noDatFiles=True)


if __name__ == '__main__':
    main()
