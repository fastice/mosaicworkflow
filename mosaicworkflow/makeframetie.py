#!/usr/bin/env python3
import utilities as u
import sys
import os
import argparse
from subprocess import call
import numpy as np
import yaml
import glob
#
# This program reads a tiefile and strips out the pieces to make a tiefile
# that can be used generate a new tiefile
# divided by frame (frame range comes from ../velocityStats/X-Y.
# Next it runs the tie_script and then uses the results with vel_thumbs
#
#  After completing, all the vel_thumbs should be in place with a common
# region for each frame range
#
#  Before running sure the following exist
#  ../velocityStats/X-Y
#  tiepoints/vel_thumb_header_XdashY (select the common region in these)
# In running it will create vel_thumb_plan_allYeardashXdashY and
# tie_planYear.XdashY
#


def getSensorInfo(tiePlanFile):
    cwd = os.getcwd()
    multiTrack = False
    projectDir = os.path.dirname(os.path.dirname(cwd))
    projectYaml = os.path.join(projectDir, 'project.yaml')
    sensorYaml = os.path.join(projectDir, 'sensor.yaml')
    if os.path.isfile(projectYaml):
        yamlFile = projectYaml
    elif os.path.isfile(sensorYaml):
        u.mywarning('sensor.yaml found but project.yaml expected — rename to project.yaml')
        yamlFile = sensorYaml
    else:
        yamlFile = None
    if yamlFile is not None:
        with open(yamlFile, 'r') as f:
            sensorData = yaml.safe_load(f)
        yamlSensor = sensorData.get('sensor', None)
        framePattern = sensorData.get('framePattern', '*')
        print(f'\033[1;31mDEBUG: read framePattern={framePattern!r} from {yamlFile}\033[0m')
        if yamlSensor is not None:
            if yamlSensor in ['TSX', 'CSK']:
                trackdirs = sorted(d.rstrip('/') for d in glob.glob('../track-*/'))
                if len(trackdirs) > 1:
                    trackUsed = tiePlanFile.split('-')[-1]
                    print(trackUsed)
                    for trackdir in trackdirs:
                        track = trackdir.split('-')[-1]
                        if track == trackUsed:
                            print('multiTrack')
                            multiTrack = True
                            break
                    if not multiTrack:
                        u.myerror('tieplan20XX-track does not match one of '
                                  'multiple tracks')
                else:
                    trackdir = trackdirs[0]
                    track = trackdir.split('-')[-1]
                    trackUsed = track
                dirs = [x.rstrip('/').split('/')[-1] for x in
                        sorted(glob.glob(f'../track-{trackUsed}/velocityStats/*-*/'))]
                print(f'Track used {trackUsed}')
                print(dirs)
            else:
                dirs = sorted(d.rstrip('/') for d in glob.glob('../velocityStats/*-*/'))
                if len(dirs) == 0:
                    u.myerror(f'no ../velocityStats/*-* {dirs}')
                trackdir = cwd.split('/')[-2]
                track = trackdir.split('-')[-1]
                if trackdir.find('track') < 0:
                    u.myerror('error: not subdir of track-x')
            return yamlSensor, dirs, trackdir, track, multiTrack, framePattern
    if 'TSX' in cwd or 'CSK' in cwd:
        if 'TSX' in cwd:
            sensor = 'TSX'
        else:
            sensor = 'CSK'
        # all frames
        #dirs = ['/1-1000']
        # track dir
        trackdirs = sorted(d.rstrip('/') for d in glob.glob('../track-*/'))
        if len(trackdirs) > 1:
            trackUsed = tiePlanFile.split('-')[-1]
            print(trackUsed)
            for trackdir in trackdirs:
                track = trackdir.split('-')[-1]
                if track == trackUsed:
                    print('multiTrack')
                    multiTrack = True
                    break
            if len(trackdirs) > 1 and not multiTrack:
                u.myerror('tieplan20XX-track does not match one of '
                          'multiple tracks')
        else:
            trackdir = trackdirs[0]
            track = trackdir.split('-')[-1]
            trackUsed = track
        dirs = [x.rstrip('/').split('/')[-1] for x in
                sorted(glob.glob(f'../track-{trackUsed}/velocityStats/*-*/'))]
        print(f'Track used {trackUsed}')
        print(dirs)
        # print('t ', trackdir)
    elif 'Sentinel' in cwd or 'ALOS2' in cwd:
        if 'Sentinel' in cwd:
            sensor = 'Sentinel1'
        else:
            sensor = 'ALOS2'
        # find frame ranges
        dirs = sorted(d.rstrip('/') for d in glob.glob('../velocityStats/*-*/'))
        if len(dirs) == 0:
            u.myerror(f'no ../velocityStats/*-* {dirs}')
        #
        trackdir = os.getcwd().split('/')[-2]
        track = trackdir.split('-')[-1]
        if trackdir.find('track') < 0:
            u.myerror('error: not subdir of track-x')
    else:
        u.myerror('invalid sensor make sure directory has TSX, or '
                  'Sentinel in path')
    return sensor, dirs, trackdir, track, multiTrack, '*'


def main():
    parser = argparse.ArgumentParser(
        description='Use tiefiles to generate vel_thumbs for velocityStats',
        epilog=(
            'Sensor and frame configuration:\n'
            '  Reads project.yaml (two directories up, i.e. the project root)\n'
            '  for "sensor" and "framePattern". Falls back to sensor.yaml with\n'
            '  a warning if project.yaml is not found.\n\n'
            'Notes for Sentinel 1 / NISAR:\n'
            '  1) Before running ensure the following exist:\n'
            '     a) ../velocityStats/X-Y\n'
            '     b) tiepoints/vel_thumb_header_XdashY '
            '(selects the common region)\n'
            '  2) Creates vel_thumb_plan_allYeardashXdashY and '
            'tie_planYear.XdashY\n'
            '  3) tie_script and vel_thumbs are called to generate '
            'velocity maps\n\n'
            'Notes for TSX:\n'
            '  1) If there are multiple tracks, the format is '
            'tie_plan20XX-TTT, otherwise -TTT is omitted'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'tiePlan',
        metavar='tie_plan',
        help='Tie plan file (e.g. tie_plan20XX or tie_plan20XX-TTT for TSX)'
    )
    parser.add_argument(
        '--overWrite', action='store_true',
        help='Pass --overWrite to vel_thumbs to rerun existing products'
    )
    parser.add_argument(
        '--keepVz', action='store_true',
        help='Pass --keepVz to vel_thumbs to retain .vz and .vz.geodat files'
    )

    args = parser.parse_args()
    tiePlanFile = args.tiePlan
    overWriteFlag = ' --overWrite' if args.overWrite else ''
    keepVzFlag = ' --keepVz' if args.keepVz else ''

    headers = ['extraties', 'DEM', 'track_root', 'base_nlooks',
               'default_nDays', 'extra_flags']
    #
    # Figure out sensor
    sensor, dirs, trackdir, track, multiTrack, framePattern = \
        getSensorInfo(tiePlanFile)
    framePrefix = framePattern.split('?')[0] if '?' in framePattern else ''
    print('sensor = ', sensor)
    print('trackdir = ', trackdir)
    #
    #
    # loop over dirs
    #
    for d in dirs:
        frameRange = d.split('/')[-1]
        f1 = int(frameRange.split('-')[0][len(framePrefix):])
        f2 = int(frameRange.split('-')[1][len(framePrefix):])
        print(frameRange, f1, f2, d)
        try:
            fin = open(tiePlanFile, 'r')
        except Exception:
            u.myerror('could not open input file ' + tiePlanFile +
                      ' specified on the command line ')
        frameRange = frameRange.replace('-', 'dash')
        # cycle through lines in field
        all = []
        # this is used to determine where to put the start(before first subdir)
        firstSub = True
        #
        # New frame so new tie file
        #
        tiefile = tiePlanFile + '.' + frameRange
        print(f'**** processing {tiefile}')
        try:
            fout = open(tiefile, 'w')
        except Exception:
            u.myerror('could not open output file ' + tiefile + ': check path  ')
        tieString = ''
        subString = ''
        for line in iter(fin):
            # set this flag to true when a line processed, so it can be
            # skipped in subsequent steps
            processed = False
            # pass through header lines
            if any(s in line for s in headers):
                print(line.format('%s'), end='', file=fout)
                if line.find('default_nDays') != -1:
                    nDays = line.split('=')[-1]
                processed = True
            # pass through blank lines (only 1 in a sequence)
            # get ndays
            if line.find('subdir') != -1:
                if firstSub:
                    print('start', file=fout)
                    firstSub = False
                if line.find('nDays') != -1:
                    nDays = int(line.split('=')[-1])
                    # print('\n'+line, end='', file=fout)
                    # save subdir string and output later if a frame is output
                    subString = '\n' + line
                    processed = True
                else:
                    # print('\n'+line.rstrip()+'nDays='+nDays,end='',file=fout)
                    # save subdir string and output later if a frame is output
                    subString = '\n' + line.rstrip() + ' nDays=' + nDays
                    processed = True
                # save tiefile string and output later if a frame is output
                tieString = '\n\ttiefile ' + frameRange + 'd' + str(nDays) + '\n'
                # print('\n\ttiefile '+frameRange+'d'+str(nDays),
                # end='\n', file=fout)
                # proc
            if not processed:
                pieces = line.split()
                if len(pieces) >= 3:
                    pieces[0].strip()
                    frame = pieces[1]
                    if len(pieces[0]) == 1 and pieces[0].find('n') != -1:
                        frame_str = pieces[2].split('_')[-1][len(framePrefix):]
                        if not frame_str:
                            continue  # frame shorter than prefix, not a virtual-frame entry
                        frame = int(frame_str)
                        # only print frame info its in the range
                        if frame >= f1 and frame <= f2:
                            if len(subString) > 0:
                                print(subString, end='', file=fout)
                                subString = ''
                            if len(tieString) > 0:
                                print(tieString, end='', file=fout)
                                s = tieString.replace('tiefile', '').strip()
                                # this saves the tie string for output with all
                                all.append(s)
                                tieString = ''
                            print('\t\t' + line.lstrip(), end='', file=fout)
                if len(pieces) == 2:
                    if (pieces[0].find('tiefile') != -1 and
                            pieces[1].find('all') != -1):
                        allRoot = pieces[1].strip()
        # print orbit data
        print(f'\n\ttiefile {allRoot}dash{frameRange}', file=fout)
        count = 0
        for ll in all:
            print(f'\t\tc {track} {ll}', file=fout)
            count += 1
        print('end', file=fout)
        fout.close()
        fin.close
        #
        #
        #
        if count > 0:
            thumbfile = f'vel_thumb_plan_{allRoot}dash{frameRange}'
            if multiTrack:
                fin = open(f'vel_thumb_header_{frameRange}-{track}', 'r')
                thumbfile += '-' + track
            else:
                fin = open(f'vel_thumb_header_{frameRange}', 'r')

            fout = open(thumbfile, 'w')
            for line in fin:
                print(line, end='', file=fout)
            print('start', file=fout)
            print(f'{trackdir} {allRoot}dash{frameRange}', file=fout)
            print('end', file=fout)
            fin.close()
            fout.close()
            #
            call('tie_script  ' + tiefile, shell=True)
            call('vel_thumbs' + overWriteFlag + keepVzFlag + ' ' + thumbfile, shell=True)
        else:
            u.myalert('No data for ' + tiefile)


if __name__ == '__main__':
    main()
