#!/usr/bin/env python3
import datetime
import os
import sys
import argparse
import yaml
from subprocess import call, check_output
import utilities as u

#
# check the tie_plans to see if data were generated - return list with valid
# data check tie_plans to see if a subdir exists
#


def getAll(tieDir, track):
    '''
    Check the tie_plans to see if data were generated.
    Return list with valid  data check tie_plans to see if a subdir exists
    '''
    if len(track) > 0:
        files = u.dols('ls -d '+tieDir+'/tie_plan20??-'+track)
    else:
        files = u.dols('ls -d '+tieDir+'/tie_plan20??')
    #
    years = []
    for file in files:
        fPlan = open(file, 'r')
        for line in fPlan:
            # if subdir found, extract date, and add it to the list
            if 'subdir' in line:
                year = file.split('/')[-1]
                year = year.split('tie_plan')[-1]
                year = year.split('-')[0]
                years.append(year)
                break
        fPlan.close()
    return years

#
# make tie_planAll
#


def makeAll(tieAll, tieDir, tracknum, mtFlag):
    '''
     Open header, copy to fAll
    '''
    fHeader = open(tieDir+'/tie_plan_header', 'r')
    fAll = open(tieAll, 'w')
    for line in fHeader:
        print(line, end='', file=fAll)
    print('start', file=fAll)
    print('#\n#', file=fAll)
    print('subdir track-'+str(tracknum), file=fAll)
    print('\ttiefile all', file=fAll)
    fHeader.close()
    #
    # Get all tie_plans (this is so if a subset of years requested, all prior
    # years get inlcuded to)
    #
    if mtFlag:
        allEntries = getAll(tieDir, str(tracknum))
    else:
        allEntries = getAll(tieDir, '')
    # add each year to the tie_plan
    for allEntry in allEntries:
        print('\t\tc '+str(tracknum)+' all'+str(allEntry), file=fAll)
    print('end\n', file=fAll)
    fAll.close()


def makeTiesArgs(year):
    maxYear = int(datetime.date.today().strftime("%Y"))
    defaultYears = list(range(year, maxYear+1))
    # used for TSX, keeping in for future use
    multiTrack = ['Jak', 'Kang', 'Helheim', 'Zach', 'DJ', 'south1', 'Leo',
                  'Lakes', 'Taku', 'Hubbard', 'Columbia', 'PIG']

    parser = argparse.ArgumentParser(
        description='Make TSX and Sentinel tie plans',
        epilog=(
            'Notes:\n'
            '  1) Works for TSX and Sentinel1 - run from track-XXX dir\n'
            '     (assumes TSX or Sentinel1 in path name)\n'
            f'  2) For {[x for x in multiTrack]}, run for each track directory'
            ' (i.e., Jak/track-5)\n'
            f'  3) For {[x for x in multiTrack]}, appends track number to names'
            ' to keep multiple tracks separate\n'
            '  4) --winter overrides extraties in header file with winterYear\n'
            '     version in Sentinel1/tiepoints/track-winter\n'
            '  5) Add an "Exclude" file to a directory to skip it\n'
            '     (will ## entry in tiefile)\n'
            '  6) If flagnnnn (where nnnn=track) exists, adds flags to use\n'
            '     statement (e.g. rparams_flags=-quadB)\n'
            '  7) Checks size and existence of offset files - if error, will\n'
            '     ## with error message\n'
            '  8) Handle special cases with tie_planSpecial[-track] and add\n'
            '     empty "Special" file to affected directory\n'
            '\nPart of the mosaicworkflow package.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'years', nargs='*', type=int,
        metavar='YEAR',
        help=f'Years to process (default: {year} to present)'
    )
    parser.add_argument(
        '-run', '--run', action='store_true',
        help='Run the refreshTies script to update tie points'
    )
    parser.add_argument(
        '-winter', '--winter', action='store_true',
        help='Use winter tiepoints generated from coastal pairs'
    )
    parser.add_argument(
        '-phase', '--phase', action='store_true',
        help='Create tiefiles for phase tiepoints'
    )
    parser.add_argument(
        '-tieFiles', '--tieFiles', type=str, default=None,
        metavar='FILE',
        help='YAML file with time-varying tiepoints'
    )

    args = parser.parse_args()

    years = args.years if args.years else defaultYears
    print(years, file=sys.stderr)
    return years, args.run, args.winter, multiTrack, args.phase, args.tieFiles


def getSensorTrackInfo():
    #
    # assumes in SensorXX/region/track-tracknum
    cwd = os.getcwd()
    track = cwd.split('/')[-1]
    # e.g., Jak, ikeq
    region = cwd.split('/')[-2]
    tracknum = track.split('-')[-1]
    #
    # figure out sensor: check project.yaml one level above cwd (preferred),
    # fall back to sensor.yaml, then fall back to checking the cwd path
    parentDir = os.path.dirname(cwd)
    projectYaml = os.path.join(parentDir, 'project.yaml')
    legacyYaml = os.path.join(parentDir, 'sensor.yaml')
    if os.path.isfile(projectYaml):
        sensorYaml = projectYaml
    elif os.path.isfile(legacyYaml):
        u.mywarning('sensor.yaml found but project.yaml expected — rename to project.yaml')
        sensorYaml = legacyYaml
    else:
        sensorYaml = None
    sensor = None
    if sensorYaml is not None:
        with open(sensorYaml, 'r') as f:
            sensorData = yaml.safe_load(f)
        sensor = sensorData.get('sensor', None)
        print(f'sensor from {sensorYaml}: {sensor}')
    if sensor is None:
        if 'TSX' in cwd:
            sensor = 'TSX'
        elif 'CSK' in cwd:
            sensor = 'CSK'
        elif 'Sentinel' in cwd:
            sensor = 'Sentinel1'
        elif 'ALOS2' in cwd:
            sensor = 'ALOS2'
        elif 'NISAR' in cwd:
            sensor = 'NISAR'
        else:
            u.myerror('invalid sensor make sure directory has TSX, or Sentinel1 '
                      'in path')
    #
    print('sensor = ', sensor)
    # get start year (by default run everything starting this year)
    yearDict = {'CSK': 2018, 'TSX': 2009, 'Sentinel1': 2015, 'ALOS2': 2020,
                'NISARTest': 2025, 'NISAR': 2025}
    if sensor in ['Sentinel1', 'ALOS2', 'NISARTest', 'NISAR']:
        tieDir = cwd + '/tiepoints'
    elif sensor == 'TSX' or sensor == 'CSK':
        tieDir = os.path.join(cwd.split('track')[0], 'tiepoints')
    #
    return sensor, region, track, tracknum, tieDir, yearDict[sensor]


def main():
    #
    # Figure out tie directory
    sensor, region, track, tracknum, tieDir, year = getSensorTrackInfo()
    print(sensor, region, track, tracknum, tieDir, year)
    #
    # get Args
    years, runTies, winterTies, multiTrack, phaseTies, tieFiles = \
        makeTiesArgs(year)
    # only use winter flag for coastal
    winterFlag = [' ', ' -winter '][winterTies and sensor == 'Sentinel1']
    phaseFlag = [' ', ' -phase '][phaseTies and
                                  sensor in ['Sentinel1', 'ALOS2', 'NISAR', 'NISARTest']]
    #
    if not os.path.isdir(tieDir+'/old'):
        try:
            os.mkdir(tieDir+'/old')
        except Exception:
            u.myerror(f'Error creating old directory  {tieDir}/old - make '
                      'sure running in correct directory')
    #
    # suffix for multi track (used to distinguish multiple tie files -
    # one for each track)
    suffix = ''
    if region in multiTrack:
        suffix = '-' + str(tracknum)
    #
    # use to combine tiefiles for all years
    tieAll = os.path.join(tieDir, f'tie_planAll{suffix}')
    tieSpecial = os.path.join(tieDir, f'tie_planSpecial{suffix}')
    #
    # this file is generic refreshTies - will get overwritten each time -
    # just calls different tie_scripts
    runFile = os.path.join(tieDir, f'refreshTies{suffix}')
    fRun = open(runFile, 'w')
    print('#', file=fRun)
    #
    if tieFiles is not None:
        tiesArg = f'-tieFiles={tieFiles}'
    else:
        tiesArg = ''
    #
    # Process each year in list
    #
    for year in years:
        print(year)
        if region in multiTrack:
            # note winterFlag set above, and is blank if not needed
            command = f'setuptopstie -tie_plan_suffix={tracknum} ' \
                f'{tiesArg} {winterFlag} {phaseFlag} {year}'
            #
            print(command)
            nImages = int(str(check_output(command, shell=True),
                              'utf-8').strip('\n'))
            tie_plan = f'tie_plan{year}-{tracknum}'
        else:
            # note winterFlag set above, and is blank if not needed
            command = 'setuptopstie ' \
                f'{tiesArg}  {winterFlag} {phaseFlag} {year}'
            print(command)
            try:
                rValue = str(check_output(command, shell=True),
                             'utf-8').strip('\n')
                print(rValue)
                nImages = int(rValue)
            except Exception:
                u.myerror(f'Error processing: \n {command:s}\n returned:\n '
                          f'{rValue:s}\n')
            tie_plan = 'tie_plan' + str(year)
            print('nImages -- ', nImages)
            #
        tie_plan1 = f'{tieDir}/{tie_plan}'
        tie_planOld = f'{tieDir}/old/{tie_plan}'
        print(tie_plan1)
        print(tie_plan)
        print(tie_planOld)
        if os.path.isfile(tie_plan1):
            os.rename(tie_plan1, tie_planOld)
        os.rename(tie_plan, tie_plan1)
        # only add entries for cases where there were images
        if nImages > 0:
            print('# '+str(year)+'\ntie_script ', tie_plan, file=fRun)
    #
    # make tie_planAll
    #
    makeAll(tieAll, tieDir, tracknum, region in multiTrack)
    #
    print('# All \ntie_script '+tieAll.split('/')[-1], file=fRun)
    #
    print(tieSpecial)
    if os.path.isfile(tieSpecial):
        print('# Special \ntie_script '+tieSpecial, file=fRun)
        #
    print('# end', file=fRun)
    fRun.close()
    #
    # Run script if --run flag used
    #
    if runTies:
        u.pushd(tieDir)
        call(f'csh {runFile}', shell=True)
        u.popd()


if __name__ == '__main__':
    main()
