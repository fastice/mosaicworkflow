#!/usr/bin/env python
import os
import sys
import glob
import argparse
from subprocess import call
import utilities as u
import threading
import datetime
import yaml


def getRefreshArgs():
    defaultTracks = ['track-26', 'track-74', 'track-90',
                     'track-112', 'track-141', 'track-170']
    parser = argparse.ArgumentParser(
        description='Make and run ties and makeframetie.py for multiple tracks',
        epilog=(
            'Notes:\n'
            '  1) Run in main directory above tracks; operates on each directory\n'
            '  2) Executes maketies.py -run year1 year2 ... for each track\n'
            '  3) --winter uses tiepoints in .../Sentinel1/tiepoints/'
            'track-winter/tiepoints-winter-Year\n'
            '  4) --tieFiles uses quarterly tiepoint files in tiefile.yaml\n'
            '\nPart of the mosaicworkflow package.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'years', nargs='*', type=int,
        metavar='YEAR',
        help='Years to process (default: current year)'
    )
    parser.add_argument(
        '-toRun', '--toRun', type=str, default=None,
        metavar='TRACKS',
        help=f'Python-formatted list of track directories, e.g. '
             f'"[\'track-26\',\'track-74\']" (default: {defaultTracks})'
    )
    parser.add_argument(
        '-tieFiles', '--tieFiles', type=str, default=None,
        metavar='FILE',
        help='YAML file with time-varying tiepoints'
    )
    parser.add_argument(
        '-tiesOnly', '--tiesOnly', action='store_true',
        help='Only do tiepoints'
    )
    parser.add_argument(
        '-velthumbsOnly', '--velthumbsOnly', action='store_true',
        help='Only do velthumbs'
    )
    parser.add_argument(
        '-phase', '--phase', action='store_true',
        help='Phase tiepoints (no thumbs; implies --tiesOnly)'
    )
    parser.add_argument(
        '-phaseAndOffsets', '--phaseAndOffsets', action='store_true',
        help='NISAR phase+offsets mode: pass -phase to maketies but still run makeframetie.py'
    )
    parser.add_argument(
        '-winter', '--winter', action='store_true',
        help='Use winter tiepoints'
    )
    parser.add_argument(
        '-noPrompt', '--noPrompt', action='store_true',
        help='Run without prompt'
    )
    parser.add_argument(
        '-noQuadFit', '--noQuadFit', action='store_true',
        help='Pass --noQuadFit to tie_script to skip the -deltaBQ quadratic '
             'baseline correction estimate'
    )
    parser.add_argument(
        '--yaml', action='store_true',
        help='Pass --yaml to maketies/tie_script to write rBaseline files in YAML format'
    )
    parser.add_argument(
        '--overWrite', action='store_true',
        help='Pass --overWrite to makeframetie.py to rerun existing products'
    )
    parser.add_argument(
        '--keepVz', action='store_true',
        help='Pass --keepVz to makeframetie.py to retain .vz and .vz.geodat files'
    )
    parser.add_argument(
        '--useSquint', action='store_true',
        help='Pass --useSquint to makeframetie.py (squint heading correction in mosaic3d '
             'and tiepoints -motion)'
    )
    parser.add_argument(
        '-nThreads', '--nThreads', type=int, default=24,
        metavar='N',
        help='Maximum tracks processed in parallel [24]'
    )

    args = parser.parse_args()

    # Validate year range
    for y in args.years:
        if y < 1990 or y > 2030:
            parser.error(f'Invalid year: {y}')

    years = args.years if args.years else [datetime.datetime.now().year]

    toRun = eval(args.toRun) if args.toRun is not None else defaultTracks

    tiesOnly = args.tiesOnly
    velthumbsOnly = args.velthumbsOnly

    flags = {'phaseFlag': '', 'winterFlag': '', 'noQuadFitFlag': '', 'yamlFlag': ''}
    if args.winter:
        flags['winterFlag'] = ' -winter '
    if args.phase:
        tiesOnly = True  # Always noThumbs with phase
        flags['phaseFlag'] = ' -phase '
    if args.phaseAndOffsets:
        # NISAR mode: maketies gets -phase but makeframetie.py still runs
        flags['phaseFlag'] = ' -phase '
    if args.noQuadFit:
        flags['noQuadFitFlag'] = ' -noQuadFit '
    if args.yaml:
        flags['yamlFlag'] = ' --yaml '

    if velthumbsOnly and tiesOnly:
        u.myerror('Cannot use both velthumbsOnly and tiesOnly')

    usePrompt = not args.noPrompt

    return years, toRun, tiesOnly, velthumbsOnly, flags, usePrompt, args.tieFiles, args.overWrite, args.keepVz, args.useSquint, args.nThreads


def clearPendingExcludes(rundir):
    """Fresh-start the SOFT excludes for a tie refresh: remove every
    Exclude.pending under rundir's frame dirs so this refresh re-evaluates all
    frames from scratch. tieScript re-creates Exclude.pending only for frames
    whose baseline fit fails again this run (build_tie_file). Hard Exclude
    files are deliberately left untouched -- those are permanent."""
    removed = 0
    for pend in glob.glob(os.path.join(rundir, '*_*', 'Exclude.pending')):
        try:
            os.remove(pend)
            removed += 1
        except OSError:
            pass
    if removed:
        print(f'refreshties: cleared {removed} stale Exclude.pending under {rundir}')


def runTies(rundir, years, tiesOnly, velthumbsOnly, tieFiles, flags, overWrite, keepVz,
           useSquint=False):
    cwd = os.getcwd()
    # Clear stale soft-excludes before maketies/setuptopstie build the tie plan,
    # so a refresh always starts clean and only genuinely-failing frames end up
    # marked pending again this run.
    clearPendingExcludes(rundir)
    sensor = None
    projectYaml = os.path.join(cwd, 'project.yaml')
    legacyYaml = os.path.join(cwd, 'sensor.yaml')
    if os.path.isfile(projectYaml):
        sensorYaml = projectYaml
    elif os.path.isfile(legacyYaml):
        u.mywarning('sensor.yaml found but project.yaml expected — rename to project.yaml')
        sensorYaml = legacyYaml
    else:
        sensorYaml = None
    if sensorYaml is not None:
        with open(sensorYaml, 'r') as f:
            sensorData = yaml.safe_load(f)
        sensor = sensorData.get('sensor', None)
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
            u.myerror('refreshties.py: unsuported sensor')
    try:
        fout = open(f'{rundir}/stdout', 'w')
        ferr = open(f'{rundir}//stderr', 'w')
    except Exception:
        u.myerror(f'could not open {rundir}/stdout - '
                  'check at right directory level')
    command = 'pushd ' + rundir + ' '
    print(sensor, rundir)
    if not velthumbsOnly:
        command += ' ; maketies -run '
        #
        if tieFiles is not None:
            command += f' -tieFiles={tieFiles} '
        for flag in flags:
            command += flags[flag]
        #
        for year in years:
            command += ' ' + str(year)
    if not tiesOnly:
        if sensor in ['Sentinel1', 'ALOS2', 'NISARTest', 'NISAR']:
            command += ' ; cd tiepoints '
            suffix = ''
        else:
            fullRunDir = f'{cwd}/{rundir}'
            glacierID = fullRunDir.split('track')[0].strip('/').split('/')[-1]
            if len(rundir.split('/')) > 1:
                runTop = rundir.split('/')[0]
            else:
                runTop = '.'
            #
            multiTrack = ['Jak', 'Kang', 'Helheim', 'Zach', 'DJ', 'south1',
                          'Leo', 'Lakes', 'Taku', 'Hubbard', 'Columbia', 'PIG']
            suffix = ''
            # for tsx runTop is the name (eg. Jak), so check if in multiTrack
            # if it is, use the track for the suffix.
            if glacierID in multiTrack:
                suffix = f'-{rundir.split("-")[-1]}'  # track num as suffix
            command += f' ; popd ; cd {runTop}/tiepoints '
        overWriteFlag = ' --overWrite' if overWrite else ''
        keepVzFlag = ' --keepVz' if keepVz else ''
        yamlFlag = flags['yamlFlag']
        squintFlag = ' --useSquint' if useSquint else ''
        for year in years:
            command += f'; makeframetie.py{overWriteFlag}{keepVzFlag}{squintFlag}{yamlFlag} tie_plan{year}{suffix}'
    print(command)
    # print('suffix:',suffix)
    call(command, shell=True, executable='/bin/csh', stdout=fout, stderr=ferr)
    fout.close()
    ferr.close()


def main():
    years, toRun, tiesOnly, velthumbsOnly, flags, usePrompt, tieFiles, overWrite, keepVz, \
        useSquint, nThreads = getRefreshArgs()
    print(years)
    threads = []

    for runfile in toRun:
        thread = threading.Thread(target=runTies,
                                  args=[runfile, years, tiesOnly,
                                        velthumbsOnly, tieFiles, flags, overWrite, keepVz,
                                        useSquint])
        threads.append(thread)
    #
    # prompt to run jobs
    #
    u.runMyThreads(threads, nThreads, 'Refresh Ties ', prompt=usePrompt)


if __name__ == '__main__':
    main()
