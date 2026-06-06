#!/usr/bin/env python3
import utilities as u
import os
import argparse
from datetime import datetime
from subprocess import check_output
import sys
import glob
import yaml
import json

def _readOrbitSensorYaml(orbitDir):
    """Read sensor.*.yaml from an orbit_frame directory; return dict or {}."""
    yamlFiles = glob.glob(os.path.join(orbitDir, 'sensor.*.yaml'))
    if yamlFiles:
        with open(yamlFiles[0]) as f:
            return yaml.safe_load(f) or {}
    return {}


def setuptopstiesUsage():
    """Print usage/help for the script using argparse formatting."""
    parser = argparse.ArgumentParser(
      prog='setuptopstie.py',
      description='Make a tops tiepoint file for specified year',
      epilog='Part of the mosaicworkflow package.')
    parser.add_argument('-winter', '--winter', action='store_true',
                help='Uses winter tiepoints generated from coastal pairs')
    parser.add_argument('-phase', '--phase', action='store_true',
                help='Creates interferogram tiefile')
    parser.add_argument('-tieFiles', '--tieFiles',
                help='YAML file with time varying tiepoints')
    parser.add_argument('-tie_plan_suffix', '--tie_plan_suffix',
                help='Suffix to append to tie_plan filename')
    parser.add_argument('year', nargs='?', help='Year to process')
    parser.print_help()
    exit()


def getSTTArgs():
    """Get arguments using argparse. Accept single-dash long options
    for backward compatibility (e.g. -winter, -tieFiles=...).
    Returns: year (string), winterTies (bool), tie_suffix (str),
             phaseTies (bool), tieFiles (str or None)
    """
    # Preprocess sys.argv so that single-dash long options work with argparse
    argv = []
    for a in sys.argv[1:]:
        if a.startswith('--'):
            argv.append(a)
        elif a.startswith('-') and not (len(a) == 2):
            # treat single-dash long option like --option
            argv.append('--' + a.lstrip('-'))
        else:
            argv.append(a)

    parser = argparse.ArgumentParser(
        prog='setuptopstie.py',
        description='Make a tops tiepoint file for specified year',
      epilog='Part of the mosaicworkflow package.')
    parser.add_argument('-winter', '--winter', action='store_true',
                        help='Uses winter tiepoints generated from coastal pairs')
    parser.add_argument('-phase', '--phase', action='store_true',
                        help='Creates interferogram tiefile')
    parser.add_argument('-tieFiles', '--tieFiles',
                        help='YAML file with time varying tiepoints')
    parser.add_argument('-tie_plan_suffix', '--tie_plan_suffix',
                        help='Suffix to append to tie_plan filename')
    parser.add_argument('year', nargs='?', help='Year to process')

    args = parser.parse_args(argv)

    year = args.year if args.year is not None else -1
    winterTies = bool(args.winter)
    phaseTies = bool(args.phase)
    tie_suffix = args.tie_plan_suffix or ''
    tieFiles = args.tieFiles

    # Validate year
    try:
        if int(year) < 2008 or int(year) > 2035:
            print(f"\n\t\t\033[1;31m *** Invalid Year {year} not in (2008,2035) *** \033[0m")
            setuptopstiesUsage()
    except Exception:
        print(f"\n\t\t\033[1;31m *** Invalid Year {year} not in (2008,2035) *** \033[0m")
        setuptopstiesUsage()

    return str(year), winterTies, tie_suffix, phaseTies, tieFiles


def getOffsize(myPath, myDatFile):
    myDatPath = f'{myPath}/{myDatFile}'
    myVrtPath = myDatPath.replace('.dat', '.vrt')
    # proceed as long as dat or vrt exists
    if os.path.isfile(myDatPath) or os.path.isfile(myVrtPath):
        if not os.path.exists(myVrtPath):
            myOffsets = u.offsets(datFile=myDatFile, myPath=myPath,
                                  verbose=False)
        else:
            myOffsets = u.offsets(vrtFile=myVrtPath, myPath=myPath,
                                  verbose=False)
       
        myOffsets.readOffsetsDat()
        #u.myerror('ffff')
        # fDat = open(myDatFile, 'r')
        # line = fDat.readline().split()
        # fDat.close()
        # ncol, nrow = int(line[2]), int(line[3])
        return myOffsets.nr, myOffsets.na
    # return error if file did not exist
    return -1, '## Missing azimuth.offsets.dat '


def checkOffsets(imageDir):
    ''' Check that the offsets files exists and are non-zero '''
    #
    # Check data files files exist
    #
    azimuthFile = f'{imageDir}/azimuth.offsets'
    azimuthVrt = f'{imageDir}/azimuth.offsets.vrt'
    rangeFile = f'{imageDir}/range.offsets'
    rangeVrt = f'{imageDir}/range.offsets.vrt'
    if not os.path.isfile(azimuthFile) and not os.path.isfile(azimuthVrt):
        return f'## {imageDir} no azimuth file '
    if not os.path.isfile(rangeFile) and not os.path.isfile(rangeVrt):
        return f'## {imageDir} no range file '
    # If only .vrt files present, skip size checks
    if not os.path.isfile(azimuthFile) and not os.path.isfile(rangeFile):
        return ''
    #
    # Now find data file sizes
    #
    azimuthFileSize = os.stat(azimuthFile if os.path.isfile(azimuthFile) else azimuthVrt).st_size
    rangeFileSize = os.stat(rangeFile if os.path.isfile(rangeFile) else rangeVrt).st_size
     
    #
    # Now check the files are the expected size (use size info from azimuth)
    #
    
    ncol, nrow = getOffsize(imageDir, 'azimuth.offsets.dat')
   
    # nrow=err msg if ncol < 0
    if ncol < 0:
        return nrow
    #
    nrowFile = int(azimuthFileSize/(ncol*4))
    # check calculated row size same as expected
    if nrow != nrowFile:
        return f'## invalid azimuth.file size {nrow})'
    nrowFile = int(rangeFileSize/(ncol*4))
    # do check for range now too
    if nrow != nrowFile:
        return(f'##invalid range.file size {nrow}')
    #
    # Success return empty string as blank prefix
    return('')


def getFlags(track, header, year):
    ''' If there are any flags in a flagstrack file, then read them and return
    as flags '''
#
    tiedir = header.split('tie_plan_header')[0]
    flags = []
    flagfile = tiedir+'flags'+track.split('-')[-1]
    if not os.path.isfile(flagfile):
        flagfile += '-' + str(year)

    if os.path.isfile(flagfile):
        flagFile = open(flagfile, 'r')
        for line in flagFile:
            if not line.isspace() and len(line) > 0:
                flags.append(line.strip('\n'))
        flagFile.close()
    return flags


def getMySensorAndTrack():
    ''' Get sensor, track, and directory info'''
    # First, try to read sensor and region from ../project.yaml (preferred) or
    # ../sensor.yaml (legacy fallback) and use them if present
    sensor = None
    region = None
    project_yaml = os.path.abspath(os.path.join(os.getcwd(), '..', 'project.yaml'))
    legacy_yaml  = os.path.abspath(os.path.join(os.getcwd(), '..', 'sensor.yaml'))
    if os.path.isfile(project_yaml):
        sensor_yaml = project_yaml
    elif os.path.isfile(legacy_yaml):
        u.mywarning('sensor.yaml found but project.yaml expected — rename to project.yaml')
        sensor_yaml = legacy_yaml
    else:
        sensor_yaml = None
    try:
        if sensor_yaml is not None:
            with open(sensor_yaml, 'r') as sf:
                sdata = yaml.load(sf, Loader=yaml.FullLoader)
            if isinstance(sdata, dict):
                for k, v in sdata.items():
                    lk = str(k).lower()
                    if lk == 'sensor':
                        sensor = str(v)
                    elif lk == 'region':
                        region = str(v)
    except Exception:
        # If anything goes wrong reading the file, fall back to path-based detection
        sensor = None
        region = None
    # default Sentinel unless other sensor name in path (only if yaml didn't provide)
    pwd = os.getcwd()
    if sensor is None:
        sensor = 'Sentinel'
        if 'TSX' in pwd:
            sensor = 'TSX'
        # pieces = pwd.split('/')
        elif 'CSK' in pwd:
            sensor = 'CSK'
        elif 'ALOS2' in pwd:
            sensor = 'ALOS2'
        elif 'NISAR' in pwd:
            sensor = 'NISARTest'
        elif not ('Sentinel' in pwd):
            u.myerror('Not in a Sentinel or TSX directory')
    
    # Not sure why all this is necessary, but try to find root dir for tiepoints based on path (looking for sensor name in path
    rootdir = ''
    pieces = pwd.split('/')
    for piece in pieces:
        rootdir = rootdir + piece + '/'
        found = False
        for sname in ['TSX', 'Sentinel', 'CSK', 'ALOS', 'NISARTest']:
            if sname in piece:
                found = True
                break
        if found:
            break
    if not found:
        rootdir = os.path.dirname(pwd) + '/'
    # Get track and tracknum from path
    track = os.getcwd().split('/')[-1]
    tracknum = track.split('-')[-1]
    print(track, tracknum, file=sys.stderr)
    #
    # If region was not supplied via sensor.yaml, check for a region file,
    # otherwise default to 'greenland'
    if region is None:
        if os.path.exists(f'{rootdir}/region'):
            with open(f'{rootdir}/region') as fp:
                region = fp.readline().strip()
        else:
            region = 'greenland'
    #
    return sensor, track, tracknum, rootdir, region, sensor_yaml


def getSensorInfoAndHeader(sensor, year, sensor_yaml=None):
    ''' Get sensor specific info '''
    # Read maxDays, header, and framePattern from project.yaml (or sensor.yaml fallback) if available
    yaml_maxDays = None
    yaml_header = None
    framePattern = '*'
    if sensor_yaml is not None and os.path.isfile(sensor_yaml):
        try:
            with open(sensor_yaml, 'r') as sf:
                sdata = yaml.load(sf, Loader=yaml.FullLoader)
            if isinstance(sdata, dict):
                for k, v in sdata.items():
                    lk = str(k).lower()
                    if lk == 'maxdays' or lk == 'max_days':
                        try:
                            yaml_maxDays = int(v)
                        except Exception:
                            yaml_maxDays = None
                    elif lk == 'header':
                        yaml_header = str(v)
                    elif lk == 'framepattern' or lk == 'frame_pattern':
                        framePattern = str(v)
        except Exception:
            yaml_maxDays = None
            yaml_header = None
            framePattern = '*'
    try:
        # determine if header exist for different dir structures
        if yaml_header is not None:
            header = yaml_header
            maxDays = 37  # default; overridden below if yaml_maxDays set
        elif sensor in ['Sentinel', 'NISARTest', 'NISAR']:
            header = 'tiepoints/tie_plan_header'
            maxDays = 37
        elif sensor == 'ALOS2':
            header = 'tiepoints/tie_plan_header'
            maxDays = 42
        elif sensor == 'TSX' or sensor == 'CSK':
            header = str(check_output('pwd', shell=True), 'utf-8').strip('\n')
            header = header.split('track')[0]+'/tiepoints/tie_plan_header'
            maxDays = 34
        else:
            raise NameError()
        # no header file
        print(header, file=sys.stderr)
        if os.path.isfile(header+str(year)):
            header = header+str(year)
        elif not os.path.isfile(header):
            raise NameError()
    except Exception:
        u.myerror(' Error NO tie_plan_header file in tiepoints - check in '
                  f'right directory {sensor}')
    return header, yaml_maxDays if yaml_maxDays is not None else maxDays, framePattern


def readGeoDate(geoFile):
    ''' Parse date from geodatnxa.in'''
    with open(geoFile, 'r') as fp:
        if '.in' in geoFile:
            for line in fp:
                if 'date' in line:
                    myDate = datetime.strptime(line.split(":")[1].strip(),
                                               "%d %b %Y")
                    return myDate
        elif '.geojson' in geoFile:
            geojson = json.load(fp)
            date = geojson["properties"]["Date"]
            return datetime.strptime(date, '%Y-%m-%d')
    u.myerror('readGeodate: date not found')


def parsePairInfo(myPairInfo):
    ''' read a pairinfo file '''
    fp = open(myPairInfo)
    pieces = fp.readline().split()
    orbit1, orbit2 = pieces[0:2]
    date1 = datetime.strptime(pieces[2], "%Y-%m-%d")
    date2 = datetime.strptime(pieces[3], "%Y-%m-%d")
    nlr, nla = pieces[4:6]
    fp.close()
    return date1, date2, orbit1, orbit2, nlr, nla


def writePairInfo(myDir, orbit1, orbit2, date1, date2, nlr, nla):
    ''' write a pairinfo file '''
    fp = open(os.path.join(myDir, f'{orbit1}.{orbit2}.pairinfo'), 'w')
    print(f'{orbit1}  {orbit2}  {date1.strftime("%Y-%m-%d")} '
          f'{date2.strftime("%Y-%m-%d")} {nlr} {nla}', file=fp)
    fp.close


def parsePairFromStrack(strackFile, myDir):
    ''' Use strack input file to find pair data '''
    # print(strackFile)
    try:
        if os.stat(strackFile).st_size == 0:
            u.myerror(f'Zero length strack file {strackFile}')
        fp = open(strackFile)
        for line in fp:
            if ';' in line:
                continue
            if 'image1 ' in line:
                image1 = line.split()[3]
                orbit1 = int(os.path.basename(image1).split('_')[0])
            elif 'image2 ' in line:
                image2 = line.split()[3]
                myDir2 = os.path.dirname(image2).replace('../', '')
                orbit2 = int(os.path.basename(image2).split('_')[0])
            elif 'intgeodat' in line:
                geodat = line.split()[3]
                break
        fp.close()
    except Exception:
        u.myerror('parsePairFromStrack: Could not parse info from strack '
                  f'file: {strackFile}')
    # finish parsing
    # print(geodat)
    suffix = f".{geodat.split('.')[-1]}"
    pieces = geodat.split('geodat')[1].split(suffix)[0].split('x')
    nlr, nla = int(pieces[0]), int(pieces[1])
    # dates
    date1 = readGeoDate(os.path.join(myDir, os.path.basename(geodat)))
    date2 = readGeoDate(os.path.join(myDir2, os.path.basename(geodat)))
    #
    return date1, date2, orbit1, orbit2, nlr, nla


def getPairs(year, phaseTies, winterTies, region, sensor, maxDays, framePattern='*'):
    ''' '''
    myDirs = sorted(glob.glob(f'*_{framePattern}'))
    nDaysList = []
    pairs = {}
    # default
    firstDate, lastDate = datetime(year, 1, 1), datetime(year, 12, 31)
    if phaseTies:
        # Split to get all phaseties from same winter
        if region in ['greenland']:
            firstDate, lastDate = datetime(year-1, 6, 1), datetime(year, 5, 31)
    #
    # Look through data
    for myDir in myDirs:
        myPair = sorted(glob.glob(f'{myDir}/*.pairinfo'))
        # get pairinfo from pairinfo file if it exists
        frame = myDir.split('_')[-1]
        if len(myPair) > 0:
            date1, date2, orbit1, orbit2, nlr, nla = parsePairInfo(myPair[0])
        else:
            strackFile = sorted(glob.glob(f'{myDir}/strackin.*'))
            if len(strackFile) > 0:
                # Print to sys.stderr to avoid error in maketies
                print(strackFile[0], file=sys.stderr)
                date1, date2, orbit1, orbit2, nlr, nla = parsePairFromStrack(
                    strackFile[0], myDir)
                writePairInfo(myDir, orbit1, orbit2, date1, date2, nlr, nla)
            else:
                date1, date2, orbit1, orbit2, nlr, nla = 6*[None]
        # skip this dir if no date
        if date1 is None:
            continue
        #
        if date1 >= firstDate and date1 <= lastDate and date1 != date2:
            nD = (date2-date1).days
            if nD < 0:
                nD = -nD
                print('**** warning: negative nDays - changing sign***',
                      file=sys.stderr)
            # for NISAR, read per-product maxDays from sensor.*.yaml; others use global
            if 'NISAR' in sensor:
                orbitParams = _readOrbitSensorYaml(myDir)
                effectiveMaxDays = orbitParams.get('maxDays', maxDays)
            else:
                effectiveMaxDays = maxDays
            if nD > effectiveMaxDays:
                continue
            if nD not in nDaysList:
                nDaysList.append(nD)
                pairs.update({nD: []})
            pairs[nD].append({'date1': date1, 'date2': date2, 'orbit1': orbit1,
                              'orbit2': orbit2, 'frame': frame, 'nlr': nlr,
                              'nla': nla})
    # sort by date
    for nD in nDaysList:
        pairs[nD] = sorted(pairs[nD], key=lambda k: k['date1'])
    #
    return pairs, nDaysList


def copyHeader(fout, header, winterTies, phaseTies, sensor, rootdir, year):
    fhead = open(header, 'r')
    for line in fhead:
        if 'extraties' in line and '#' not in line:
            if winterTies:
                if 'TSX' in sensor or 'CSK' in sensor:
                    u.myerror('-winter does not work for TSX at this time')
                #
                pieces = line.split('=')
                tiedir = f'{rootdir}/tiepoints/track-winter/tiepoints-' \
                    f'winter-{year}.culled'
                if not os.path.isfile(tiedir):
                    u.myerror(f'Invalid tiedir {tiedir}--for winter in {year}')
                line = pieces[0] + ' = ' + tiedir+'\n'
            elif phaseTies:
                line = line.replace('.year', f'.{year}')
            extraTies = line.split('=')[1].strip()
        print(line, file=fout, end='')
    print('', file=fout)
    fhead.close()
    print('start', file=fout)
    return extraTies

def getPrefix(imageDir, phaseTies):
    prefix = ''
    noUse = False
    if os.path.isfile(imageDir+'/Exclude'):
        prefix = '## Exclude '
    # if not Exclude, check offsets for existence/size
    elif os.path.isfile(imageDir+'/Special'):
        prefix = '## Special '
    else:
        if not phaseTies:
            prefix = checkOffsets(imageDir)
        else:
            # currently no checking for phase files, could add later
            prefix = ''
    #
    # print warning if prefix has text
    if len(prefix) > 1:
        noUse = True
        print(f'*** Warning: {prefix[2:]} --- {imageDir} ***', file=sys.stderr)
    #
    return prefix, noUse


def tieCodes(phaseTies, track, sensor=''):
    ''' codes for phase and offsets'''
    if phaseTies and 'NISAR' in sensor:
        # yaml flat-earth mode: tiepoints -yaml, phase + offsets, no changeflat
        myCodes = {'usePre': 'pys', 'useSuf': '', 'tiePre': 'py',
                   'tieSuf': '', 'track': f' {track} '}
    elif phaseTies:
        myCodes = {'usePre': 'pps', 'useSuf': 'uw', 'tiePre': 'pn',
                   'tieSuf': 'uw', 'track': f' {track} '}
    else:
        myCodes = {'usePre': '', 'useSuf': '', 'tiePre': 'n', 'tieSuf': '',
                   'track': ''}
    return myCodes


def getTieInterval(pair, tieIntervalCodes, year):
    ''' Get interval code based on pair date '''
    if tieIntervalCodes is None:
        return f'1y{year}'
    #
    for key in tieIntervalCodes:
        if pair['date1'] >= tieIntervalCodes[key]['firstDate'] and \
                pair['date1'] <= tieIntervalCodes[key]['lastDate']:
            return key


def processNDays(fout, pairs, track, tracknum, nDays, year, flags, phaseTies,
                 tieIntervalCodes, sensor=''):
    ''' write tie_plan for each use,  .... tiefile orbitdnDays....'''
    #
    codes = tieCodes(phaseTies, tracknum, sensor=sensor)
    #
    # setup new subdir
    print(f'\nsubdir {track} nDays={nDays}', file=fout)
    # print(f'\ttiefile 1y{year}', file=fout)
    # print(f'\t\t{tracknum} bedrock', file=fout)
    # print(f'\n\tuse 1y{year}', file=fout)
    # print(f'\n\tuse 1y{year}', file=fout)
    #
  
    #
    # write each pair to be processed
    currentTieInterval = None
    for pair in pairs:
        prefix, noUse = getPrefix(f'{pair["orbit1"]}_{pair["frame"]}',
                                  phaseTies)
        #print(pair)
        #u.myerror('fff=f')
        if noUse:
            pair['orbit1'] = None
        #
       
        intervalCode = getTieInterval(pair, tieIntervalCodes, year)
        
        if intervalCode != currentTieInterval:
            currentTieInterval = intervalCode
            print(f'\n\tuse {currentTieInterval}', file=fout)
            # print flags
            if len(flags) > 0:
                for flag in flags:
                    print('\t'+flag, file=fout)

        # print(pair)
        print(f'\t\t{prefix}{codes["usePre"]}{codes["track"]} '
              f'{pair["orbit1"]}_{pair["frame"]} {codes["useSuf"]}', file=fout)
    print(f'\n#\n#  {nDays}Day  Pairs\n#', file=fout)
    #
   
    lastOrb = 0
    orbGroups = []
    #
    # build the track-orbitdnDays
    for pair in pairs:
        if pair['orbit1'] is not None:
            if pair['orbit1'] != lastOrb:
                print(f'\ttiefile {pair["orbit1"]}d{nDays}', file=fout)
                orbGroups.append(f'{pair["orbit1"]}d{nDays}')
            print(f'\t\t{codes["tiePre"]} '
                  f'{tracknum} {pair["orbit1"]}_{pair["frame"]} '
                  f'{codes["tieSuf"]}', file=fout)
            lastOrb = pair["orbit1"]
    return orbGroups


def processAll(fout, allOrbGroups, year, tracknum):
    ''' Create all tiefile with all data for year'''
    print(f'\n\ttiefile all{year}', file=fout)
    count = 0
    for orbGroups in allOrbGroups:
        for orbGroup in orbGroups:
            print(f'\t\tc {tracknum} {orbGroup}', file=fout)
            count += 1
    print('end', file=fout)
    return count


def setupExtraTies(fout, tieFilesData, track, tracknum, year, extraTies):
    '''
    Create the bedrock ties for multiple periods
    '''
    if tieFilesData is None:
        print(f'\nsubdir {track}', file=fout)
        print(f'\textraties={extraTies}', file=fout)
        print(f'\ttiefile 1y{year}', file=fout)
        print(f'\t\t{tracknum} bedrock', file=fout)
        
        print(f'\nsubdir {track}', file=fout)
        print(f'\textraties={extraTies}.480', file=fout)
        print(f'\ttiefile 1y{year}_480', file=fout)
        print(f'\t\t{tracknum} bedrock', file=fout)
        print('# force extra ties reset', file=fout)
        print(f'\nsubdir {track}', file=fout)
        print(f'\textraties={extraTies}', file=fout)
        return None
    # Handle case where there are multiple tie files
    yearStart = datetime(int(year), 1, 1, 0, 0, 0)
    yearEnd = datetime(int(year), 12, 31, 12, 59, 59)
    tieIntervalCodes = {}
    n = 1
    for key in tieFilesData:
        myCode = f'{n}y{year}'
        if tieFilesData[key]['firstDate'] >= yearStart and \
                tieFilesData[key]['firstDate'] <= yearEnd:
            # Keep dates for code
            # print('track ', track)
            tieIntervalCodes[myCode] = \
                {'firstDate': tieFilesData[key]['firstDate'],
                 'lastDate': tieFilesData[key]['lastDate']}
            #
            print(f'\nsubdir {track}', file=fout)
            print(f'\textraties={tieFilesData[key]["tiefile"]}', file=fout)
            print(f'\ttiefile {n}y{year}', file=fout)
            print(f'\t\t{tracknum} bedrock', file=fout)
            n += 1
    return tieIntervalCodes


#
def readTieFiles(tieFiles):
    if tieFiles is not None:
        with open(tieFiles, 'r') as fp:
            tieFilesData = yaml.load(fp, Loader=yaml.FullLoader)
        for key in tieFilesData:
            tieFilesData[key]['firstDate'] = \
                datetime.strptime(tieFilesData[key]['firstDate'], '%Y-%m-%d')
            tieFilesData[key]['lastDate'] = \
                datetime.strptime(tieFilesData[key]['lastDate'], '%Y-%m-%d')
    else:
        tieFilesData = None
    return tieFilesData
#
# ****************************************
#


def main():
    #
    # get sensor info
    sensor, track, tracknum, rootdir, region, sensor_yaml = getMySensorAndTrack()
    print(sensor, track, tracknum, rootdir, region, sensor_yaml, file=sys.stderr)
   
    #
    # Usage
    # ------ modify this for other sensors
    year, winterTies, tie_suffix, phaseTies, tieFiles = getSTTArgs()
    print(year, winterTies, tie_suffix, phaseTies, tieFiles, file=sys.stderr)
    tieFilesData = readTieFiles(tieFiles)
  
    #
    # get tie points header
    header, maxDays, framePattern = getSensorInfoAndHeader(sensor, year, sensor_yaml)
    #
    # get flags
    flags = getFlags(track, header, year)
    #
    # find pairs for years, sort into dicts, with indxing by nDays
    pairs, nDays = getPairs(int(year), phaseTies, winterTies, region,
                            sensor, maxDays, framePattern)
    #print(pairs, nDays)
    #u.myerror('ffff')
    # print(nDays)
    #
    # Create tie file
    if len(tie_suffix) < 1:
        outfile = f'tie_plan{year}'
    else:
        outfile = f'tie_plan{year}-{tie_suffix}'
    fout = open(outfile, 'w')
    #
    # Copy tie_plan_header
    extraTies = copyHeader(fout, header, winterTies, phaseTies, sensor,
                           rootdir, year)
    #
    #print(pairs)
    if bool(pairs):
        tieIntervalCodes = setupExtraTies(fout, tieFilesData, track, tracknum,
                                          year, extraTies)
    #
    allOrbGroups = []
   # print(tieIntervalCodes)
    #u.myerror('ffff')
    for nDay in nDays:
        orbGroups = processNDays(fout, pairs[nDay], track, tracknum, nDay,
                                 year, flags, phaseTies, tieIntervalCodes,
                                 sensor=sensor)
        
        allOrbGroups.append(orbGroups)
    #
    #
    count = processAll(fout, allOrbGroups, year, tracknum)
    # output this for maketies to use
    print(count)
    fout.close()


if __name__ == '__main__':
    main()
