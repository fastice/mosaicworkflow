#!/usr/bin/env python3
import glob
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from shutil import copyfile
from subprocess import call
import utilities as u
from datetime import datetime
import threading
import yaml

# A vel_thumb_header is named vel_thumb_header_<firstFrame>dash<lastFrame>.
# Anything else sharing the prefix (editor backups, .save copies) is not a
# header -- matching on the full name keeps a stray file from taking the whole
# run down, since buildResolutionDictionary() runs over every track before any
# velocity directory is touched.
headerPattern = re.compile(r'^vel_thumb_header_(\d+)dash(\d+)$')


LOGFILE = 'makevelnoclean.log'
logLock = threading.Lock()


def logMessage(message):
    ''' Append a timestamped line to LOGFILE in the run directory. Excluded
    frames are recorded here rather than on the console: there can be many of
    them and they are an expected outcome, not something to warn about. The
    lock keeps concurrent setup threads from interleaving lines. '''
    with logLock:
        with open(LOGFILE, 'a') as fp:
            print(f'{datetime.now()}: {message}', file=fp)


def frameExcluded(velDir):
    ''' True when the frame holding velDir carries a hard Exclude. The file
    lives in the frame dir, one level above velocity/ (Exclude.pending is a
    soft flag and is deliberately not honoured here). '''
    return os.path.exists(f'{os.path.dirname(velDir)}/Exclude')


def removeProducts(velDir):
    ''' Delete the mosaic products in velDir, so an excluded frame cannot leave
    a stale velocity behind for the mosaic to pick up. Returns the number
    removed. Inputs (runOff, inputFile) and the run logs are left alone. '''
    removed = 0
    for product in glob.glob(f'{velDir}/mosaicOffsets.*'):
        try:
            os.remove(product)
            removed += 1
        except OSError as error:
            logMessage(f'{velDir}: could not remove {product}: {error}')
    return removed


def fixError(saptr, srptr):
    ''' Remove prior link to clean out errors,
    a new pointer will be generated'''
    for link in [saptr, srptr]:
        if os.path.islink(link):
            os.remove(link)


def getOffsetName(mainDir, suffix):
    name = {'da': '*.interp.da', 'dr': '*.interp.dr'}
    altName = {'da': '*.da.interp', 'dr': '*.dr.interp'}
    files = u.globOffsetProducts(f'{mainDir}/{name[suffix]}')
    if len(files) < 1:
        # pre-2020 name order, never written in tiff mode
        files = sorted(glob.glob(f'{mainDir}/{altName[suffix]}'))
    return files


def subOffsets(line, mainDir, fast=False):
    ''' handle offsets'''
    if not fast:
        # hack to deal with file naming started in F 2018 (July first month)
        az = getOffsetName(mainDir, 'da')
        # removing hack because all data should have been culled and renamed
        # if len(az) == 0:
        #     u.myerror(f'{mainDir}: {az} missing, probably needs reculling')
        # az=u.dols('ls '+mainDir+'/*.da.interp')
        # print(mainDir,az)
        # this will filter other dirs with false matches eg test/velocity
        try:
            azname = az[0].split('/')[-1]
        except Exception:
            u.mywarning(f'Skipping {mainDir}; no azimuth file ')
            return None
        # hack to deal with file naming convention started in Jul 2018
        rg = getOffsetName(mainDir, 'dr')
        # if len(rg) == 0: Removed
        #     rg=u.dols('ls '+mainDir+'/*.dr.interp')
        rgname = rg[0].split('/')[-1]
        vrtFileInterp = az[0].replace('.da', '.vrt')
        # If not vrt, copy .dat files
        if not os.path.exists(vrtFileInterp):
            copyfile(f'{mainDir}/azimuth.offsets.dat', az[0]+'.dat')
            copyfile(f'{mainDir}/azimuth.offsets.dat', rg[0]+'.dat')
            #
            # sa = azname.replace('.cull.da.interp','.sa')
            # sr = rgname.replace('.cull.dr.interp','.sr')
            saptr = azname + '.sa'
            srptr = rgname + '.sr'
            u.pushd(mainDir)
            fixError(saptr, srptr)
            if not os.path.exists(saptr):
                os.symlink('azimuth.offsets.slow.sa', saptr)
                os.symlink('range.offsets.slow.sr', srptr)
                # call('ln -s '+sa+' '+saptr, shell=True)
                # call('ln -s '+sr +' '+srptr, shell=True)
            u.popd()
    else:

        azname = 'fast/azimuth.offsets.fast'
        rgname = 'fast/range.offsets.fast'
        if not os.path.exists(f'{mainDir}/fast/azimuth.offsets.fast.dat') and \
                not os.path.exists(f'{mainDir}/fast/offsets.fast.vrt'):
            copyfile(mainDir+'/azimuth.offsets.dat',
                     mainDir+'/fast/azimuth.offsets.fast.dat')
            copyfile(mainDir+'/azimuth.offsets.dat',
                     mainDir+'/fast/range.offsets.fast.dat')
    line = line.replace('azimuth.offsets', azname)
    line = line.replace('range.offsets', rgname)
    return line


def runRunFile(rundir):
    # runOff was (re)written for this run just before this call, so its mtime
    # is a start-of-run timestamp taken from the same filesystem clock that
    # will stamp the products -- unlike time.time() here, which is the client
    # clock and can be skewed against a network volume by enough to make a
    # freshly written product look stale.
    cutoffTime = os.path.getmtime(f'{rundir}/runOff')
    fout = open(rundir+'/stdout', 'w')
    ferr = open(rundir+'/stderr', 'w')
    command = f'cd {rundir}; csh runOff'
    status = call(command, shell=True, stdout=fout, stderr=ferr)
    fout.close()
    ferr.close()
    #
    # Only clean up after a run that produced something. A failed run wrote
    # nothing, so every file predates cutoffTime and the sweep below would
    # delete the products the rerun was meant to replace. runOff ends with
    # the mosaic3d call, so csh exits with mosaic3d's own status.
    if status != 0:
        u.mywarning(f'runOff failed in {rundir} (exit {status}) - keeping '
                    'existing products')
        return status
    # error and vz bands are never kept; both the binary (.ex) and GeoTIFF
    # (.ex.tif) spellings match. An infinite cutoff removes every match
    # regardless of age, which is what the old `rm *.e* *.vz*` did.
    for pattern in ('mosaicOffsets.e*', 'mosaicOffsets.vz*'):
        u.removeStaleFiles(f'{rundir}/{pattern}', float('inf'))
    # Products of an earlier run that this one did not overwrite -- chiefly
    # the other output format (binary mosaicOffsets.vx/.geodat left beside a
    # new mosaicOffsets.vx.tif/.vrt), plus anything else stale such as a
    # legacy mosaicOffsets.report.
    u.removeStaleFiles(f'{rundir}/mosaicOffsets.*', cutoffTime, verbose=True)
    return status


def makevelnocleanUsage():
    ''' Usage statement '''
    print('\n\t\t\033[1;34m Run Nocull Velocities (or redo velocity)\033[0m\n')
    print('\033[1m\n Reproduce velocity directories with no culling in '
          'corresponding velocity_nocull directory \n\n\tUsage: ')
    print('\t\tmakevelnoclean.py -help -reset  -threads=threads -noprompt '
          '-useHeader '
          '-quad -deltaBp \n\t\t-firstdate=YYYY:MM:DD -lastdate=YYYY:MM:DD '
          '-redoculled -reSize -resolution="X0 Y0 XS YS DX DY" '
          '-toRun="[\'track-X\', \'track-Y\',...]"\n\twhere\n')
    print('\t\t-reset\t\tforces reprocessing [default is only create  '
          'in new cases]')
    print('\t\t-threads\tnumber of threads (<= 30) to use [12]')
    print('\t\t-noprompt\trun with no prompt [False]')
    print('\t\t-useHeader\tGet size from tiepoints dir header [False]')
    print('\t\t-SVAlongTrack\tcreate velocity with azimuth varying '
          'corrections to SV baseline/offsets (if available)')
    print('\t\t-SVconst\tcreate velocity_with range and azimuth constant '
          'corrections to SV derived results (if available)')
    print('\t\t-firstdate\trun only those with velocity/mosaicOffsets.vx'
          ' >= firstdate [1990:12:31]')
    print('\t\t-lastdate\trun only those with velocity/mosaicOffsets.vx <= '
          'lastdate [2100:12:31]')
    print('\t\t-redoculled\tforces a rerun for all -- culled -- velocity '
          'directories')
    print('\t\t\t\t(note: -redoculled rebuilds velocity/ ONLY; run again '
          'without it, with -reset, to rebuild velocity_nocull/)')
    print('\t\t-tiff\t\tforce GeoTIFF output (mosaic3d -GTiff) regardless of '
          'project.yaml velThumbOutput')
    print('\t\t-reSize\t\tforces files to be run with natural bounding '
          'box new makevelstatsregions.py can run')
    print('\t\t-resolution\tOverrides resolution with new string [None]')
    print('\t\t-toRun\t\tPython formated list of tracks directories to '
          'process [track-26,74,90,112,141,170]\n')
    print('\n\tNote:\n\t\t 1) Run in main track directory and it will '
          'operate on orbit_frame')
    print('\t\t 2) Requires an  existing orbit_frame/velocity dir '
          '(makeframe.py)')
    print('\t\t 3) SVAlongTrack uses az.est.svlinear and rBaseline.quad '
          'if available. Reverts to azest and rBaseline if not')
    print('\t\t 4) SVConst uses az.est.const and rBaseline.deltaBp if '
          'available. Reverts to azest and rBaseline if not')
    print('\t\t \033[0m\n')
    print('Part of the mosaicworkflow package.')
    exit()


def processDate(arg):
    try:
        firstDate = datetime.strptime(arg.split('=')[-1], "%Y:%m:%d")
    except Exception:
        u.myerror(f'invalid firstdate {arg}')
    return firstDate


def getMakeVelNocleanArgs():
    args = sys.argv[1:]
    reset = False
    redoCulled = False
    maxThreads = 12
    usePrompt = True
    flavor = None
    resetSize = False
    resolution = None
    useHeader = False
    forceTiff = False
    # default firstdate ensures all are processed
    firstDate = datetime(1990, 12, 31)
    lastDate = datetime(2100, 12, 31)
    # Default: every track-* directory here, numerically sorted (track-2 before
    # track-10). Was a hardwired six-track list from an old project, which
    # silently processed the wrong subset in any other tree.
    toRun = sorted([d for d in glob.glob('track-*') if os.path.isdir(d)],
                   key=lambda p: int(re.search(r'track-(\d+)', p).group(1)))
    if 'track' in os.getcwd():
        toRun = ['.']
    if len(sys.argv) > 1:
        for arg in args:
            if '-help' in arg:
                makevelnocleanUsage()
            elif '-reset' in arg:
                reset = True
            elif '-redoculled' in arg:
                redoCulled = True
            elif '-SVAlongTrack' in arg:
                flavor = 'SVAlongTrack'
            elif '-SVConst' in arg:
                flavor = 'SVConst'
            elif '-toRun' in arg:
                print(arg)
                toRun = eval(arg.split('=')[-1])
            elif '-firstdate' in arg:
                firstDate = processDate(arg)
            elif '-lastdate' in arg:
                lastDate = processDate(arg)
            elif '-noprompt' in arg:
                usePrompt = False
            elif '-useHeader' in arg:
                useHeader = True
            elif '-tiff' in arg:
                forceTiff = True
            elif '-reSize' in arg:
                resetSize = True
            elif 'resolution' in arg:
                resolution = arg.split('=')[-1]
            elif '-threads' in arg:
                maxThreads = int(arg.split('=')[-1])
                if maxThreads > 30:
                    maxThreads = 30
            else:
                u.mywarning('Invalid argument '+arg)
                makevelnocleanUsage()
    if len(toRun) < 1:
        u.myerror('no track-* directories found -- run from the project root '
                  '(above the tracks), from inside a track dir, or pass -toRun')
    print('reset = \t\t', reset)
    print('redoCulled = \t', redoCulled)
    print('toRun = \t\t', toRun)
    print('threads = \t\th', maxThreads)
    return reset, redoCulled, toRun, maxThreads, \
        usePrompt, firstDate, lastDate, flavor, resetSize, resolution, \
        useHeader, forceTiff


def getMetaDate(myV):
    metaFile = f'{myV}/mosaicOffsets.meta'
    try:
        fp = open(metaFile, 'r')
    except Exception:
        u.mywarning(f'could not open {metaFile}')
        return None
    #
    for line in fp:
        if "First Image Date" in line:
            fp.close()
            return datetime.strptime(line.split('=')[1].strip(), "%b:%d:%Y")
    u.mywarning(f'problem parsing date {metaFile}')
    return None


def velDirDate(myV, vx_suffixes):
    ''' (velDir, first image date) when the dir holds a velocity product, else
    None. Reads only, so these run concurrently. '''
    vxFile = next((os.path.join(myV, s) for s in vx_suffixes
                   if os.path.exists(os.path.join(myV, s))), None)
    if vxFile is None:
        return None
    return myV, getMetaDate(myV)


def getVelDirs(toRun, firstDate, lastDate, tiff_output=False, maxThreads=12):
    ''' get orginal velocity and screen by date directories.

    The per-track globs and the per-frame metadata reads are both NFS-bound
    (thousands of round trips), so each is run through a thread pool. Results
    are reassembled in the original track order, and glob order within a track,
    so the run list stays deterministic. '''
    veldirs = []
    myDates = []
    print(f'Only files modified after {firstDate} will run')
    # With tiff output selected, a directory still holding only the binary
    # product is a conversion candidate, not something to skip -- accept either
    # spelling so -redoculled/-reset can rewrite it as GeoTIFF.
    vx_suffixes = ['mosaicOffsets.vx.tif', 'mosaicOffsets.vx'] if tiff_output \
        else ['mosaicOffsets.vx']
    with ThreadPoolExecutor(max_workers=max(min(maxThreads, len(toRun)), 1)) \
            as pool:
        perTrack = list(pool.map(
            lambda trackDir: sorted(glob.glob(trackDir + '/*/velocity')),
            toRun))
    candidates = [myV for trackVels in perTrack for myV in trackVels]
    with ThreadPoolExecutor(max_workers=max(maxThreads, 1)) as pool:
        found = list(pool.map(lambda myV: velDirDate(myV, vx_suffixes),
                              candidates))
    # filter by date
    for entry in found:
        if entry is None:
            continue
        myV, myDate = entry
        # getMetaDate warns and returns None when the meta file will not parse
        if myDate is not None and firstDate <= myDate <= lastDate:
            veldirs.append(myV)
            myDates.append(myDate)
    if len(veldirs) < 1:
        u.myerror('no velocity dirs: In directory above tracks ? '
                  'Correct tracks specified ?')
    #
    return veldirs, myDates


def addFlagToRun(destrun, flavor, tiff_output=False):
    flags = {'SVConst': ' -SVConst ', ' -SVAlongTrack ': ' -SVAlongTrack ',
             'None': ' -SVConst '}  # Force all to SVConst
    fpIn = open(destrun, 'r')
    lines = []
    for line in fpIn:
        lines.append(line)
    fpIn.close()
    fpOut = open(destrun, 'w')
    for line in lines:
        if 'mosaic3d' in line:
            if '-vzFlag' not in line:
                line = line.replace('mosaic3d ', 'mosaic3d -vzFlag 3 ')
            # Output format follows project.yaml velThumbOutput, so a runOff
            # inherited from an earlier run in the other format is rewritten
            # here rather than silently reproducing that format. Kept
            # symmetric: dropping -GTiff matters as much as adding it, since
            # the existence checks above look for the format the project asks
            # for and would otherwise rebuild the same dir on every pass.
            if tiff_output and '-GTiff' not in line:
                line = line.replace('mosaic3d ', 'mosaic3d -GTiff ')
            elif not tiff_output and '-GTiff' in line:
                line = line.replace(' -GTiff', '')
            line = line.replace('mosaic3d ', 'mosaic3d ' + flags[flavor])
        fpOut.write(line)
    fpOut.close()


def resetInputFileSize(velDir, resolution=None):
    ''' Set size params for inputFile to 0 0 to force natural size'''
    fpIn = open(f'{velDir}/inputFile', 'r')
    fpOut = open(f'{velDir}/inputFile.tmp', 'w')
    #
    notDone = True
    for line in fpIn:
        if ';' not in line and notDone:
            pieces = line.split()
            if len(pieces) != 6:  # swap size for 0 0
                u.myerror(f'problem with inputfile: {velDir}/inputFile')
            if resolution is not None:
                line = f'{resolution}\n'
            else:
                pieces[2:4] = ['0.0', '0.0']
                line = ' '.join(pieces) + '\n'
        notDone = False
        print(line, file=fpOut, end='')
    fpOut.close()
    fpIn.close()
    #
    if os.path.exists(f'{velDir}/inputFile.old'):  # Remove old if present
        os.remove(f'{velDir}/inputFile.old')
    os.rename(f'{velDir}/inputFile', f'{velDir}/inputFile.old')
    os.rename(f'{velDir}/inputFile.tmp', f'{velDir}/inputFile')


def createNoCullInput(inputFile, velDirNoCull, mainDir, fastOffs):
    ''' Make an input files for no cull case with both slow and fast offsets'''
    fin = open(inputFile, 'r')
    fout = open(velDirNoCull+'/inputFile', 'w')
    #
    lines, nLines, nOrig = [], 0, 0
    # for each line copy, if azimuth in the line, duplicate
    # with noclean.fast. versions
    # build a list of lines
    for line in fin:
        if 'azimuth' in line:
            line1 = subOffsets(line, mainDir)
            if line1 is not None:
                lines.append(line1)
                nLines += 1
                nOrig += 1
                if os.path.exists(fastOffs):
                    line2 = subOffsets(line, mainDir, fast=True)
                    lines.append(line2)
                    nLines += 1
        else:
            lines.append(line)
    # print lines, & update count if fast lines have been ins
    for line in lines:
        # rough test to find nLines in the input file
        if len(line.strip()) < 5 and ';' not in line:
            try:
                nL = int(line)
                if nL == nOrig:
                    line = str(nLines) + '\n'
            except Exception:
                pass
        # insert semicolon if fast line added
        if 'fast' in line:
            print(';', file=fout)
        print(line, file=fout, end='')
    fin.close()
    fout.close()


def setFlavor(defaultFlavor, velDate, redoCulled):
    ''' For noculled cases after Dec 31, 2019, set flavor to deltabp
        Update: April 28, 2025, changing to do for all years'''
    # return old cases, or redoculled cases
    #if redoCulled or velDate < datetime(2000, 1, 1):
    #    return defaultFlavor
    # This will force all non-culled cases after Jan-1-2020 to use 'SVConst'
    return 'SVConst'


def getResolution(header):
    ''' get resolution and frame range from vel_thumb_header. Returns {} for a
    file that is not a <first>dash<last> header, or for a header with no
    resolution line yet (makevelstatsregions has not been run on it)'''
    resDict = {}
    m = headerPattern.match(os.path.basename(header))
    if m is None:
        return resDict
    frame1, frame2 = int(m.group(1)), int(m.group(2))
    res = None
    with open(header, 'r') as fp:
        for line in fp:
            if 'resolution' in line:
                res = line.split('=')[-1].replace('\"', '').replace("\'", '')
                res = res.strip()
                break
    if res is None:
        u.mywarning(f'no resolution line in {header} - skipping')
        return resDict
    for frame in range(frame1, frame2+1):
        resDict[frame] = res
    return resDict


def buildResolutionDictionary(toRun):
    ''' for list of tracks toRun, create a dictionary by frame that includes
    the reoslution string from the vel_thumb_header files'''
    print('Reading Headers from vel_thumb_headers:')
    resolutionDict = {}  # Master dict
    for track in toRun:
        resolutionDict[track] = {}  # track specific dict
        headers = sorted(glob.glob(f'{track}/tiepoints/vel_thumb_header_*'))
        print(track)
        for header in headers:
            if not headerPattern.match(os.path.basename(header)):
                continue
            print(header)
            resolutionDict[track].update(getResolution(header))
    return resolutionDict


def main():
    reset, redoCulled, toRun, maxThreads, usePrompt, firstDate, lastDate,\
        defaultFlavor, resetSize, resolution, useHeader, forceTiff = \
        getMakeVelNocleanArgs()
    nNotRun = 0

    tiff_output = False
    if os.path.exists('project.yaml'):
        with open('project.yaml') as _yf:
            _proj = yaml.safe_load(_yf)
        if isinstance(_proj, dict) and _proj.get('velThumbOutput') == 'tiff':
            tiff_output = True
            print('velThumbOutput: tiff (from project.yaml)')
    # -tiff forces GeoTIFF whether or not project.yaml says so, matching the
    # --tiff that setupS1Tracks passes to refreshties/vel_thumbs. Without it a
    # project missing velThumbOutput would get a GeoTIFF velocity/ and a binary
    # velocity_nocull/, since the nocull runOff is copied from velocity/ and
    # addFlagToRun would then strip the inherited -GTiff.
    if forceTiff:
        tiff_output = True
        print('velThumbOutput: tiff (from -tiff)')

    # get date filtered list of files (default is all)
    veldirs, velDates = getVelDirs(toRun, firstDate, lastDate, tiff_output,
                                   maxThreads=maxThreads)
    #for veldir, velDate in zip(veldirs, velDates):
    #    if velDate > datetime(2022, 1, 1):
    #        print(veldir, velDate)
    if useHeader:
        resolutionDictionary = buildResolutionDictionary(toRun)
    #
    # Loop through velocity directories
    velToRun = []

    def setupRedoCulled(veldir):
        ''' -redoculled work for one directory: refresh the runOff flags and,
        where a thumbnail resolution applies, the input file size. Per-frame
        file I/O on distinct paths with no shared state, so it parallelises.
        Returns the dir to run, or None when the frame is excluded. '''
        # An excluded frame is not rerun; its products are removed so a stale
        # velocity cannot reach the mosaic. Logged, not warned about.
        if frameExcluded(veldir):
            logMessage(f'{veldir}: Exclude present - not run, removed '
                       f'{removeProducts(veldir)} product files')
            return None
        # this will add vzFlag 3 to runOff for Knut's work
        addFlagToRun(f'{veldir}/runOff', 'None', tiff_output=tiff_output)
        myResolution = resolution
        if useHeader:
            frame = int(veldir.split('/')[1].split('_')[-1])
            track = veldir.split('/')[0]
            myResolution = resolutionDictionary[track][frame]
        if resetSize or myResolution is not None:
            resetInputFileSize(veldir, resolution=myResolution)
        return veldir  # just add prior veldir to run list

    if redoCulled:
        # Thousands of NFS round trips with nothing shared between frames;
        # pool.map keeps velToRun in the original order.
        with ThreadPoolExecutor(max_workers=max(maxThreads, 1)) as pool:
            velToRun = [v for v in pool.map(setupRedoCulled, veldirs)
                        if v is not None]
    # Serial: createNoCullInput -> subOffsets uses pushd/popd, which chdirs the
    # whole process and so cannot be threaded.
    for veldir, velDate in zip(veldirs, velDates):
        if redoCulled:
            break
        inputFile = veldir + '/inputFile'
        # setup no cull dir
        mainDir = "/".join(veldir.split('/')[0:2])
        if defaultFlavor is None:
            velDirNoCull = f'{veldir}_nocull'
        else:
            velDirNoCull = f'{veldir}_{defaultFlavor}'
        #
        # Excluded frames are not built; drop anything a prior run left behind.
        if frameExcluded(veldir):
            logMessage(f'{velDirNoCull}: Exclude present - not run, removed '
                       f'{removeProducts(velDirNoCull)} product files')
            continue
        nocull_vx = velDirNoCull + ('/mosaicOffsets.vx.tif' if tiff_output else '/mosaicOffsets.vx')
        if not reset and os.path.exists(nocull_vx):
            nNotRun += 1
            continue
        # mkcull directory if needed
        if not os.path.isdir(velDirNoCull):
            os.mkdir(velDirNoCull)
        fastOffs = f'{mainDir}/fast/azimuth.offsets.noclean.fast'
        # if velocity has valid input file clone input file
        if os.path.exists(inputFile):
            createNoCullInput(inputFile, velDirNoCull, mainDir, fastOffs)
            #
            flavor = setFlavor(defaultFlavor, velDate, redoCulled)
            #
            srcrun = veldir+'/runOff'
            destrun = velDirNoCull+'/runOff'
            if os.path.exists(srcrun):
                copyfile(srcrun, destrun)
                velToRun.append(velDirNoCull)
            if flavor is not None:
                addFlagToRun(destrun, flavor, tiff_output=tiff_output)
            else:
                addFlagToRun(destrun, 'None', tiff_output=tiff_output)
    u.mywarning(f'{nNotRun} products already exist so skipping (set '
                'reset flag to rebuild)')
    #
    # now set up threads
    threads = []
    for runfile in velToRun:
        # print(runfile)
        thread = threading.Thread(target=runRunFile, args=[runfile])
        threads.append(thread)
    # run threads
    u.runMyThreads(threads, maxThreads, 'Make Velocity', prompt=usePrompt)


if __name__ == '__main__':
    main()
