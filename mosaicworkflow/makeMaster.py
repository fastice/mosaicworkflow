#!/usr/bin/env python3
"""
Assemble master inputFile(s) for mosaic3d from several source regions.

Two modes:

  config mode   makeMaster [masterInput.yaml]
      A YAML file names the sources (project roots, flat tiepoint directories,
      or already-assembled master files) and the named outputs built from them.
      Any number of outputs can be defined over any subset of sources without
      touching this code.

  legacy mode   makeMaster --projectDir DIR [--inputPaths PATH ...] [--outputPath PATH]
      The original NISAR behaviour, kept verbatim so existing callers
      (nisargrimpworkflow.autoupdate) keep working unchanged.  Equivalent to a
      config with one trackTree source per --inputPaths entry.

Data lines are verified before they are written: every referenced file must
exist (accepting a .vrt sidecar for a bare name, as getROffsets.c does), any
frame carrying an Exclude marker is dropped, and -- where enabled -- each range
baseline's recorded offsetCorrectionFile must agree with the range.offsets
VRT's ionosphereRangeOffsetCorrection metadata, the consistency
mosaic3d/checkForIonosphereCorrection() enforces at runtime.  Entries failing
any check are skipped and summarised at the end; --verbose lists them.

Replaces both nisargrimpworkflow/makeMaster.py and the Sentinel-1
Release/masterInput/createMaster.ipynb notebook.
"""

import argparse
import os
import re
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from functools import lru_cache

import utilities as u
import yaml
from osgeo import gdal

#
# The nine fields of a mosaic3d input data line, in order.
#
FIELDKEYS = ['phase', 'geodat', 'baseline', 'nDays', 'weight',
             'azOff', 'azParams', 'rgOff', 'rParams']

#
# Per-source options and their defaults.  Anything here may be set in the
# config's top-level `defaults:` block or overridden on an individual source.
#
SOURCEDEFAULTS = {
    'layout': None,             # trackTree | flat | verbatim
    'root': None,               # project root (trackTree) or tiepoints dir (flat)
    'file': None,               # verbatim: one already-assembled master
    'files': None,              # verbatim: several
    'filePattern': None,        # {track} / {year} substituted, globs allowed
    'years': None,              # flat: {first: , last: } or an explicit list
    'tracks': None,             # trackTree: optional whitelist
    'minFields': 9,             # drop lines with fewer fields than this
    'lineFilter': 'fields',     # fields | heuristic
    'normalizeSlashes': True,   # '//' -> '/' and re-join fields on one space
    'excludeMarker': 'Exclude',  # marker filename beside the geodat; None disables
    'excludeLevels': 1,         # directory levels above the geodat to test
    'checkFields': ['geodat', 'azOff', 'azParams', 'rgOff', 'rParams'],
    'nonFileTokens': ['none', 'nophase', 'nil'],
    'vrtFallback': True,        # accept <tok>.vrt when <tok> is absent
    'checkIonosphere': False,   # baseline <-> range.offsets VRT consistency
    'groupComment': None,       # None | track | file
    'frameStatPattern': None,   # e.g. 'track-*/*_0???/motion'
}

#
# Per-output options and their defaults.
#
OUTPUTDEFAULTS = {
    'path': None,               # directory the inputFile is written into
    'sources': None,            # source or group names, in emission order
    'fileName': 'inputFile',
    'headerComments': True,     # emit '; Output area' / '; Number of inputs'
    'ordering': 'bySource',     # bySource | byTrack
    'outputArea': None,         # 'xo yo xs ys deltaX deltaY', or a list
}

#
# Legacy (NISAR) mode: the settings that reproduce nisargrimpworkflow's
# makeMaster, except that 'nophase' is now correctly treated as a placeholder.
# The old tool read it as a file path and so dropped every offsets-only
# product -- defeating the ExcludePhase.pending demotion, whose whole point is
# "phase dropped, offsets kept".  See Documents/makeMaster.md.
#
LEGACYSOURCE = {
    'layout': 'trackTree',
    'filePattern': 'track-{track}/tiepoints/track-{track}/inputFile-{track}-all20??',
    'lineFilter': 'heuristic',
    'normalizeSlashes': False,
    'excludeMarker': None,
    'checkFields': 'all',
    'vrtFallback': True,
    'checkIonosphere': True,
    'groupComment': 'track',
    'frameStatPattern': 'track-*/*_0???/motion',
}

LEGACYOUTPUT = {
    'headerComments': False,
    'ordering': 'byTrack',
    'outputArea': '0 0 0 0 .2 .2',
}

RED = '\033[1;31m'
REDBG = '\033[1;41m'
YELLOWBG = '\033[1;43m'
RESET = '\033[0m'


def parseArgs():
    '''Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mAssemble master inputFile(s) for mosaic3d '
        'from several source regions\033[0m\n\n',
        epilog='Part of the mosaicworkflow package.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('config', type=str, nargs='?', default=None,
                        help='YAML config describing sources and outputs '
                        '[masterInput.yaml in the current directory]')
    parser.add_argument('--outputs', nargs='+', metavar='NAME', default=None,
                        help='Build only these named outputs [all of them]')
    parser.add_argument('--sources', nargs='+', metavar='NAME', default=None,
                        help='Restrict every output to these named sources '
                        '[all the sources it lists]')
    parser.add_argument('--workers', type=int, default=8,
                        help='Threads used for the file-existence checks [8]')
    parser.add_argument('--check', action='store_true',
                        help='Run every check and print the summary, but do '
                        'not write any output file')
    parser.add_argument('--verbose', action='store_true',
                        help='List every missing file path in the skip summary')
    # Legacy (NISAR) flags -- defaults are None so an explicit value is
    # distinguishable from an omitted one.
    parser.add_argument('--projectDir', type=str, default=None,
                        help='Legacy mode: root for anchoring a relative '
                        '--outputPath [.]. Ignored for track discovery when '
                        '--inputPaths is given.')
    parser.add_argument('--inputPaths', nargs='+', metavar='PATH', default=None,
                        help='Legacy mode: one or more directories to search '
                        'for track-N subdirectories (overrides --projectDir '
                        'for input). Example: . ../anotherDir')
    parser.add_argument('--outputPath', default=None, metavar='PATH',
                        help='Legacy mode: directory to write inputFile '
                        '[Release/masterInput, relative to --projectDir '
                        'unless absolute]')
    return parser.parse_args()


# --------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------

def readConfig(configPath):
    '''Load and validate the YAML config, anchoring every path to the config
    file's own directory.  Modelled on reprocessS1.loadConfig.'''
    if not os.path.exists(configPath):
        u.myerror(f'{configPath}: config file not found')
    with open(configPath) as fp:
        config = yaml.safe_load(fp) or {}
    for key in ['sources', 'outputs']:
        if key not in config:
            sys.exit(f'{configPath}: missing required key {key}')
    config['configPath'] = os.path.abspath(configPath)
    config['configDir'] = os.path.dirname(config['configPath'])
    return config


def anchorPath(path, configDir):
    '''Absolute paths are used as given; relative ones hang off the config.'''
    if os.path.isabs(path):
        return os.path.normpath(path)
    return os.path.normpath(os.path.join(configDir, path))


def mergeDefaults(entry, defaults, allowed):
    '''Merge an entry over the config defaults over the built-in defaults,
    keeping only the keys this kind of entry understands.'''
    merged = dict(allowed)
    for key, value in defaults.items():
        if key in merged:
            merged[key] = value
    for key, value in entry.items():
        if key not in merged and key != 'name':
            u.mywarning(f'{entry.get("name", "?")}: unknown key {key}, ignored')
        merged[key] = value
    return merged


def resolveSources(config):
    '''Return {name: source} for every usable source, in config order.'''
    defaults = config.get('defaults') or {}
    configDir = config['configDir']
    sources = {}
    for index, entry in enumerate(config['sources']):
        name = entry.get('name')
        if not name:
            sys.exit(f'{config["configPath"]}: every source needs a name')
        if name in sources:
            sys.exit(f'{config["configPath"]}: duplicate source name {name}')
        source = mergeDefaults(entry, defaults, SOURCEDEFAULTS)
        source['name'] = name
        source['index'] = index
        if source['layout'] not in ['trackTree', 'flat', 'verbatim']:
            sys.exit(f'{name}: layout must be trackTree, flat or verbatim')
        if source['layout'] in ['trackTree', 'flat']:
            if not source['root']:
                sys.exit(f'{name}: layout {source["layout"]} needs a root')
            source['root'] = anchorPath(source['root'], configDir)
            if not os.path.isdir(source['root']):
                u.mywarning(f'source root not found, skipping: '
                            f'{source["root"]}')
                continue
        else:
            files = source['files'] or ([source['file']] if source['file']
                                        else None)
            if not files:
                sys.exit(f'{name}: layout verbatim needs file or files')
            source['files'] = [anchorPath(f, configDir) for f in files]
        sources[name] = source
    return sources


def resolveOutputs(config, sources):
    '''Expand groups, validate source names, and return the output specs.'''
    defaults = config.get('defaults') or {}
    configDir = config['configDir']
    groups = config.get('groups') or {}
    outputs = []
    for entry in config['outputs']:
        name = entry.get('name')
        if not name:
            sys.exit(f'{config["configPath"]}: every output needs a name')
        output = mergeDefaults(entry, defaults, OUTPUTDEFAULTS)
        output['name'] = name
        if not entry.get('path'):
            sys.exit(f'{name}: output needs a path')
        output['path'] = anchorPath(entry['path'], configDir)
        if output['outputArea'] is None:
            sys.exit(f'{name}: no outputArea, and none in defaults')
        if output['ordering'] not in ['bySource', 'byTrack']:
            sys.exit(f'{name}: ordering must be bySource or byTrack')
        # Expand group names, drop duplicates, keep order.
        wanted = []
        for item in (entry.get('sources') or []):
            for sourceName in (groups[item] if item in groups else [item]):
                if sourceName not in wanted:
                    wanted.append(sourceName)
        if not wanted:
            sys.exit(f'{name}: output lists no sources')
        missing = [s for s in wanted if s not in sources]
        if missing:
            u.mywarning(f'{name}: skipping unavailable source(s) '
                        f'{" ".join(missing)}')
        output['sources'] = [s for s in wanted if s in sources]
        outputs.append(output)
    return outputs


def legacyConfig(args):
    '''Build the in-memory config that reproduces the original NISAR tool.'''
    projectDir = args.projectDir if args.projectDir is not None else '.'
    outputPath = (args.outputPath if args.outputPath is not None
                  else 'Release/masterInput')
    inputPaths = args.inputPaths if args.inputPaths else [projectDir]
    if not os.path.isabs(outputPath):
        outputPath = os.path.join(projectDir, outputPath)
    sources = {}
    for index, root in enumerate(inputPaths):
        source = dict(SOURCEDEFAULTS)
        source.update(LEGACYSOURCE)
        source['name'] = root
        source['index'] = index
        source['root'] = root
        sources[root] = source
    output = dict(OUTPUTDEFAULTS)
    output.update(LEGACYOUTPUT)
    output['name'] = 'legacy'
    output['path'] = outputPath
    output['sources'] = list(sources.keys())
    return sources, [output]


# --------------------------------------------------------------------------
# Discovery -- every layout returns a list of groups
#   {'label': str, 'trackNum': int|None, 'files': [path, ...]}
# --------------------------------------------------------------------------

def discoverTrackTree(source, stats):
    '''root/track-N/tiepoints/track-N/<filePattern>, one group per track.'''
    wanted = None
    if source['tracks']:
        wanted = set(int(str(t).replace('track-', '')) for t in source['tracks'])
    trackNums = []
    for trackDir in u.globSort(os.path.join(source['root'], 'track-*')):
        match = re.fullmatch(r'track-(\d+)', os.path.basename(trackDir))
        if match:
            trackNums.append(int(match.group(1)))
    groups = []
    for trackNum in trackNums:
        if wanted is not None and trackNum not in wanted:
            continue
        stats['nTrackEntries'] += 1
        pattern = source['filePattern'].format(track=trackNum)
        files = u.globSort(os.path.join(source['root'], pattern))
        if not files:
            stats['tracksNoInputFile'].append(trackNum)
            continue
        groups.append({'label': f'track-{trackNum}', 'trackNum': trackNum,
                       'files': files})
    return groups


def discoverFlat(source, stats):
    '''root/<filePattern with {year}>, one group per year file.'''
    years = source['years'] or {'first': 2014, 'last': 2029}
    if isinstance(years, dict):
        years = range(int(years['first']), int(years['last']) + 1)
    thisYear = datetime.now().year
    groups = []
    for year in years:
        path = os.path.join(source['root'],
                            source['filePattern'].format(year=year))
        if not os.path.exists(path):
            # Years still in the future are simply not there yet.
            if int(year) <= thisYear:
                u.mywarning(f'{source["name"]}: no {path}')
            continue
        stats['fileTimes'].append(
            (path, datetime.fromtimestamp(os.path.getmtime(path))))
        groups.append({'label': os.path.basename(path), 'trackNum': None,
                       'files': [path]})
    return groups


def discoverVerbatim(source, stats):
    '''Already-assembled master file(s), spliced in as-is.'''
    groups = []
    for path in source['files']:
        for match in (u.globSort(path) if '*' in path or '?' in path
                      else [path]):
            if not os.path.exists(match):
                u.mywarning(f'{source["name"]}: no {match}')
                continue
            stats['fileTimes'].append(
                (match, datetime.fromtimestamp(os.path.getmtime(match))))
            groups.append({'label': os.path.basename(os.path.dirname(match)),
                           'trackNum': None, 'files': [match]})
    return groups


def discoverFiles(source, stats):
    '''Dispatch on layout.'''
    if source['layout'] == 'trackTree':
        return discoverTrackTree(source, stats)
    if source['layout'] == 'flat':
        return discoverFlat(source, stats)
    return discoverVerbatim(source, stats)


# --------------------------------------------------------------------------
# Parsing and checking
# --------------------------------------------------------------------------

def isDataLine(line):
    '''Heuristic filter: skip blanks, ; comments and pure-number lines.'''
    stripped = line.strip()
    if not stripped:
        return False
    if stripped.startswith(';'):
        return False
    try:
        [float(tok) for tok in stripped.split()]
        return False
    except ValueError:
        pass
    return True


def readDataLines(filePath, source):
    '''Return (rawLine, fields) for every data line in one inputFile.'''
    lines = []
    with open(filePath) as fp:
        for line in fp:
            if source['lineFilter'] == 'heuristic':
                if not isDataLine(line):
                    continue
                raw = line.rstrip()
                fields = raw.split()
            else:
                # The notebook's filter: any ';' anywhere disqualifies a line.
                if ';' in line:
                    continue
                raw = line.rstrip()
                if source['normalizeSlashes']:
                    raw = raw.replace('//', '/')
                fields = raw.split()
                if len(fields) < source['minFields']:
                    continue
            lines.append((raw, fields))
    return lines


def renderLine(raw, fields, source):
    '''The text actually written out for a kept line.'''
    if source['normalizeSlashes']:
        return ' '.join(fields)
    return raw


def isFileToken(tok, nonFileTokens):
    '''True for tokens that are file paths, not numbers or placeholders.'''
    if tok.lower() in nonFileTokens:
        return False
    try:
        float(tok)
        return False
    except ValueError:
        return True


def tokenExists(tok, vrtFallback):
    '''True if tok exists as given, or as tok + '.vrt' -- inputFile templates
    write bare names like 'range.offsets'/'azimuth.offsets' unconditionally,
    but projects that have moved to VRT-only storage rely on rparams/azparams
    (getROffsets.c) falling back to the .vrt sidecar for the bare name.'''
    if os.path.exists(tok):
        return True
    return vrtFallback and os.path.exists(tok + '.vrt')


@lru_cache(maxsize=None)
def findIonCorrection(rangeOffsetsPath):
    '''If rangeOffsetsPath's .vrt carries an ionosphereRangeOffsetCorrection
    tag (written by SetupNISAR.globalFillIonosphere), return the resolved path
    to that correction file -- relative to rangeOffsetsPath's own directory,
    matching how checkForIonosphereCorrection() resolves it.  None when there
    is no .vrt or no such tag (nothing to verify, not an error).'''
    vrtPath = rangeOffsetsPath + '.vrt'
    if not os.path.exists(vrtPath):
        return None
    ds = gdal.Open(vrtPath)
    if ds is None:
        return None
    ionName = ds.GetMetadataItem('ionosphereRangeOffsetCorrection')
    ds = None
    if not ionName:
        return None
    return os.path.join(os.path.dirname(rangeOffsetsPath), ionName)


@lru_cache(maxsize=None)
def baselineOffsetCorrection(baselinePath):
    '''Basename of the range baseline's offsetCorrectionFile entry, or None if
    it records no correction.  'nil' and a missing/blank entry both map to
    None, matching getRParams()/readOffsets.c.'''
    try:
        with open(baselinePath) as fp:
            for line in fp:
                if line.strip().startswith('offsetCorrectionFile:'):
                    val = line.split(':', 1)[1].strip()
                    if not val or val == 'nil':
                        return None
                    return os.path.basename(val)
    except OSError:
        return None
    return None


def checkLineIonConsistency(fields, source):
    '''Mismatch messages when a line's range baseline records an
    offsetCorrectionFile whose basename disagrees with the range.offsets VRT's
    ionosphereRangeOffsetCorrection metadata -- the same inconsistency
    mosaic3d/checkForIonosphereCorrection() aborts on.  Empty when consistent
    or not applicable.'''
    rangeOffsetsTok = None
    baselineTok = None
    for tok in fields:
        if not isFileToken(tok, source['nonFileTokens']):
            continue
        base = os.path.basename(tok)
        if base == 'range.offsets':
            rangeOffsetsTok = tok
        elif base.startswith('rBaseline') and base.endswith('.yaml'):
            baselineTok = tok
    if rangeOffsetsTok is None or baselineTok is None:
        return []
    blName = baselineOffsetCorrection(baselineTok)
    if blName is None:
        return []
    ionPath = findIonCorrection(rangeOffsetsTok)
    vrtName = os.path.basename(ionPath) if ionPath is not None else None
    if blName != vrtName:
        return [f'{baselineTok}: offsetCorrectionFile {blName} != VRT '
                f'ionosphereRangeOffsetCorrection {vrtName}']
    return []


def tokensToCheck(fields, source):
    '''The tokens whose existence this source wants verified.'''
    checkFields = source['checkFields']
    if checkFields == 'all':
        candidates = fields
    elif not checkFields:
        return []
    else:
        candidates = []
        for key in checkFields:
            index = FIELDKEYS.index(key)
            if index < len(fields):
                candidates.append(fields[index])
    return [tok for tok in candidates
            if isFileToken(tok, source['nonFileTokens'])]


def checkLineFiles(fields, source):
    '''Missing file paths for a single data line.'''
    missing = []
    for tok in tokensToCheck(fields, source):
        if not tokenExists(tok, source['vrtFallback']):
            missing.append(tok)
        elif os.path.basename(tok) == 'range.offsets' and \
                source['checkIonosphere']:
            ionPath = findIonCorrection(tok)
            if ionPath is not None and not os.path.exists(ionPath):
                missing.append(ionPath)
    return missing


def isExcluded(fields, source):
    '''True if the frame carries the source's Exclude marker.  The marker sits
    beside the geodat (excludeLevels 1); mosaic3d itself never sees it.'''
    marker = source['excludeMarker']
    if not marker:
        return False
    index = FIELDKEYS.index('geodat')
    if index >= len(fields):
        return False
    path = fields[index]
    for _ in range(int(source['excludeLevels'])):
        path = os.path.dirname(path)
    return os.path.exists(os.path.join(path, marker))


def fileTypeLabel(path):
    '''Categorise a missing file path into a short human-readable label.'''
    base = os.path.basename(path)
    if 'rBaseline' in base:
        return 'rBaseline'
    if 'range.offsets' in base:
        return 'range.offsets'
    if 'azimuth.offsets' in base:
        return 'azimuth.offsets'
    if 'ionosphere' in base.lower():
        return 'ionosphere correction'
    if base.endswith('.geojson') or 'geodat' in base:
        return 'geodat'
    return base


def orbitFrameFromLine(fields, source):
    '''Extract the orbit_frame token (e.g. '4040_0000') from a data line.'''
    for tok in fields:
        if not isFileToken(tok, source['nonFileTokens']):
            continue
        for part in tok.replace('\\', '/').split('/'):
            if re.fullmatch(r'\d+_\d+', part):
                return part
    return None


# --------------------------------------------------------------------------
# Collection
# --------------------------------------------------------------------------

def newStats():
    '''The per-source counters carried through collection and reporting.'''
    return {'nTotal': 0, 'nSkipped': 0, 'nExcluded': 0, 'nMismatch': 0,
            'nTrackEntries': 0, 'skipCounts': Counter(),
            'skipFiles': defaultdict(list), 'mismatchMsgs': [],
            'tracksNoInputFile': [], 'excludedDirs': [], 'fileTimes': [],
            'framesInOutput': set(), 'nKept': 0}


def collectSource(source, workers, verbose):
    '''Read, check and keep every data line of one source.

    Returns (groups, stats) where each group gains a 'blocks' entry: one list
    of rendered output lines per input file.  Called once per source however
    many outputs use it.
    '''
    stats = newStats()
    groups = discoverFiles(source, stats)
    # GDAL is kept on one thread; pure stat() work is worth spreading over NFS.
    threaded = workers > 1 and not source['checkIonosphere']
    executor = ThreadPoolExecutor(max_workers=workers) if threaded else None
    for group in groups:
        group['blocks'] = []
        for filePath in group['files']:
            lines = readDataLines(filePath, source)
            stats['nTotal'] += len(lines)
            kept = []
            for raw, fields in lines:
                if isExcluded(fields, source):
                    stats['nExcluded'] += 1
                    if verbose:
                        index = FIELDKEYS.index('geodat')
                        stats['excludedDirs'].append(
                            os.path.dirname(fields[index]))
                    continue
                kept.append((raw, fields))
            if executor is not None and kept:
                # executor.map preserves input order, so the result is
                # identical to the serial path.
                missingList = list(executor.map(
                    lambda item: checkLineFiles(item[1], source), kept))
            else:
                missingList = [checkLineFiles(f, source) for _, f in kept]
            block = []
            for (raw, fields), missing in zip(kept, missingList):
                mismatches = (checkLineIonConsistency(fields, source)
                              if source['checkIonosphere'] else [])
                if missing or mismatches:
                    stats['nSkipped'] += 1
                    seenTypes = set()
                    for path in missing:
                        label = fileTypeLabel(path)
                        if label not in seenTypes:
                            stats['skipCounts'][label] += 1
                            seenTypes.add(label)
                        if verbose:
                            stats['skipFiles'][label].append(path)
                    if mismatches:
                        stats['nMismatch'] += 1
                        stats['mismatchMsgs'].extend(mismatches)
                    continue
                frame = orbitFrameFromLine(fields, source)
                if frame:
                    stats['framesInOutput'].add(frame)
                block.append(renderLine(raw, fields, source))
            stats['nKept'] += len(block)
            group['blocks'].append(block)
    if executor is not None:
        executor.shutdown()
    return groups, stats


def emitGroup(group, source):
    '''Render one group as output entries: an optional group comment, then the
    data lines of each file with a ';' between consecutive files.'''
    entries = []
    comment = source['groupComment']
    if comment == 'track' and group['trackNum'] is not None:
        entries.append(f'; {group["label"]}')
    elif comment == 'file':
        entries.append(f'; {group["label"]}')
    for index, block in enumerate(group['blocks']):
        if index > 0:
            entries.append(';')
        entries.extend(block)
    return entries


def orderGroups(collected, output, sources):
    '''Return [(group, source), ...] in the output's requested order.

    bySource concatenates the sources in the order the output lists them.
    byTrack merges every group into one sequence sorted by track number, ties
    broken by source order -- which is what interleaves two project roots that
    both hold the same track.
    '''
    pairs = []
    for sourceName in output['sources']:
        source = sources[sourceName]
        for group in collected[sourceName][0]:
            pairs.append((group, source))
    if output['ordering'] == 'byTrack':
        # Stable, so equal track numbers keep source order.
        pairs.sort(key=lambda pair: (pair[0]['trackNum'] is None,
                                     pair[0]['trackNum'] or 0))
    return pairs


def writeMaster(entries, output, check):
    '''Write the assembled entries to output['path']/output['fileName'].

    Data lines are separated by a ';' line, matching what both the NISAR tool
    and the notebook produced; mosaic3d's getDataString() skips any line
    holding a ';'.
    '''
    nData = sum(1 for entry in entries if not entry.startswith(';'))
    outFile = os.path.join(output['path'], output['fileName'])
    if check:
        print(f'{output["name"]}: would write {nData} data lines to {outFile}')
        return nData
    if not os.path.isdir(output['path']):
        os.makedirs(output['path'], exist_ok=True)
        print(f'Created {output["path"]}')
    # A string is written verbatim, so a config can pin the exact literal
    # (e.g. '.2' rather than the '0.2' a YAML float would render as).
    area = output['outputArea']
    if not isinstance(area, str):
        area = ' '.join(str(value) for value in area)
    with open(outFile, 'w') as fp:
        if output['headerComments']:
            fp.write('; Output area\n')
        fp.write(f'{area}\n')
        if output['headerComments']:
            fp.write('; Number of inputs\n')
        else:
            fp.write(';\n')
        fp.write(f'{nData}\n')
        fp.write(';\n')
        for index, entry in enumerate(entries):
            fp.write(entry + '\n')
            if index < len(entries) - 1:
                fp.write(';\n')
    print(f'Wrote {nData} data lines to {outFile}')
    return nData


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def printSkipBreakdown(stats, verbose):
    '''The shared skip/mismatch tail of both summary styles.'''
    if stats['nSkipped']:
        print(f'{YELLOWBG}Skip breakdown:{RESET}')
        for label, count in sorted(stats['skipCounts'].items(),
                                   key=lambda item: -item[1]):
            print(f'  {count:4d} missing {label}')
            if verbose:
                for path in sorted(set(stats['skipFiles'][label])):
                    print(f'    {path}')
    if stats['nMismatch']:
        print(f'{REDBG}Baseline/VRT ionosphere mismatches '
              f'(re-fit baseline): {stats["nMismatch"]}{RESET}')
        for msg in sorted(set(stats['mismatchMsgs'])):
            print(f'    {msg}')


def mergeStats(statsList):
    '''Combine per-source stats into one set of totals.'''
    total = newStats()
    for stats in statsList:
        for key in ['nTotal', 'nSkipped', 'nExcluded', 'nMismatch',
                    'nTrackEntries', 'nKept']:
            total[key] += stats[key]
        total['skipCounts'].update(stats['skipCounts'])
        for label, paths in stats['skipFiles'].items():
            total['skipFiles'][label].extend(paths)
        for key in ['mismatchMsgs', 'tracksNoInputFile', 'excludedDirs',
                    'fileTimes']:
            total[key].extend(stats[key])
        total['framesInOutput'].update(stats['framesInOutput'])
    return total


def printSourceSummary(sources, collected, verbose):
    '''Per-source table, printed once before the outputs are written.'''
    print(f'\n{"source":<16} {"files":>6} {"lines":>8} {"kept":>8} '
          f'{"excluded":>9} {"skipped":>8}')
    for name, source in sources.items():
        groups, stats = collected[name]
        nFiles = sum(len(group['files']) for group in groups)
        print(f'{name:<16} {nFiles:>6} {stats["nTotal"]:>8} '
              f'{stats["nKept"]:>8} {stats["nExcluded"]:>9} '
              f'{stats["nSkipped"]:>8}')
        if stats['tracksNoInputFile']:
            print(f'  {len(stats["tracksNoInputFile"])} of '
                  f'{stats["nTrackEntries"]} tracks have no inputFile'
                  + (f': {stats["tracksNoInputFile"]}' if verbose else ''))
        if verbose:
            for path, when in stats['fileTimes']:
                print(f'  {when:%Y-%m-%d %H:%M}  {path}')
    print()


def printLegacySummary(sources, collected, args):
    '''Reproduce the original NISAR tool's stdout, line for line.'''
    stats = mergeStats([collected[name][1] for name in sources])
    trackNumbers = sorted(set(group['trackNum']
                              for name in sources
                              for group in collected[name][0]
                              if group['trackNum'] is not None)
                          | set(stats['tracksNoInputFile']))
    print(f'Found tracks: {trackNumbers}')
    motionDirs = []
    for source in sources.values():
        if source['frameStatPattern']:
            motionDirs.extend(u.globSort(os.path.join(
                source['root'], source['frameStatPattern'])))
    print(f'Found {len(motionDirs)} track-*/*_0???/motion directories')
    print(f'Tracks with no inputFile (tiepoints not yet run): '
          f'{len(stats["tracksNoInputFile"])} of {stats["nTrackEntries"]}')
    if args.verbose and stats['tracksNoInputFile']:
        print(f'  tracks: {stats["tracksNoInputFile"]}')
    print(f'Data lines found in inputFiles: {stats["nTotal"]}  '
          f'kept: {stats["nTotal"] - stats["nSkipped"]}  '
          f'skipped (missing files): {stats["nSkipped"]}')
    print(f'Unique frames in output: {len(stats["framesInOutput"])} '
          f'of {len(motionDirs)} motion dirs')
    printSkipBreakdown(stats, args.verbose)


# --------------------------------------------------------------------------

def main():
    '''main for makeMaster'''
    args = parseArgs()
    legacy = (args.config is None
              and (args.projectDir is not None or args.inputPaths is not None
                   or args.outputPath is not None
                   or not os.path.exists('masterInput.yaml')))
    if legacy:
        sources, outputs = legacyConfig(args)
    else:
        config = readConfig(args.config or 'masterInput.yaml')
        sources = resolveSources(config)
        outputs = resolveOutputs(config, sources)
        if args.outputs:
            outputs = [o for o in outputs if o['name'] in args.outputs]
            if not outputs:
                u.myerror(f'no output matches {" ".join(args.outputs)}')
        if args.sources:
            for output in outputs:
                output['sources'] = [s for s in output['sources']
                                     if s in args.sources]
    # Collect each source once, however many outputs use it.
    needed = [name for name in sources
              if any(name in output['sources'] for output in outputs)]
    collected = {name: collectSource(sources[name], args.workers, args.verbose)
                 for name in needed}

    if legacy:
        printLegacySummary({n: sources[n] for n in needed}, collected, args)
    else:
        printSourceSummary({n: sources[n] for n in needed}, collected,
                           args.verbose)

    for output in outputs:
        entries = []
        for group, source in orderGroups(collected, output, sources):
            entries.extend(emitGroup(group, source))
        writeMaster(entries, output, args.check)

    if not legacy:
        stats = mergeStats([collected[name][1] for name in needed])
        printSkipBreakdown(stats, args.verbose)
        print('\nCheck the result with:')
        for output in outputs:
            print(f'  checkinput.py '
                  f'{os.path.join(output["path"], output["fileName"])}')


if __name__ == '__main__':
    main()
