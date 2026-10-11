#!/usr/bin/env python3
"""
checkExclude -- review Exclude.pending (and annotate hard Exclude) frames
against the track's velocityStats reference, then, on approval, turn each
pending frame into a hard `Exclude` (FAIL) or a `checked` marker (PASS) with
the statistics written inside the file.

Runs from a project root (Sentinel-1 prime or clone, or NISAR). The comparison
uses velocityStatsMode from project.yaml: XY (vx/vy, ex/ey) for Sentinel-1, RA
(line of sight: vr/va, er/ea) for NISAR.

Per pending frame:
  * tie/fit stats: motion/az.est.yaml and the range-baseline yaml named in
    velocity/inputFile (sigma, tiepoints used/given);
  * the frame's velocity/mosaicOffsets.<c1,c2>.tif against
    velocityStats/<range>/velocity.<c1,c2,e1,e2,navg>.tif over the pixels both
    cover (navg >= 1), with sigma = max(e, sigmaFloor):
      within3Sigma  fraction of pixels within 3 sigma on both components
      medNorm_<c>   median |v - mean| / sigma per component
      bias_<c>      median (v - mean) per component (m/yr)
  * verdict:
      FAIL        no velocity data (no baseline solution), or the test fails,
                  or (Sentinel-1) no azimuth solution -- step 7 of
                  setupS1Tracks would hard-Exclude it at the next update anyway
      UNTESTABLE  fewer than minOverlap pixels, no reference, or grid mismatch
                  (left pending)
      PASS        within3Sigma >= minWithin3Sigma and every medNorm <= maxMedNorm
Hard Exclude frames are ANNOTATED only (stats appended, status never changed).

`checked` is conditional: tieScript and setupS1Tracks step 7 re-flag a checked
frame if its tiepoints later drop or its fit loses its solution
(utilities.checkedStillValid).

    checkExclude --dryRun
    checkExclude --tracks track-26 track-90
    checkExclude --yes
    checkExclude --replay [latest | list | 20261010T1513]
"""
import argparse
import datetime as dt
import glob
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import yaml
from osgeo import gdal

import utilities as u
from utilities.checkedMarker import CHECKED_NAME
from mosaicworkflow.autocleanNISAR import findVelStatsDir

gdal.UseExceptions()

FRAME_RE = re.compile(r'^\d+_\d+$')
DEFAULTS = {'minOverlap': 1000, 'minWithin3Sigma': 0.80, 'maxMedNorm': 1.5,
            'sigmaFloor': 5.0, 'minTiepoints': 4}
COMPONENTS = {'XY': (('vx', 'ex'), ('vy', 'ey')),
              'RA': (('vr', 'er'), ('va', 'ea'))}


def readProject(project):
    path = os.path.join(project, 'project.yaml')
    if not os.path.exists(path):
        u.myerror(f'checkExclude: no project.yaml in {project}')
    with open(path) as fp:
        return yaml.safe_load(fp) or {}


def readFit(path):
    """(sigma, nUsed, nGiven) from an az.est / rBaseline yaml; Nones if absent."""
    vals = {}
    try:
        with open(path) as fp:
            for line in fp:
                m = re.match(r'(sigma|nTiepointsUsed|nTiepointsGiven):\s*([-\d.eE+]+)',
                             line)
                if m:
                    vals[m.group(1)] = float(m.group(2))
    except OSError:
        pass
    return (vals.get('sigma'), vals.get('nTiepointsUsed'),
            vals.get('nTiepointsGiven'))


def rangeBaselineFile(frameDir):
    """The range-baseline yaml the frame's velocity run uses (last token of the
    inputFile line), else the usual names."""
    try:
        with open(os.path.join(frameDir, 'velocity', 'inputFile')) as fp:
            for line in fp:
                tok = line.split()
                if tok and tok[-1].endswith('.yaml') and 'rBaseline' in tok[-1]:
                    return tok[-1]
    except OSError:
        pass
    for name in ('rBaseline.deltabp.yaml', 'rBaseline.yaml'):
        path = os.path.join(frameDir, 'motion', name)
        if os.path.exists(path):
            return path
    return None


def pairDates(frameDir):
    for pi in glob.glob(os.path.join(frameDir, '*.pairinfo')):
        try:
            tok = open(pi).read().split()
            return tok[2], tok[3]
        except (OSError, IndexError):
            pass
    return '?', '?'


def markerReason(path):
    """First line of a marker, minus the '<datetime>: ' prefix."""
    try:
        with open(path) as fp:
            line = fp.readline().strip()
    except OSError:
        return ''
    m = re.match(r'\d{4}-\d\d-\d\d \d\d:\d\d:\d\d(\.\d+)?: (.*)', line)
    reason = m.group(2) if m else line
    # absolute paths -> basenames, so the table shows the message
    return re.sub(r'/\S*/(\S+)', r'\1', reason)


def readBand(path):
    """Float array with no-data (NaN, < -1e6, nodata value) as NaN, or None."""
    if not os.path.exists(path):
        return None
    ds = gdal.Open(path)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float64)
    nd = band.GetNoDataValue()
    if nd is not None:
        arr[arr == nd] = np.nan
    arr[arr < -1e6] = np.nan
    return arr


def compareVelocity(frameDir, refDir, comps, sigmaFloor):
    """Comparison stats (dict) of the frame velocity against the reference, or
    a dict with 'untestable'/'noVelocity' explaining why it could not be run."""
    vel = [readBand(os.path.join(frameDir, 'velocity', f'mosaicOffsets.{c}.tif'))
           for c, _ in comps]
    if any(v is None for v in vel) or not any(np.isfinite(v).any() for v in vel):
        return {'noVelocity': 'no velocity data (no usable baseline solution)'}
    mean = [readBand(os.path.join(refDir, f'velocity.{c}.tif')) for c, _ in comps]
    sig = [readBand(os.path.join(refDir, f'velocity.{e}.tif')) for _, e in comps]
    navg = readBand(os.path.join(refDir, 'velocity.navg.tif'))
    if any(a is None for a in mean + sig) or navg is None:
        return {'untestable': f'incomplete reference in {refDir}'}
    if len({a.shape for a in vel + mean + sig + [navg]}) != 1:
        return {'untestable': f'grid mismatch: frame {vel[0].shape} vs '
                              f'reference {mean[0].shape}'}
    ok = (navg >= 1)
    for a in vel + mean + sig:
        ok &= np.isfinite(a)
    n = int(ok.sum())
    stats = {'overlapPixels': n}
    if n == 0:
        return stats
    within = np.ones(n, bool)
    for (c, _), v, m, s in zip(comps, vel, mean, sig):
        d = v[ok] - m[ok]
        norm = np.abs(d) / np.maximum(s[ok], sigmaFloor)
        within &= norm <= 3
        stats[f'medNorm_{c}'] = round(float(np.median(norm)), 3)
        stats[f'bias_{c}'] = round(float(np.median(d)), 2)
    stats['within3Sigma'] = round(float(within.mean()), 3)
    return stats


def evaluateFrame(job, cfg):
    trackDir, frameDir, kind = job
    track, frame = os.path.basename(trackDir), os.path.basename(frameDir)
    d1, d2 = pairDates(frameDir)
    sA, nA, gA = readFit(os.path.join(frameDir, 'motion', 'az.est.yaml'))
    rFile = rangeBaselineFile(frameDir)
    sR, nR, gR = readFit(rFile) if rFile else (None, None, None)
    rec = {'track': track, 'frame': frame, 'kind': kind, 'date1': d1, 'date2': d2,
           'reason': markerReason(os.path.join(frameDir, kind)),
           'stats': {'nTiepointsAz': nA, 'nTiepointsAzGiven': gA, 'sigmaAz': sA,
                     'nTiepointsRange': nR, 'nTiepointsRangeGiven': gR,
                     'sigmaRange': sR}}
    frameNum = int(frame.split('_')[1][len(cfg['framePrefix']):])
    refName = findVelStatsDir(frameNum, cfg['framePattern'], trackDir)
    refDir = os.path.join(trackDir, 'velocityStats', refName)
    if not os.path.isdir(refDir):
        cmp = {'untestable': f'no velocityStats range for frame {frameNum}'}
    else:
        cmp = compareVelocity(frameDir, refDir, cfg['comps'], cfg['sigmaFloor'])
        rec['stats']['reference'] = f'velocityStats/{refName}'
    noVelocity = cmp.pop('noVelocity', None)
    untestable = cmp.pop('untestable', None)
    rec['stats'].update(cmp)
    if kind == 'Exclude':
        rec['verdict'], rec['why'] = 'ANNOTATE', noVelocity or untestable or ''
        return rec
    # Fit problems decide on their own -- no amount of overlap rescues a fit
    # with no solution or no redundancy -- so they are checked first.
    fitProblems = []
    if cfg['requireAzSolution'] and (sA is None or sA < 0):
        fitProblems.append('no azimuth solution (setupS1Tracks step 7 would '
                           'hard-Exclude it)')
    fits = [('range', nR)] + ([('azimuth', nA)] if cfg['requireAzSolution'] else [])
    for label, n in fits:
        if n is None or n < cfg['minTiepoints']:
            fitProblems.append(f'{label} fit from {fmt(n)} tiepoints '
                               f'(< {cfg["minTiepoints"]}, no redundancy)')
    if noVelocity:
        rec['verdict'], rec['why'] = 'FAIL', noVelocity
    elif fitProblems:
        rec['verdict'], rec['why'] = 'FAIL', '; '.join(fitProblems)
    elif untestable:
        rec['verdict'], rec['why'] = 'UNTESTABLE', untestable
    elif cmp['overlapPixels'] < cfg['minOverlap']:
        rec['verdict'] = 'UNTESTABLE'
        rec['why'] = f'only {cmp["overlapPixels"]} overlapping pixels'
    else:
        problems = []
        if cmp['within3Sigma'] < cfg['minWithin3Sigma']:
            problems.append(f'within3Sigma {cmp["within3Sigma"]:.2f} < '
                            f'{cfg["minWithin3Sigma"]}')
        for c, _ in cfg['comps']:
            if cmp[f'medNorm_{c}'] > cfg['maxMedNorm']:
                problems.append(f'medNorm_{c} {cmp[f"medNorm_{c}"]:.2f} > '
                                f'{cfg["maxMedNorm"]}')
        rec['verdict'] = 'FAIL' if problems else 'PASS'
        rec['why'] = '; '.join(problems) if problems else 'agrees with velocityStats'
    return rec


def fmt(v, spec=''):
    if v is None:
        return '-'
    if isinstance(v, float) and v.is_integer() and not spec:
        return str(int(v))
    return format(v, spec) if spec else str(v)


def reportLines(recs, cfg, project, stamp):
    c1, c2 = (c for c, _ in cfg['comps'])
    head = (f'| track | frame | dates | reason | nAz used/given | nR used/given | '
            f'sigmaAz | sigmaR | overlap px | within 3σ | medNorm {c1}/{c2} | '
            f'bias {c1}/{c2} (m/yr) | verdict | why |')
    lines = [f'# checkExclude {stamp}', '',
             f'- project: `{project}`  mode: {cfg["mode"]}',
             f'- thresholds: minOverlap {cfg["minOverlap"]}, minWithin3Sigma '
             f'{cfg["minWithin3Sigma"]}, maxMedNorm {cfg["maxMedNorm"]}, sigmaFloor '
             f'{cfg["sigmaFloor"]} m/yr, minTiepoints {cfg["minTiepoints"]}, '
             f'requireAzSolution {cfg["requireAzSolution"]}', '']
    for verdict in ('PASS', 'FAIL', 'UNTESTABLE', 'ANNOTATE'):
        sel = [r for r in recs if r['verdict'] == verdict]
        lines += [f'## {verdict} ({len(sel)})', '']
        if not sel:
            lines += ['none', '']
            continue
        lines += [head, '|---' * 14 + '|']
        for r in sel:
            s = r['stats']
            lines.append(
                f'| {r["track"]} | {r["frame"]} | {r["date1"]}..{r["date2"]} | '
                f'{r["reason"][:60]} | {fmt(s["nTiepointsAz"])}/{fmt(s["nTiepointsAzGiven"])} | '
                f'{fmt(s["nTiepointsRange"])}/{fmt(s["nTiepointsRangeGiven"])} | '
                f'{fmt(s["sigmaAz"], ".3g")} | {fmt(s["sigmaRange"], ".3g")} | '
                f'{fmt(s.get("overlapPixels"))} | {fmt(s.get("within3Sigma"), ".2f")} | '
                f'{fmt(s.get(f"medNorm_{c1}"), ".2f")}/{fmt(s.get(f"medNorm_{c2}"), ".2f")} | '
                f'{fmt(s.get(f"bias_{c1}"), ".1f")}/{fmt(s.get(f"bias_{c2}"), ".1f")} | '
                f'{r["verdict"]} | {r["why"]} |')
        lines.append('')
    return lines


def statsBlock(rec, cfg, now):
    stats = {'checkedDate': now.strftime('%Y-%m-%d'), 'verdict': rec['verdict'],
             'why': rec['why'] or '-', 'mode': cfg['mode']}
    stats.update({k: v for k, v in rec['stats'].items() if v is not None})
    stats.update({f'threshold_{k}': cfg[k] for k in DEFAULTS})
    return stats


def applyVerdicts(recs, cfg, project, changeLog, stamp):
    """Write the markers, and one line per change to changeLog: the record of
    what this run did to each frame (before -> after, and why)."""
    now = dt.datetime.now()
    action = {'FAIL': 'Exclude.pending -> Exclude', 'PASS': 'Exclude.pending -> checked',
              'ANNOTATE': 'Exclude annotated (status unchanged)'}
    with open(changeLog, 'a') as log:
        print(f'# checkExclude report checkExclude-{stamp}.md applied '
              f'{now:%Y-%m-%d %H:%M:%S} {project}', file=log)
        for r in recs:
            if r['verdict'] in action:
                print(f'{r["track"]}\t{r["frame"]}\t{r["date1"]}..{r["date2"]}\t'
                      f'{action[r["verdict"]]}\t{r["why"] or "-"}', file=log)
    for r in recs:
        frameDir = os.path.join(project, r['track'], r['frame'])
        pending = os.path.join(frameDir, 'Exclude.pending')
        stats = statsBlock(r, cfg, now)
        if r['verdict'] == 'FAIL':
            original = open(pending).read() if os.path.exists(pending) else ''
            with open(os.path.join(frameDir, 'Exclude'), 'w') as fp:
                fp.write(original)
                print(f'{now}: checkExclude FAIL: {r["why"]}', file=fp)
                u.writeStatsBlock(fp, stats)
            os.remove(pending)
        elif r['verdict'] == 'PASS':
            with open(os.path.join(frameDir, CHECKED_NAME), 'w') as fp:
                print(f'{now}: checkExclude PASS', file=fp)
                u.writeStatsBlock(fp, stats)
            os.remove(pending)
        elif r['verdict'] == 'ANNOTATE':
            with open(os.path.join(frameDir, 'Exclude'), 'a') as fp:
                print(f'{now}: checkExclude stats (status unchanged)', file=fp)
                u.writeStatsBlock(fp, stats)


def replay(project, which):
    """Re-show a past run as it was printed at the start (its report) and the
    changes it applied, from logs/checkExclude-<stamp>.md and
    logs/checkExclude-changes.log. which: 'latest', 'list', or a stamp."""
    logs = os.path.join(project, 'logs')
    reports = sorted(glob.glob(os.path.join(logs, 'checkExclude-2*.md')))
    if not reports:
        u.myerror(f'checkExclude --replay: no reports in {logs}')
    stampOf = lambda p: os.path.basename(p)[len('checkExclude-'):-len('.md')]
    changeLog = os.path.join(logs, 'checkExclude-changes.log')
    blocks = {}
    if os.path.exists(changeLog):
        key = None
        for line in open(changeLog):
            if line.startswith('#'):
                m = re.search(r'(\d{8}T\d{4})', line)
                if m:
                    key = m.group(1)
                else:
                    # header with only the apply time: the run is the latest
                    # report written at or before it
                    t = re.search(r'(\d{4}-\d\d-\d\d) (\d\d):(\d\d)', line)
                    when = f'{t.group(1).replace("-", "")}T{t.group(2)}{t.group(3)}' \
                        if t else ''
                    earlier = [stampOf(r) for r in reports if stampOf(r) <= when]
                    key = earlier[-1] if earlier else None
                blocks.setdefault(key, [])
            elif key:
                blocks[key].append(line.rstrip('\n'))
    if which == 'list':
        for r in reports:
            st = stampOf(r)
            head = open(r).read()
            counts = ' '.join(f'{v} {n}' for v, n in
                              re.findall(r'^## (\w+) \((\d+)\)', head, re.M))
            applied = f'{len(blocks[st])} changes applied' if st in blocks \
                else 'not applied (dry run or declined)'
            print(f'{st}  {counts}  -- {applied}')
        return 0
    path = reports[-1] if which == 'latest' else \
        os.path.join(logs, f'checkExclude-{which}.md')
    if not os.path.exists(path):
        u.myerror(f'checkExclude --replay: no report {path} (try --replay list)')
    st = stampOf(path)
    text = open(path).read()
    print(text.rstrip('\n'))
    counts = dict(re.findall(r'^## (\w+) \((\d+)\)', text, re.M))
    print(f'report: {path}')
    print('  '.join(f'{v} {counts.get(v, 0)}' for v in
                    ('PASS', 'FAIL', 'UNTESTABLE', 'ANNOTATE')))
    if st in blocks:
        print(f'\nchanges applied ({len(blocks[st])}, from {changeLog}):')
        for line in blocks[st]:
            print('  ' + line.replace('\t', '  '))
    else:
        print('\nno changes recorded for this run (dry run, declined, or before '
              'the change log existed)')
    return 0


def parseArgs():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog='Part of the mosaicworkflow package.')
    ap.add_argument('--project', default='.', help='project root [.]')
    ap.add_argument('--tracks', nargs='*', default=None,
                    help='track dirs (default every track-*)')
    ap.add_argument('--dryRun', action='store_true',
                    help='report only; never change markers')
    ap.add_argument('--yes', action='store_true',
                    help='apply without the confirmation prompt')
    ap.add_argument('--noAnnotate', action='store_true',
                    help='do not append stats to hard Exclude files')
    for k, v in DEFAULTS.items():
        ap.add_argument(f'--{k}', type=type(v), default=None,
                        help=f'override project.yaml checkExclude.{k} [{v}]')
    ap.add_argument('--firstDate', default=None,
                    help='only frames whose pair starts on/after YYYY-MM-DD')
    ap.add_argument('--lastDate', default=None,
                    help='only frames whose pair starts on/before YYYY-MM-DD')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--replay', nargs='?', const='latest', default=None,
                    metavar='STAMP|list',
                    help='show a past run as it was printed, plus the changes it '
                         'applied, instead of running: latest (default), list, '
                         'or a stamp such as 20261010T1513')
    return ap.parse_args()


def main():
    args = parseArgs()
    project = os.path.abspath(args.project)
    if args.replay:
        return replay(project, args.replay)
    proj = readProject(project)
    mode = str(proj.get('velocityStatsMode', 'XY')).upper()
    framePattern = proj.get('framePattern', '???')
    cfg = dict(DEFAULTS)
    cfg.update(proj.get('checkExclude') or {})
    for k in DEFAULTS:
        if getattr(args, k) is not None:
            cfg[k] = getattr(args, k)
    cfg.update({'mode': mode, 'comps': COMPONENTS[mode],
                'framePattern': framePattern,
                'framePrefix': framePattern.split('?')[0] if '?' in framePattern else '',
                'requireAzSolution': cfg.get('requireAzSolution',
                                             str(proj.get('sensor', '')).upper() != 'NISAR')})

    trackDirs = ([os.path.join(project, t) for t in args.tracks] if args.tracks
                 else sorted(glob.glob(os.path.join(project, 'track-*')),
                             key=lambda p: int(re.sub(r'\D', '', os.path.basename(p)) or 0)))
    jobs = []
    for td in trackDirs:
        for fd in sorted(glob.glob(os.path.join(td, '*_*'))):
            if not (os.path.isdir(fd) and FRAME_RE.match(os.path.basename(fd))):
                continue
            if args.firstDate or args.lastDate:
                d1 = pairDates(fd)[0]
                if d1 == '?' or (args.firstDate and d1 < args.firstDate) or \
                        (args.lastDate and d1 > args.lastDate):
                    continue
            exclude = os.path.join(fd, 'Exclude')
            if os.path.exists(exclude):
                # annotate once: a file this tool already wrote or annotated
                # carries its stats already
                if not args.noAnnotate and 'checkExclude' not in open(exclude).read():
                    jobs.append((td, fd, 'Exclude'))
            elif os.path.exists(os.path.join(fd, 'Exclude.pending')):
                jobs.append((td, fd, 'Exclude.pending'))
    print(f'checkExclude: {len(jobs)} frame(s) in {len(trackDirs)} track(s), '
          f'mode {mode}')
    with ThreadPoolExecutor(max(1, args.threads)) as ex:
        recs = list(ex.map(lambda j: evaluateFrame(j, cfg), jobs))

    stamp = dt.datetime.now().strftime('%Y%m%dT%H%M')
    lines = reportLines(recs, cfg, project, stamp)
    out = os.path.join(project, 'logs', f'checkExclude-{stamp}.md')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as fp:
        fp.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    counts = {v: sum(r['verdict'] == v for r in recs)
              for v in ('PASS', 'FAIL', 'UNTESTABLE', 'ANNOTATE')}
    print(f'report: {out}')
    print(f'PASS {counts["PASS"]}  FAIL {counts["FAIL"]}  '
          f'UNTESTABLE {counts["UNTESTABLE"]}  ANNOTATE {counts["ANNOTATE"]}')

    if args.dryRun or not (counts['PASS'] or counts['FAIL'] or counts['ANNOTATE']):
        return 0
    if not args.yes:
        answer = input(f'Apply: {counts["FAIL"]} -> Exclude, {counts["PASS"]} -> checked, '
                       f'{counts["ANNOTATE"]} Exclude annotated? [y/N] ')
        if answer.strip().lower() not in ('y', 'yes'):
            print('nothing changed')
            return 0
    changeLog = os.path.join(project, 'logs', 'checkExclude-changes.log')
    applyVerdicts(recs, cfg, project, changeLog, stamp)
    print(f'applied; changes appended to {changeLog}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
