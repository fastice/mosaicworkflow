#!/usr/bin/env python3
'''
Figure for Documents/s1PhaseShift2026.md: the 2026 Sentinel-1 constellation
change over our Greenland tracks, and what it does to the mosaic date grid.

Reads satellite, date and absolute orbit from the file names in the S1 archive
(no downloads), derives each pass's track with fileS1.computeTrack, and keeps
the tracks that have a track-N directory in any of the given project roots.

  python3 s1PhaseShift2026.py [--archive DIR] [--projects DIR ...] [--out PNG]
'''
import argparse
import glob
import os
import re
from collections import Counter, defaultdict
from datetime import datetime, timedelta

import matplotlib
matplotlib.use('Agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import mosaicfunc as mosf
import numpy as np
from asfsearchdownload.fileS1 import computeTrack

parser = argparse.ArgumentParser(description='S1 2026 phase-shift figure')
parser.add_argument('--archive', default='/Volumes/insar4/ian/Data/S1-Greenland')
parser.add_argument('--projects', nargs='+',
                    default=['/Volumes/insar9/ian/Sentinel1',
                             '/Volumes/insar10/ian/Sentinel1'])
parser.add_argument('--out', default='s1PhaseShift2026.png')
args = parser.parse_args()

EPOCH = datetime(2015, 1, 1)
T0, T1 = datetime(2026, 1, 1), datetime(2026, 9, 22)
REF_END = datetime(2026, 4, 1)          # each track's reference phase: Jan-Mar

SAT = {'S1A': dict(color='#2a78d6', marker='o', label='S1A'),
       'S1C': dict(color='#eb6834', marker='s', label='S1C'),
       'S1D': dict(color='#1baf7a', marker='^', label='S1D')}
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
UNCHANGED, SHIFTED, OTHER = '#d4d3cf', '#4a3aa7', '#eda100'

ours = {int(os.path.basename(t).split('-')[1])
        for p in args.projects for t in glob.glob(f'{p}/track-*')}
pat = re.compile(r'^(S1[ACD])_IW_SLC__\w{4}_(\d{8})T\d{6}_\d{8}T\d{6}_'
                 r'(\d{6})_\w+_\w{4}\.zip(\.1)?$')
acq = set()
for monthDir in glob.glob(f'{args.archive}/20??-??'):
    for name in os.listdir(monthDir):
        m = pat.match(name)
        if not m:
            continue
        sat, date, orbit = m.group(1), m.group(2), int(m.group(3))
        date = datetime.strptime(date, '%Y%m%d')
        track = computeTrack(orbit, sat)
        if T0 <= date <= T1 and track in ours:
            acq.add((sat, track, date, orbit))
acq = sorted(acq)
tracks = sorted({a[1] for a in acq})
ypos = {t: i for i, t in enumerate(tracks)}

# ---- reference phase per track, and each pass's offset from it ----------
def res(d):
    return (d - EPOCH).days % 6
ref = {}
for t in tracks:
    c = Counter(res(d) for s, tt, d, o in acq if tt == t and d < REF_END)
    ref[t] = c.most_common(1)[0][0]
def offset(t, d):
    k = (res(d) - ref[t]) % 6
    return k - 6 if k > 3 else k         # show as -2..3, expecting 0 and +1

# ---- events, measured from the data --------------------------------------
firstD = min(d for s, t, d, o in acq if s == 'S1D')
lastA = max(d for s, t, d, o in acq if s == 'S1A')
man = min((d for s, t, d, o in acq if s == 'S1C' and o > 8018), default=None)
preMan = max((d for s, t, d, o in acq if s == 'S1C' and o <= 8018), default=None)

# ---- 12-day mosaic windows and each track's first day in each ------------
byTrack = defaultdict(list)
for s, t, d, o in acq:
    byTrack[t].append(d)


def windowState(wins):
    '''Per track and window: 0 first pass on the same day of the window as in
    Jan-Mar, 1 a day later, 2 a day earlier, 3 other; NaN no pass.'''
    state = np.full((len(tracks), len(wins)), np.nan)
    for i, t in enumerate(tracks):
        ks = []
        for w in wins:
            ds = [d for d in byTrack[t] if w['date1'] <= d <= w['date2']]
            ks.append((min(ds) - w['date1']).days if ds else None)
        k0 = Counter(k % 6 for k, w in zip(ks, wins)
                     if k is not None and w['date1'] < REF_END).most_common(1)
        if not k0:
            continue
        for j, k in enumerate(ks):
            if k is not None:
                state[i, j] = {0: 0, 1: 1, 5: 2}.get((k - k0[0][0]) % 6, 3)
    return state


def grid(**kw):
    return [w for w in mosf.standardDates(**kw)
            if w['date1'] >= datetime(2025, 12, 28) and w['date2'] <= T1]


winsOld, winsNew = grid(phaseShifts=[]), grid()
stateOld, stateNew = windowState(winsOld), windowState(winsNew)

# ---- figure ---------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'axes.edgecolor': GRID, 'axes.labelcolor': INK2,
                     'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.titlecolor': INK, 'text.color': INK})
fig, (a1, a2, a3, a4) = plt.subplots(
    4, 1, figsize=(13, 17.5), sharex=True, facecolor=SURF,
    gridspec_kw=dict(height_ratios=[4.2, 1.5, 2.6, 2.6], hspace=0.16))
fig.subplots_adjust(top=0.905)
for ax in (a1, a2, a3, a4):
    ax.set_facecolor(SURF)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

# events on every panel
def events(ax, label=False):
    ax.axvspan(firstD, lastA, color='#efeeea', zorder=0, lw=0)
    ax.axvspan(preMan + timedelta(days=1), man, color='#dcdad4', zorder=0.5,
               lw=0)
    for d in (firstD, lastA):
        ax.axvline(d, color=INK2, lw=1, ls=(0, (4, 3)), zorder=1)
    if label:
        tr = ax.get_xaxis_transform()
        ax.text(firstD, 1.005, f'first S1D\n{firstD:%b %d}', transform=tr,
                ha='center', va='bottom', fontsize=8.5, color=INK2)
        ax.text(man, 1.005, f'S1C manoeuvre: no S1C\n'
                f'{preMan + timedelta(days=1):%b %d} – {man - timedelta(days=1):%b %d}  ',
                transform=tr, ha='right', va='bottom', fontsize=8.5,
                color=INK2)
        ax.text(lastA, 1.005, f'  last S1A\n  {lastA:%b %d}', transform=tr,
                ha='left', va='bottom', fontsize=8.5, color=INK2)

# A: every pass on every track
events(a1, label=True)
for s, st in SAT.items():
    pts = [(d, ypos[t]) for ss, t, d, o in acq if ss == s]
    a1.scatter([p[0] for p in pts], [p[1] for p in pts], s=16,
               color=st['color'], marker=st['marker'], lw=0.4,
               edgecolor=SURF, zorder=3, label=st['label'])
a1.set_yticks(range(len(tracks)))
a1.set_yticklabels([str(t) for t in tracks], fontsize=7)
a1.set_ylim(-0.8, len(tracks) - 0.2)
a1.set_ylabel('track (relative orbit)')
a1.grid(axis='x', color=GRID, lw=0.6)
a1.set_title('A   Every 2026 pass over our 42 Greenland tracks', loc='left',
             fontsize=10.5, fontweight='bold', pad=40)
h, l = a1.get_legend_handles_labels()
h += [Patch(color='#efeeea', label='S1A, S1C and S1D all acquiring')]
a1.legend(handles=h, loc='lower right', bbox_to_anchor=(1.0, 1.085), ncol=4,
          frameon=False, handletextpad=0.3, columnspacing=1.2, fontsize=9,
          markerscale=1.6)

# B: phase offset from each track's Jan-Mar phase
events(a2)
jit = {'S1A': -0.17, 'S1C': 0.0, 'S1D': 0.17}
for s, st in SAT.items():
    pts = [(d, offset(t, d) + jit[s]) for ss, t, d, o in acq if ss == s]
    a2.scatter([p[0] for p in pts], [p[1] for p in pts], s=14,
               color=st['color'], marker=st['marker'], lw=0.3,
               edgecolor=SURF, zorder=3)
a2.set_yticks([0, 1])
a2.set_yticklabels(['on the old\ngrid (0)', '+1 day'])
a2.set_ylim(-0.6, 1.6)
a2.grid(axis='y', color=GRID, lw=0.6)
a2.set_title('B   Acquisition day relative to each track\'s Jan–Mar phase '
             '(date mod 6)', loc='left', fontsize=10.5, fontweight='bold')
for s, y, x, dy in [('S1A', -0.17, datetime(2026, 2, 1), -14),
                    ('S1C', 0.0, datetime(2026, 3, 1), 9),
                    ('S1D', 1.17, datetime(2026, 5, 12), 9),
                    ('S1C', 1.0, datetime(2026, 8, 10), -14)]:
    a2.annotate(s, (x, y), xytext=(0, dy), textcoords='offset points',
                ha='center', va='center', fontsize=8.5, color=INK2)

# C and D: the 12-day mosaic windows, before and after the fix
COLORS = [UNCHANGED, SHIFTED, OTHER, '#e34948']
LABELS = ['same day as Jan–Mar', 'one day later', 'one day earlier', 'other']


def windowPanel(ax, wins, state, title, legend):
    x0 = mdates.date2num([w['date1'] for w in wins])
    x1 = mdates.date2num([w['date2'] + timedelta(days=1) for w in wins])
    for j in range(len(wins)):
        for i, v in enumerate(state[:, j]):
            if not np.isnan(v):
                ax.add_patch(plt.Rectangle((x0[j] + 0.35, i - 0.42),
                                           x1[j] - x0[j] - 0.7, 0.84,
                                           color=COLORS[int(v)], lw=0))
    events(ax)
    ax.set_ylim(-0.8, len(tracks) - 0.2)
    ax.set_yticks(range(len(tracks)))
    ax.set_yticklabels([str(t) for t in tracks], fontsize=6)
    ax.set_ylabel('track')
    ax.set_title(title, loc='left', fontsize=10.5, fontweight='bold')
    if legend:
        present = set(np.concatenate([stateOld[~np.isnan(stateOld)],
                                      stateNew[~np.isnan(stateNew)]]).astype(int))
        ax.legend(handles=[Patch(color=COLORS[k], label=LABELS[k])
                           for k in sorted(present)]
                  + [Patch(facecolor=SURF, edgecolor=GRID, label='no pass')],
                  loc='upper left', bbox_to_anchor=(0, -0.09), ncol=5,
                  frameon=False, fontsize=9)


long = next(w for w in winsNew if (w['date2'] - w['date1']).days == 12)
windowPanel(a3, winsOld, stateOld,
            'C   Old grid: day within each 12-day window on which the track '
            'is first acquired', legend=False)
windowPanel(a4, winsNew, stateNew,
            f"D   New grid, one 13-day window ({long['date1']:%b %d} – "
            f"{long['date2']:%b %d}): every window from Jul 9 on is back in "
            'sync', legend=True)
a4.axvspan(long['date1'], long['date2'] + timedelta(days=1), ymin=0, ymax=1,
           fill=False, edgecolor=INK, lw=1.2, zorder=4)
a4.xaxis.set_major_locator(mdates.MonthLocator())
a4.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
a4.set_xlim(T0, T1)

fig.suptitle('Sentinel-1 over Greenland in 2026: the constellation moved one '
             'day later', x=0.125, y=0.99, ha='left', fontsize=14,
             fontweight='bold')
fig.text(0.125, 0.972,
         'S1D arrived one day behind S1C; S1C was then manoeuvred onto S1D\'s '
         'phase and S1A retired, so from July every track is acquired one day '
         'later than before.\nThe 12-day mosaic grid did not move, so each '
         'window caught its tracks a day later (C) until one window was made '
         '13 days long (D).  Source: file names in '
         f'{args.archive}.',
         ha='left', va='top', fontsize=9.5, color=INK2)
out = args.out
fig.savefig(out, dpi=130, facecolor=SURF, bbox_inches='tight')
print('wrote', out)
print(f'first S1D {firstD:%Y-%m-%d}  last S1C<=8018 {preMan:%Y-%m-%d}  '
      f'first S1C>8018 {man:%Y-%m-%d}  last S1A {lastA:%Y-%m-%d}')
print('old grid states:', Counter(stateOld[~np.isnan(stateOld)].astype(int).tolist()))
print('new grid states:', Counter(stateNew[~np.isnan(stateNew)].astype(int).tolist()))
