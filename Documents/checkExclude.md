# checkExclude

Reviews frames marked `Exclude.pending` by comparing each frame's velocity with
its track's velocityStats reference, then — once you approve — turns every
pending frame into a hard `Exclude` (FAIL) or a `checked` marker (PASS), with the
statistics written inside the file. Hard `Exclude` frames are only annotated.

Works on Sentinel-1 (XY: `vx/vy`, `ex/ey`) and NISAR (RA, line of sight:
`vr/va`, `er/ea`), chosen by `velocityStatsMode` in `project.yaml`.

```
checkExclude --dryRun                        # report only
checkExclude                                 # report, then ask before changing anything
checkExclude --tracks track-26 track-90 --yes
checkExclude --replay                        # re-show the latest run + its changes
checkExclude --replay list                   # every run: verdict counts, applied or not
checkExclude --replay 20261010T1513          # one run by its stamp
```

`--replay` changes nothing and runs no comparison: it prints the run's report
exactly as it was shown at the start (`logs/checkExclude-<stamp>.md`), its verdict
counts, and then the changes that run applied, from `logs/checkExclude-changes.log`.

Run from a project root: a Sentinel-1 prime, a Sentinel-1 clone (`Sentinel1-S1C`
etc. — its `velocityStats` is the prime's, through the symlink), or a NISAR project.

| Option | Default | Meaning |
|---|---|---|
| `--project DIR` | `.` | project root |
| `--tracks` | all `track-*` | restrict to these tracks |
| `--dryRun` | off | report only, never change markers |
| `--yes` | off | apply without the confirmation prompt |
| `--noAnnotate` | off | leave hard `Exclude` files alone (a file already carrying checkExclude stats is never annotated twice) |
| `--firstDate` / `--lastDate` | all | only frames whose pair starts in this window (from `*.pairinfo`) |
| `--replay [latest\|list\|STAMP]` | — | show a past run and the changes it applied, instead of running |
| `--minOverlap N` | 1000 | fewer overlapping pixels → UNTESTABLE |
| `--minWithin3Sigma F` | 0.80 | PASS needs at least this fraction within 3σ |
| `--maxMedNorm F` | 1.5 | PASS needs every component's median \|v−mean\|/σ at or below this |
| `--sigmaFloor S` | 5.0 m/yr | σ = max(reference sigma, this) |
| `--minDays N` | 2 | pairs shorter than this are SKIPPED: listed, never judged, markers never changed (the 1-day S1C–S1D tandem pairs are too noisy for offsets) |
| `--minTiepoints N` | 4 | PASS needs at least this many tiepoints in the range fit (and the azimuth fit when `requireAzSolution`); 3 points fit the model exactly, so their sigma of 0 means no redundancy |

The thresholds can also be set per project in `project.yaml`; flags win:

```yaml
checkExclude:
  minWithin3Sigma: 0.80
  maxMedNorm: 1.5
  requireAzSolution: true     # default true for Sentinel-1, false for NISAR
```

## What it measures

For every frame with `Exclude.pending` (and, unless `--noAnnotate`, every hard
`Exclude`):

- **Tiepoints and fits** — `motion/az.est.yaml` and the range-baseline yaml the
  frame's `velocity/inputFile` uses (`rBaseline.deltabp.yaml` / `rBaseline.yaml`):
  sigma, tiepoints used, tiepoints given.
- **Velocity agreement** — `velocity/mosaicOffsets.<c>.tif` against
  `velocityStats/<range>/velocity.<c>.tif`, over pixels where both are valid and
  the reference has `navg >= 1`:
  - `within3Sigma`: fraction within 3σ on both components;
  - `medNorm_<c>`: median \|v − mean\| / σ per component;
  - `bias_<c>`: median (v − mean) per component, m/yr.

## Verdicts

| verdict | when | applied as |
|---|---|---|
| FAIL | no velocity data (the frame has no usable baseline solution, so mosaic3d wrote a blank map); or the agreement test fails; or, for Sentinel-1, the azimuth fit has no solution (setupS1Tracks step 7 would hard-Exclude it at the next update anyway) | `Exclude` = the original pending text + `<datetime>: checkExclude FAIL: <why>` + stats; `Exclude.pending` removed |
| PASS | overlap ≥ minOverlap, within3Sigma ≥ minWithin3Sigma, every medNorm ≤ maxMedNorm, and the fits use ≥ minTiepoints | `checked` = `<datetime>: checkExclude PASS` + stats; `Exclude.pending` removed |
| UNTESTABLE | too little overlap, no velocityStats range, or grid mismatch | left pending, listed |
| SKIPPED | pair shorter than `minDays` | nothing |
| ANNOTATE | hard `Exclude` | a dated stats block appended; status never changed (an empty hand-placed file, or setupisceuw's `AUTO-EXCLUDED` tag, is kept as is) |

The stats are measured before any conversion, while the pending frame's velocity
still exists — `makevelnoclean` later removes a hard-Excluded frame's products,
so the numbers in the `Exclude` file are the only record of why.

The report goes to stdout and `logs/checkExclude-<stamp>.md`, grouped by verdict,
with tiepoints used/given for both fits. Every applied change is also appended to
`logs/checkExclude-changes.log` (track, frame, dates, before -> after, why) — the
running record of what the tool changed.

## `checked` is conditional

`checked` passes a frame *for the fit it was checked against*. Its stats block
records `nTiepointsAz`, `nTiepointsRange` and the sigmas. The writers of
`Exclude.pending` consult it through `utilities.checkedStillValid`:

- **insarworkflow `tieScript`** (range baseline) — does not re-flag a checked
  frame while its range fit still has a solution and at least as many tiepoints
  as when it was checked;
- **s1setup `setupS1Tracks --fullUpdate` step 7** (azimuth baseline) — same rule
  for the azimuth fit.

**Legacy `checked` files.** Empty `checked` files placed by hand before this tool
(insar9 has 241, from 2021–2025) carry no counts to compare against. They are
honoured while the current fit still has a solution and at least 4 tiepoints
(`LEGACY_MIN_TIEPOINTS`); below that the frame is flagged again.

If a later cleaning cuts tiepoints, or a fit loses its solution, the frame is
flagged again, with `[checked <date>: nR 5->3 (nR dropped)]` appended to the
reason, and the next `checkExclude` re-tests it. A missing offsets file is always
flagged.

## Marker file format

```
2026-10-10 10:26:00.123456: checkExclude PASS
  checkedDate: 2026-10-10
  verdict: PASS
  why: agrees with velocityStats
  mode: XY
  nTiepointsAz: 5.0
  nTiepointsAzGiven: 38466.0
  sigmaAz: 0.0509
  nTiepointsRange: 5.0
  ...
  overlapPixels: 35642
  medNorm_vx: 0.65
  ...
  within3Sigma: 0.98
  threshold_minWithin3Sigma: 0.8
```

## Part of the mosaicworkflow package.
