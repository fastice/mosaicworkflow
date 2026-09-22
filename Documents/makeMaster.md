# makeMaster

Assembles the master `inputFile`(s) that `mosaic3d` reads, gathering data lines
from several source regions and verifying every file they reference.

Replaces both `nisargrimpworkflow/makeMaster.py` (NISAR) and the Sentinel-1
`Release/masterInput/createMaster.ipynb` notebook.

```
makeMaster [config.yaml] [--outputs NAME ...] [--sources NAME ...]
           [--workers N] [--check] [--verbose]

makeMaster --projectDir DIR [--inputPaths PATH ...] [--outputPath PATH]
```

## Two modes

**Config mode** reads a YAML file (default `masterInput.yaml` in the current
directory) that names the sources and the outputs built from them. Any number
of outputs can be defined over any subset of sources without code changes.

**Legacy mode** is the original NISAR behaviour, kept unchanged so existing
callers keep working. It is selected when any of `--projectDir`,
`--inputPaths` or `--outputPath` is given, or when no config is named and no
`masterInput.yaml` exists. `nisargrimpworkflow.autoupdate` uses this path via
the `masterInputPaths` / `masterOutputPath` keys of `autoupdate.yaml`.

Both modes produce byte-identical output to the tools they replace — see
*Verification* below.

## Config file

```yaml
defaults:        # merged under every source and output; overridable per entry
sources:         # ordered list, each with a unique name
  - {name: ..., layout: trackTree|flat|verbatim, root|file: ..., ...}
groups:          # optional named bundles of source names
  offsets: [S1-insar9, S1A-insar9, ...]
outputs:
  - {name: ..., path: ..., sources: [name-or-group, ...]}
```

Relative paths are resolved against the config file's own directory.

### Layouts

| Layout | Discovers | Used for |
|---|---|---|
| `trackTree` | `root/track-N/tiepoints/track-N/<filePattern>`, one group per track | Sentinel-1 and NISAR project roots |
| `flat` | `root/<filePattern with {year}>` over a year range | TSX / TSXwinter tiepoint directories |
| `verbatim` | an already-assembled master, spliced in as-is | the phase master on insar11 |

A `trackTree` track with no matching inputFile is warned about and skipped —
the normal case, not an edge case (insar10 `Sentinel1-S1D` currently has none
of its 30 tracks tied).

### Source keys

| Key | Default | Meaning |
|---|---|---|
| `filePattern` | — | `{track}` / `{year}` substituted; globs allowed |
| `years` | `{first: 2014, last: 2029}` | `flat` only; or an explicit list |
| `tracks` | all | `trackTree` only; optional whitelist |
| `minFields` | `9` | drop lines with fewer fields |
| `lineFilter` | `fields` | `fields` (the notebook's: any `;` disqualifies) or `heuristic` (the NISAR tool's) |
| `normalizeSlashes` | `true` | `//` → `/`, then re-join fields on single spaces |
| `excludeMarker` | `Exclude` | marker file beside the geodat; `null` disables |
| `excludeLevels` | `1` | directory levels above the geodat to test |
| `checkFields` | `[geodat, azOff, azParams, rgOff, rParams]` | which fields get an existence check; `all` scans every file-like token; `[]` checks nothing |
| `nonFileTokens` | `[none, nophase, nil]` | tokens that are placeholders, not paths |
| `vrtFallback` | `true` | accept `<tok>.vrt` when `<tok>` is absent |
| `checkIonosphere` | `false` | baseline ↔ `range.offsets` VRT consistency |
| `groupComment` | `null` | `track` emits `; track-N`, `file` emits `; <basename>` |
| `frameStatPattern` | `null` | e.g. `track-*/*_0???/motion`, for the frame count |

### Output keys

| Key | Default | Meaning |
|---|---|---|
| `path` | required | directory the `inputFile` is written into |
| `sources` | required | source or group names, in emission order |
| `fileName` | `inputFile` | |
| `outputArea` | required | the six numbers of the first data line |
| `headerComments` | `true` | emit `; Output area` / `; Number of inputs` |
| `ordering` | `bySource` | `bySource` concatenates sources in order; `byTrack` merges every source into one track-number-sorted sequence, ties broken by source order |

`outputArea` given as a **string** is written verbatim. Use that when the exact
literal matters — a YAML float `.2` would otherwise render as `0.2`.

## Checks

A data line is dropped, and counted in the summary, if:

- its frame carries the `excludeMarker` file beside the geodat. On Sentinel-1
  this is a safety net rather than the primary filter: the tie scripts already
  omit Excluded frames when they build `inputFile-N-all`, so a full run
  currently reports zero exclusions even though 255 markers exist on disk;
- any checked field names a file that does not exist (accepting a `.vrt`
  sidecar for a bare name, as `getROffsets.c` resolves it);
- with `checkIonosphere`, the range baseline's `offsetCorrectionFile` disagrees
  with the `range.offsets` VRT's `ionosphereRangeOffsetCorrection` metadata —
  the consistency `mosaic3d`/`checkForIonosphereCorrection()` aborts on at
  runtime, so the line is caught here instead of failing the mosaic.

`--verbose` lists every dropped path. `--check` runs everything and prints the
summary without writing any file.

Validate a finished file with `checkinput.py <path>/inputFile`, which the tool
prints as a suggested next step. It is a separate program and is not run
automatically.

**`checkinput.py` currently crashes on a Sentinel-1 master**, and the fault is
in `checkinput.py`, not in the file. It calls `os.path.getsize()` on the bare
`azimuth.offsets` / `range.offsets` name with no `.vrt` fallback, so it dies
with `FileNotFoundError` on the first VRT-migrated frame — 10,435 of the
20,206 lines in a current full run. `mosaic3d` itself resolves those names
correctly: `readOffsets.c:926` calls `checkForOffsetsVrt()`, which appends
`.vrt` and uses it if present (`gdalIO.c:133`). `checkinput.py` also prints
spurious `Missing .../az.est.yaml.const` and `.../rBaseline.deltabp.yaml.deltabp`
warnings, because its `svConst` handling appends `.const` / `.deltabp` to names
that already carry the newer `.yaml` suffix. Both are pre-existing staleness in
that script, unrelated to this tool.

## `nophase` — a bug fixed on 2026-09-22

`nophase` is the placeholder in the `phase` field of an offsets-only product.
The tool this replaced treated it as a file path, found no such file, and threw
away **the whole data line** — so the offsets went too.

That silently defeated the `ExcludePhase.pending` demotion. When a frame's
*phase baseline* fit fails, the pipeline writes that marker and drops the phase
while keeping the offsets — the marker says so in as many words:

```
2026-09-20 09:27:49: no usable phase baseline (sigma -1.0, 0 tiepoints)
                     -- phase dropped, offsets kept
```

The offsets were then discarded at final assembly anyway. Measured on the
Antarctic master the day it was fixed, it cost 283 lines across 73 tracks, and
the excluded products were in no way defective:

- all 283 have `phase.uw.interp.vrt` built on disk — the phase exists, it is
  the *baseline* that has no solution (257 of them `sigma -1.0, 0 tiepoints`);
- 273 have `ionosphereRangeOffsetCorrection` metadata on their
  `range.offsets.vrt` and 208 a non-nil `offsetCorrectionFile`, so "no
  ionosphere correction" was never the reason;
- **zero** were true duplicates. 217 frames were absent from the master
  entirely; the other 66 appeared only from the *other* bandwidth (Antarctica40
  vs Antarctica80), which is exactly what combining the two roots is for.

`nophase` is now in the default `nonFileTokens` for every source, NISAR
included. The Antarctic master gains those 283 lines and loses nothing — the
only other change is the header count, 5502 → 5785; unique frames go 4544 →
4753. A config that wants the old behaviour back can set
`nonFileTokens: [none]`, but it should not.

## Caveats

- **`makeall.py` is not run.** The notebook regenerated the TSX
  `inputFile-all<year>` files before reading them; that script is no longer on
  `PATH`. They are treated as pre-existing input. `--verbose` prints each
  one's modification time so a stale year is visible.
- **`vrtFallback` on Sentinel-1 is why this tool keeps twice as much data as
  the notebook.** The notebook used a plain `os.path.exists` on the bare
  `range.offsets` / `azimuth.offsets` names. Since those products were migrated
  to tiff/VRT storage there is no bare file on disk any more, only a `.vrt`
  sidecar — which is what `getROffsets.c` resolves at runtime. Measured on
  2026-09-22 over the notebook's own six roots plus TSX: of 20,120 candidate
  lines the notebook keeps 9,771 and this tool keeps all 20,120. Every one of
  the 10,349 extra lines has *both* offsets present only as `.vrt`, and no line
  the notebook keeps is lost. Run today, the notebook would silently drop 51%
  of the archive.
- **The notebook's unused date lookup is gone.** It called
  `utilities.geodatrxa` once per line and discarded the result; removing it
  eliminates ~20,000 geodat parses over NFS. Its side effect of failing loudly
  on a corrupt geodat goes with it; `checkinput.py` still parses dates.
- **Ordering differs from the notebook.** The notebook listed tracks in
  collation order (`track-1, track-10, track-100, track-11`); this tool sorts
  numerically. The set of lines is unchanged.

## Verification

Both legacy consumers were reproduced byte-for-byte on 2026-09-22, output and
stdout alike:

```
makeMaster --projectDir .../Antarctica40 --inputPaths .../Antarctica40 .../Antarctica80 --outputPath <old>
makeMaster --projectDir .../newGreenlandProject --outputPath <old>
```

and the Antarctic config file reproduces its own legacy invocation exactly,
which is what proves `ordering: byTrack` interleaves two project roots holding
the same track number the way the original did.
