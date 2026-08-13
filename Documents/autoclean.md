# autoclean.py

Run from a `track-N/` directory.  Cycles through all per-segment offset directories, compares each segment's velocity map against the `velocityStats` reference mean and sigma, flags outlier pixels and ragged edge pixels, maps them back to range/azimuth offset coordinates, and writes a `badoffsets_auto.list` file that `cleanoff` can apply to the raw offsets.

---

## Usage

```
autoclean [options]
```

Must be run from the `track-N/` directory.

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--nocull` | off | Use `velocity_nocull` stats instead of `velocity` as the reference. |
| `--refresh` | off | Overwrite existing `badoffsets_auto.list` files (default: skip dirs that already have one). |
| `--remove` | off | Delete existing `badoffsets_auto.list` files and stop (combine with `--applyall` to rerun `cleanoff` afterwards). |
| `--applyall` | off | Run `cleanoff` in all offset dirs, even those not updated this pass. |
| `--noapply` | off | Do not run `cleanoff` in any directory regardless. |
| `--firstDate YYYY:MM:DD` | `1900:12:31` | Only process pairs on or after this date. |
| `--lastDate YYYY:MM:DD` | `2100:12:31` | Only process pairs on or before this date. |
| `--sigThresh N` | `3.0` | Flag pixels whose deviation from the reference exceeds `N × sigma`. |
| `--threads N` | `12` | Maximum number of parallel threads for mask computation. |
| `--offDir GLOB` | `./*_*` | Glob pattern to find offset directories (must each contain `azimuth.offsets`). |
| `--dem FILE` | auto | DEM file.  For Greenland, `gimp1` is used before 2015, `gimp2` from 2015 onward.  For other regions the region-default DEM is used. |
| `--region NAME` | `greenland` | Region name.  Overridden by `../region` file if present. |

`--applyall` and `--noapply` are mutually exclusive.

---

## Algorithm

### 1. Setup

- Reads `velocityStats/*-*/velocity{_nocull}.vx` to build per-frame mean velocity and sigma maps (`getVelRef`).
- Globs `<offDirRoot>/azimuth.offsets` to find candidate offset directories; skips any with an `Exclude` file or no velocity directory (`getOffsetDirs`).
- Detects the sensor from the working directory path (`getSensor`): `TSX`, `S1`, or `CSK`.

### 2. Per-segment mask generation (threaded)

For each offset directory, in up to `--threads` parallel threads (`processOffsets`):

1. **Read velocity** — loads `velocity_nocull/mosaicOffsets.vx/.vy` (`getVelAndFiles`).
2. **SigThresh map** — builds a float32 threshold array, optionally non-uniform where a region `.shp` file specifies scale factors (`makeSigThresh2D`).
3. **Initial mask** — marks all pixels with valid data (`makeInitialCullMask`).
4. **Edge cleaning** — morphological open (erosion then dilation) to find ragged fringe pixels that survive erosion but not the open (`cleanMaskEdges`).
5. **Difference mask** — flags pixels where `|vel - mean| > sigThresh × sigma` in either component, then unions with edge points; closes the result with a morphological close (`makeDiffMask`).
6. **lat/lon → r/a** — extracts bad-pixel XY positions from the velocity grid, converts to lat/lon via `geo.xykmtoll`, then calls the `lltora` binary to get range/azimuth pixel coordinates (`computeBadRA`).
7. **RA mask** — maps r/a coordinates into the offset grid coordinate system, sets a binary mask, dilates with a 3×3 kernel to fill gaps, and flattens to a 1D index array (`makeRAMask`).
8. **Write list** — writes flat binary int32 indices to `badoffsets_auto.list` (`writeIndexList`).

A `sigThresh.override` single-line file in the segment directory overrides `--sigThresh` for that segment only.

### 3. Apply cleanoff

After all masks are computed, `cleanoff` is run (threaded) in each directory where a new `badoffsets_auto.list` was written, unless `--noapply` is set.  With `--applyall`, `cleanoff` is run in all dirs regardless.

---

## Key files

| File | Description |
|------|-------------|
| `<offDir>/badoffsets_auto.list` | Flat binary int32 array of bad-pixel indices in offset coordinates |
| `<offDir>/velFile.report` | Per-segment report: mean absolute difference in vx and vy |
| `<offDir>/sigThresh.override` | Optional single-line float to override `--sigThresh` for this segment |

---

## Notes

- Directories containing an `Exclude` file are always skipped.
- The `velocity_nocull` map is always used as the per-segment input (even when `--nocull` is not set), because culled maps have ragged edges that cause under-culling on a second pass.  `--nocull` controls only which reference statistics are used.
- `--remove --applyall` with a date range is the recommended way to reset prior auto-cull results and rerun `cleanoff`.

---

## Key internal functions

| Function | Description |
|----------|-------------|
| `getVelRef(noCull, epsg, wktFile)` | Load velocityStats mean and sigma maps; return dicts keyed by frame number. |
| `getOffsetDirs(noCull, offDirRoot)` | Find valid offset dirs (have azimuth.offsets, no Exclude, have velocity dir). |
| `processOffsets(...)` | Main per-segment worker: read vel, build mask, convert to r/a, write list. |
| `makeSigThresh2D(vel, sigThresh, sigShape)` | Build float32 per-pixel threshold array; apply shape-file scale factors if present. |
| `makeInitialCullMask(vel, mean)` | Mark valid data pixels and compute goodpoints intersection. |
| `cleanMaskEdges(maskb, kernel1)` | Find ragged edge pixels via morphological open. |
| `makeDiffMask(vel, ...)` | Flag outliers and edges; return closed binary mask. |
| `computeBadRA(mask1, vel, geodatFile, dem, offDir)` | Convert bad-pixel XY → lat/lon → r/a via `lltora`. |
| `makeRAMask(offDir, r, a)` | Project r/a into offset grid coords; dilate; return flat int32 index array. |
| `getDem(dems, myDate)` | Select gimp1 (pre-2015) or gimp2 (2015+) DEM for Greenland. |
| `readSigThresh(offDir, sigThresh)` | Return sigThresh, overridden by `sigThresh.override` file if present. |
| `getRegion(region)` | Return region from arg; override from `../region` file if present. |
