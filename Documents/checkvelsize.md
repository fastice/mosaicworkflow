# checkvelsize

Checks that the per-frame velocity mosaics in a track share a common grid.

`setupS1Tracks --runVelstatsregions` finds a common bounding box per velocityStats range, rebuilds `velocity/` on it, and resizes `velocity_nocull/` to match. Confirming that used to be a matter of running `ls track-*/*_*/velocity/*.vx` and checking the byte sizes were all equal. That stopped working once the products became GeoTIFFs — compressed, so equal grids give unequal file sizes.

This reads the grid from the sidecar each format already carries, so no marker files are needed and it works on products that already exist, including a project part-way through the tiff migration:

| Output mode | Grid read from |
|---|---|
| tiff (`mosaic3d -GTiff`) | `mosaicOffsets.vrt` — `rasterXSize`/`rasterYSize` + `GeoTransform` |
| binary | `mosaicOffsets.<band>.geodat` — size, pixel size, origin |

Both go through `utilities.geodat`, so the two origin conventions (VRT is an upper-left pixel corner in metres, `.geodat` is a lower-left pixel centre in km) are normalised to the same lower-left pixel centre in km before comparison.

Run from the project root (the directory holding `track-*/`) or inside a single track directory.

---

## Usage

```
checkvelsize [options]
```

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--tracks track-N [...]` | all `track-*` | Restrict to these tracks |
| `--velDirs DIR [...]` | `velocity velocity_nocull` | Per-frame product dirs to check |
| `--long` | off | List every frame, not just the mismatches |
| `--missingIsError` | off | Also exit non-zero when a frame has no product in one of the `velDirs` |

---

## Output

One line per velocityStats range, giving the grid shared by most of its products, followed by any frame that disagrees:

```
track-119  350-380     229 frames   3415 x 8770  200m  (  -308.0, -3458.0)  MISMATCH (65)
    61791_362  velocity_nocull  3265 x 8620  200m  (  -298.0, -3448.0)  <- differs
    9070_352   velocity_nocull  (no product)

1 track(s), 2 range(s), 242 frame(s): 65 grid mismatch(es), 170 missing product(s)
```

Frames are grouped by the `velocityStats/<first>-<last>` burst ranges, since that is the unit `makevelstatsregions` computes a common box over. A frame whose burst number falls in no range is grouped under `unassigned`. Frames with no product in any `velDir` are skipped entirely — they were never mosaicked, which is a different problem.

Exit status is 1 if any grid disagrees (and, with `--missingIsError`, if any product is missing), so it can gate a rerun:

```
checkvelsize || setupS1Tracks --runVelstatsregions
```

---

## Examples

Whole project, summary only:
```
checkvelsize
```

One track, every frame listed:
```
checkvelsize --tracks track-119 --long
```

Culled velocity only, ignoring `velocity_nocull`:
```
checkvelsize --velDirs velocity
```

---

## Part of the mosaicworkflow package.
