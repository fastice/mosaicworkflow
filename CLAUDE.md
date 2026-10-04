# CLAUDE.md — mosaicworkflow

Produces the GrIMP velocity and image mosaics delivered to NSIDC. Orchestrates the C-language binaries in GIT64 using Python threading for parallelism. See the [packages CLAUDE.md](../CLAUDE.md) for the full pipeline context.

## Programs

| Script | Entry point | Role |
|---|---|---|
| `makemosaic` | `makemosaic.py:main()` | Time-series loop; front-end to `setupquarters` |
| `setupquarters` | `setupquarters.py:main()` | Produce one velocity mosaic for a single time period |
| `setupimagemosaic` | `setupimagemosaic.py:main()` | Image (amplitude) mosaic production |
| `makeimagemosaics` | `makeimagemosaics.py:main()` | Time-series loop for image mosaics |
| `simoffsets` | `simoffsets.py:main()` | Simulate SAR offsets (Python wrapper) |
| `prepareTSXrelease` | `prepareTSXrelease.py:main()` | TSX release preparation |
| `makeMaster` | `makeMaster.py:main()` | Assemble the master `inputFile`(s) for `mosaic3d` from several source regions, driven by a YAML config (see [doc](Documents/makeMaster.md)) |

The two most important programs are `makemosaic` → `setupquarters`.

---

## makemosaic

Front-end that loops over a date range, calling `setupquarters` once per period.

### Key arguments

```
makemosaic [--firstdate YYYY-MM-DD] [--lastdate YYYY-MM-DD]
           [--interval {s1cycle,s1-12day,nisarcycle,monthly,quarterly,annual,multiYear}]
           [--inputFile FILE] [--landsatPath DIR]
           [--mosaicMask SHAPEFILE] [--outputMask SHAPEFILE]
           [--baseFlags "..."] [--LSFitType {custom.mask,custom,static,static.mask}]
           [--template mosaic.template.yaml]
           [--mosaicsSetupFile mosaicsSetup.yaml]
           [--noReprocess] [--noLandsat] [--keepFast] [--noLabel] [--check]
```

### Date intervals

- `s1cycle` — Sentinel-1 repeat period (6-day before 2021-12-24, 12-day otherwise)
  with one 13-day window from 2026-06-02 to absorb the 2026 S1 phase shift
  (`mosaicfunc.S1_PHASE_SHIFTS`; 7-day window in the 6-day geomosaic chain; NISAR
  exempt). See [Documents/s1PhaseShift2026.md](Documents/s1PhaseShift2026.md).
- `s1-12day` — 12-day, except one 13-day window (2026-06-08 – 2026-06-20) for the same shift
- `nisarcycle` — 12-day windows on the NISAR repeat cycles (cycle 1 = 2025-09-23,
  `makeimagemosaics.NISAR_CYCLE_EPOCH`, so velocity and image cycle mosaics share windows);
  no phase-shift exceptions, TSX excluded. Pair selection is `setupquarters --overlapSelect`:
  every pair overlapping the cycle, with the true 12-day window passed to mosaic3d so
  `-timeOverlap` weights each pair by its overlap (a true 12-day average). The earlier
  `--refDateSelect` (whole pairs referenced in the cycle, `-date2` + 48 d) defeated that
  weighting; it is kept as an option but no longer used by `nisarcycle`.
- `monthly` — calendar month
- `quarterly` — 3-month season (quarters defined by region in `mosaicsSetup.yaml`)
- `annual` — calendar year (Oct–Sep for Greenland)
- `multiYear` — single product spanning the full date range

### mosaicsSetup.yaml

YAML file with seasonal configuration:
```yaml
landsatPath: /path/to/landsat/listfiles
mosaicMask: /path/to/mask.shp   # or null
intervalMasks:
  monthly: /path/to/mask_YY_SS.shp   # YY=year, SS=season code
  quarterly: ...
  annual: ...
whichMask:           # season code by month index (1-based)
  quarterly: [...]
  monthly: [...]
seasonData:
  quarterlyDates: [0,12,12,3,3,3,6,6,6,9,9,9,12]  # month → start month
  quarterlyMonth: 2
  annualFirstMonth: 12
```

### Landsat list merging

`getLists()` finds `Listfile.*.{fitType}` files in `{landsatPath}/listfiles/`. `createMergedList()` merges all years needed into a single `mergedList.YYYY-YYYY` file passed to `setupquarters`.

---

## setupquarters

Produces a single velocity mosaic. This is the core program.

### Key arguments

```
setupquarters [--template mosaic.template.yaml]
              [--firstdate YYYY-MM-DD] [--lastdate YYYY-MM-DD]
              [--inputFile FILE] [--lsFile FILE]
              [--mosaicMask SHAPEFILE] [--outputMask SHAPEFILE]
              [--dem FILE] [--baseFlags "..."]
              [--noTSX] [--noCull] [--noReprocess] [--noLabel]
              [-timeOverlapOff] [--check]
```

### mosaic.template.yaml

Defines the output grid and processing parameters. Key fields:

```yaml
regionID: greenland          # used in output filenames
xll: -650.0                  # lower-left x (km, polar stereographic)
yll: -3380.0                 # lower-left y (km)
sx: 1500                     # output size in x (pixels)
sy: 2100                     # output size in y
dx: 0.5                      # pixel spacing x (km)
dy: 0.5                      # pixel spacing y
epsg: 3413                   # or null if wktFile used
wktFile: null                # WKT projection file (alternative to epsg)
dem: /path/to/dem            # DEM for simulations
velMap: /path/to/velmap      # Reference velocity map
inputFile: /path/to/inputfile  # SAR data-take input file (can be overridden)
baseFlags: "-vhParams ... "  # flags passed directly to mosaic3d
outputMask: /path/to/mask.shp
regionFile: /path/to/region.yaml  # optional: reads region-specific settings
```

Can include `regionFile:` to pull `dem`, `velMap`, `epsg`, `wktFile` from a separate YAML.

Sections (for domain chunking) are defined by `sectionTemplate()` from `mosaicfunc` — it reads the template and returns lists of `(xdim, ydim)` dictionaries, each defining a sub-domain.

### SAR input file format

The master input file (one data take per line, length > 80 chars triggers parsing):
```
; comment
{phaseFile} {geodatFile} {date} {dT_days} {weight} {crossOrbitFlag}
...
```
- `phaseFile` — path to the range offset file (`.vrt` or binary)
- `geodatFile` — path to geodat (`.geojson` or `.in`)
- `crossOrbitFlag` — `1` for non-TSX (crossing orbit), `0` for TSX (same orbit)

Entries with an `Exclude` file two levels up from the geodat are silently skipped.

### Parallelism model

**Domain chunking**: `mosaicfunc.sectionTemplate()` splits the output domain into a grid of sub-regions (pieces). Each piece runs `mosaic3d` in a separate Python `threading.Thread`.

```python
# Per sector thread call:
command = f'cd {outPath}; mosaic3d -writeBlank -center {flags} '
          f'pieces/{inputFileName} {dem} pieces/{outFileName}'
call(command, shell=True)
```

Maximum concurrent threads: `maxThreads = 24` (hardcoded in `main()`).

The C binaries (`mosaic3d`) are each **single-threaded** — parallelism is entirely from Python threading over sectors. The OpenMP threading added to `mosaic3d` (flag `-ompThreads`) is additive.

### Processing pipeline within setupquarters

```
1. processQArgs()         — read template + flags
2. processInputFileQ()    — parse + filter data takes by date
3. processSectors()       — create per-sector input files + threads
4. runMyThreads()         — run mosaic3d on all sectors (up to 24 concurrent)
5. maskMosaics()          — apply output shapefile mask to each piece (threaded)
6. interpMosiacs()        — run intfloat on each piece/component (threaded)
7. writeQTiffs()          — write per-piece COG GeoTIFFs (threaded)
8. makeQVRTSs()           — gdalbuildvrt to assemble pieces into full-domain VRTs
9. makeShapeOutputs()     — create SAR + Landsat frame shapefiles
10. makeFinalOutput()
    → processFinalVRTs()  — gdal_translate VRT → final COG GeoTIFFs (threaded, 8 at once)
    → processFinalPreview() — HSV-coloured speed browse image + JPG quicklook
    → makeFinalShapes()   — rename shapefiles to NSIDC naming
11. makePremets()         — write NSIDC premet metadata files
12. makeSpatial()         — write .spo corner-point files
```

### Output directory structure

```
Vel-YYYY-MM-DD.YYYY-MM-DD/
├── io/                  # stdout/stderr logs for each piece
├── pieces/              # per-sector mosaic outputs + input files
├── masked/              # output-mask-applied pieces
├── interp/              # intfloat-interpolated pieces + geodats
├── tiff/                # per-piece GeoTIFFs + assembled VRTs
├── shp/                 # output mask shapefile copy
├── meta/                # SAR + Landsat frame shapefiles
└── release/             # Final NSIDC-ready products
    ├── {regionID}_vel_mosaic_{type}_{d1}_{d2}_vx_v{VV}.{sv}.tif
    ├── {regionID}_vel_mosaic_{type}_{d1}_{d2}_vy_v{VV}.{sv}.tif
    ├── {regionID}_vel_mosaic_{type}_{d1}_{d2}_vv_v{VV}.{sv}.tif
    ├── {regionID}_vel_mosaic_{type}_{d1}_{d2}_ex_v{VV}.{sv}.tif
    ├── {regionID}_vel_mosaic_{type}_{d1}_{d2}_ey_v{VV}.{sv}.tif
    ├── {regionID}_vel_mosaic_{type}_{d1}_{d2}_browse_v{VV}.{sv}.tif
    ├── {regionID}_vel_mosaic_{type}_{d1}_{d2}_browse_v{VV}.{sv}.jpg
    ├── *.premet           # NSIDC metadata
    ├── *.spo              # corner-point spatial files
    └── *.shp / .dbf / .shx / .prj   # SAR + Landsat frame shapefiles
```

### Product naming

`{regionID}_vel_mosaic_{type}_{d1}_{d2}_{component}_v{VV}.{sv}.tif`

- `type` — `Monthly`, `Quarterly`, `Annual`, `s1cycle`, `multiYear`
- `d1`, `d2` — `{day}{Mon}{YY}` format (e.g., `01Jan20`)
- `component` — `vx`, `vy`, `vv`, `ex`, `ey`, `browse`
- `VV` — version from `currentVersions` dict (Monthly=5, Quarterly=5, Annual=5, s1cycle=2)
- `sv` — sub-version from `subVersions` dict (all currently 0)

### Version numbers

Hardcoded dicts at module top:
```python
currentVersions = {'s1cycle': 2, 's1-12day': 2, 'Monthly': 5, 'Quarterly': 5,
                   'Annual': 5, 'multiYear': 1}
subVersions = {'s1cycle': 0, 'Monthly': 0, 'Quarterly': 0, 'Annual': 0,
               'multiYear': 0}
```
Update these when releasing a new version.

### TSX handling

TSX (TerraSAR-X) data is treated specially:
- `--noTSX` flag excludes all TSX from the input
- `--noCull` switches TSX to use `*.interp.da` / `*.interp.dr` files instead of culled offsets
- In non-summer, non-single-cycle products, TSX frames whose midDate falls outside the window are screened out
- Annual products weight TSX by 1/count-per-season to prevent summer TSX from dominating

### Landsat date filtering

`processLsList()` filters Landsat data takes by:
- Date overlap with the mosaic period
- At least 60% coverage (`percentCover >= 0.6`)
- Center date within half a window width of the mosaic center

## Image mosaics: geomosaic output options (`makeimagemosaics` / `setupimagemosaic`)

Added 2026-09-23 for NISAR GCOV mosaics. Each option is a template key, and the CLI flag of
the same name overrides it. `makeimagemosaics` forwards the CLI flags to `setupimagemosaic`.
With none of them set, the geomosaic command and the products are exactly as before.

| template key / CLI | geomosaic | effect |
|---|---|---|
| `int16: true` / `--int16` | `-int16` | calibrated sigma0/gamma0 as Int16 dB×100, with scale 0.01 in the file and nodata −3000. Lossless (geomosaic rounds to 0.01 dB) and about 2× smaller. Assembly uses `-vrtnodata -3000` and a predictor in the COG. Read floats with `rioxarray(..., mask_and_scale=True)` or `gdal_translate -unscale` |
| `calOutput: sigma0\|gamma0\|both` / `--calOutput` | `-calOutput` | products made, assembled, premet'd and checked by `checkProduct` [both] |
| `featherDistance: <m>` / `--featherDistance` | `-fl round(m/dx)` | feathering. Implies averaging, because geomosaic refuses `-fl` with orbit priority |
| `gcov: <yaml>` / `--gcov` | `-gcov` | NISAR GCOV HDF5 inputs, alone (`DBpaths: []`) or alongside range/Doppler inputs. Sets the NISAR platform in the premet |
| `rangeSelect: near\|far` / `--rangeSelect` | `-nearRange` / `-farRange` | keep only the near- or far-range look at each pixel [none]. See `geomosaic.md` "Near-range / far-range selection": cross-track seams become deliberate hard boundaries, and the tolerance's real job is keeping repeat passes of the same track |
| `angleTolerance: <deg>` / `--angleTolerance` | `-angleTolerance` | tolerance for the above [1.0]; an error without `rangeSelect` |
| `angleStride: <n>` / `--angleStride` | `-angleStride` | pre-pass cell size in output pixels [10] |
| `rangeSelect: near\|far` / `--rangeSelect` | `-nearRange` / `-farRange` | keep only the near- or far-range look at each pixel [none]. See `geomosaic.md` "Near-range / far-range selection" — cross-track seams become deliberate hard boundaries, and the tolerance's real job is keeping repeat passes of the same track |
| `angleTolerance: <deg>` / `--angleTolerance` | `-angleTolerance` | tolerance for the above [1.0]; an error without `rangeSelect` |
| `angleStride: <n>` / `--angleStride` | `-angleStride` | pre-pass cell size in output pixels [10] |

| `message: <str or list of lines>` / `--message` | — | text burned into the JPG quicklook with the same `setupquarters.addLabel()` credit-block code (template `xLabel`/`yLabel` placement, font auto-fit). Also stored as COG metadata `message` (lines joined with "; "). Without it, nothing is burned into image mosaics, as before |

`int16` and `calOutput` apply only to `--calibrate` runs. With `int16`, the quicklook is made
with `-ot Byte`; otherwise GDAL writes a 12-bit JPEG that PIL and most viewers cannot open.

Velocity browse images (`setupquarters`): a template `message` now replaces the default credit
block, and without one `DEFAULT_CREDIT_LABELS` is used unchanged. `addLabel` falls back through
`CREDIT_FONTS`, because the original DejaVu path does not exist on Ubuntu hosts (helheim/store).
It also now returns the image rather than `None` when `xLabel`/`yLabel` are missing.

**Feathering and tiles.** A tile cannot see data edges beyond its own grid, so feathering
alone would leave steps at tile boundaries: 1.3% of pixels, up to 0.32 dB, in a test. With
`featherDistance` set, `padInputFiles()` grows every tile by one feather length. `cropTiles()`
crops those margins into `pieces/cropped/` before `gdalbuildvrt`; the padded originals stay, so
`--noMosaic` reassembly is safe. Verified: the tiled result is identical to an untiled padded
run, with 0 of 10⁶ pixels different.

GCOV-only template: `DBpaths: []`, `tracks: []`, `trackOrder: []`, plus a `gcov:` yaml. Tile
input lists then have 0 range/Doppler products. With `nisarCycle: true`, cycle N covers
2025-09-23 + 12(N−1) and the 11 days after it; cycle 30 is 2026-09-06 to 09-17. A test setup is
in `/Volumes/insar4/ian/Data/NISAR/Antarctica-GCOV/mosaicTest/workflowTest/`.

## Squint correction (`--useSquint`)

`refreshties.py` and `makeframetie.py` both accept `--useSquint`, which threads through to `tie_script` (i.e. `insarworkflow.tieScript`), which in turn passes `-useSquint` to both `mosaic3d` and `tiepoints -motion`. Default off; can also be set project-wide via `applySquintCorrection: true` in `project.yaml` and overridden per-run by `setupNISARTracks --useSquint`/`--noUseSquint`.

See `mosaicSource/CLAUDE.md` "Squint (residual Doppler) sensitivity" for the C-side implementation and why this is off by default.
