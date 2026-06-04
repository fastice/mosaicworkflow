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

The two most important programs are `makemosaic` → `setupquarters`.

---

## makemosaic

Front-end that loops over a date range, calling `setupquarters` once per period.

### Key arguments

```
makemosaic [--firstdate YYYY-MM-DD] [--lastdate YYYY-MM-DD]
           [--interval {s1cycle,s1-12day,monthly,quarterly,annual,multiYear}]
           [--inputFile FILE] [--landsatPath DIR]
           [--mosaicMask SHAPEFILE] [--outputMask SHAPEFILE]
           [--baseFlags "..."] [--LSFitType {custom.mask,custom,static,static.mask}]
           [--template mosaic.template.yaml]
           [--mosaicsSetupFile mosaicsSetup.yaml]
           [--noReprocess] [--noLandsat] [--keepFast] [--noLabel] [--check]
```

### Date intervals

- `s1cycle` — Sentinel-1 repeat period (6-day before 2021-12-24, 12-day otherwise)
- `s1-12day` — always 12-day
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
