# setupquarters.py

Produces a single GrIMP-format velocity mosaic for one time interval. Typically called repeatedly by `makemosaic.py` but can be run standalone. Reads SAR and optional Landsat inputs, runs `mosaic3d` in parallel over spatial sectors, masks, interpolates, writes tiffs, and assembles NASA release products including premets and spatial files.

---

## Usage

```
setupquarters.py [options]
```

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--firstdate YYYY-MM-DD` | `2015-01-01` | Start of mosaic window. |
| `--lastdate YYYY-MM-DD` | today | End of mosaic window. |
| `--noReprocess` | False | Skip mosaicking; only reformat/repackage existing piece outputs. |
| `--noTSX` | False | Exclude TSX/CSK data from input. |
| `--noLabel` | False | Suppress NASA/ESA credit label on browse image. |
| `--noCull` | False | Use unculled offset files (`.interp.da`/`.interp.dr`) instead of standard culled offsets. |
| `-timeOverlapOff` | False | Disable `-timeOverlap` flag passed to mosaic3d. |
| `--baseFlags STR` | None | Override default mosaic3d flags (from template). |
| `--outputMask FILE` | None | Shapefile used to mask final output (copied into `outDir/shp/`). |
| `--mosaicMask FILE` | None | Mask passed to mosaic3d as `-shelfMask`. |
| `--check` | False | Set up inputs and check masks but do not run mosaic3d. |
| `--dem FILE` | None | Override DEM path from template. |
| `--inputFile FILE` | None | SAR input file (overrides template). |
| `--lsFile FILE` | None | Landsat list file for this mosaic period (overrides template). |
| `--template FILE` | `mosaic.template.yaml` | YAML template defining mosaic geometry and parameters. |

---

## Template file (`mosaic.template.yaml`)

Required keys:

| Key | Description |
|-----|-------------|
| `xll`, `yll` | Lower-left corner in projected coordinates (metres). |
| `sx`, `sy` | Number of pixels in x and y. |
| `dx`, `dy` | Pixel spacing in x and y (metres). |
| `regionID` | Short region identifier used in output filenames (e.g. `GRE`). |
| `dem` | Path to DEM file. |
| `baseFlags` | Default flags string passed to mosaic3d. |
| `epsg` | EPSG code for output projection (or omit if `wktFile` given). |
| `wktFile` | Path to WKT projection file (alternative to `epsg`). |
| `velMap` | Velocity map used by siminsar. |
| `sigmaShape` | Sigma shape file. |
| `xLabel`, `yLabel` | Pixel offsets for credit label in browse image. |

Optional keys (can be overridden by command line):

`inputFile`, `lsFile`, `outputMask`, `mosaicMask`, `baseFlags`, `dem`, `wktFile`, `regionFile`

---

## Output directory structure

Output is written to `Vel-<firstdate>.<lastdate>/`:

```
Vel-YYYY-MM-DD.YYYY-MM-DD/
  io/              mosaic3d stdout/stderr logs per sector
  pieces/          per-sector input files and raw mosaic3d outputs
  masked/          per-sector outputs after shapefile masking
  interp/          per-sector outputs after hole interpolation
  tiff/            per-sector COG tiffs + sector VRTs
  shp/             output mask shapefile copy
  meta/            SAR and Landsat frame shapefiles
  release/         final versioned COG tiffs, browse, premets, spatial files
```

---

## Product types and versioning

Product type is inferred from the interval duration (`d2 - d1`):

| Type | Duration | Notes |
|------|----------|-------|
| `s1cycle` | ≤ 12 days | Sentinel-1 repeat; TSX excluded by `makemosaic.py` |
| `Monthly` | 27–31 days | |
| `Quarterly` | 86–93 days | |
| `Annual` | 363–367 days | Phase window extended ±61/+121 days; TSX weights equalised by season |
| `multiYear` | > 367 days | Single product |

Current version numbers are set in `currentVersions` / `subVersions` dicts at the top of the file.

---

## Input file format (SAR)

Lines with length > 80 and no `;` are data lines. Each line contains (space-separated):

```
<phase_file>  <geodat>  <weight>  <dT_days>  <weight2>  ...
```

- Lines are filtered to those whose mid-date (`firstDate + dT/2`) falls in the mosaic window (with special rules for Annual and s1cycle).
- Lines containing `/Exclude` in their geodat path are skipped.
- For `--noCull`: data lines are replaced with `.interp.da`/`.interp.dr` unculled offsets; fast-tracked areas add a second line pointing to `fast/` subdirectory.
- For Annual products: TSX weights are rescaled to `1/count` per glacier per season to prevent dominance by dense summer acquisitions.

---

## Spatial sectors

The mosaic is broken into a grid of spatial sectors defined by `mosf.sectionTemplate()` reading the template YAML. Each sector runs `mosaic3d` as a separate thread (up to 24 threads).

Per-sector input files are written as `pieces/inputFile.XXX.YYY` with only the data takes that have any temporal overlap with the window.

---

## Freestanding programs called

### `mosaic3d`

The core mosaicker. Called once per spatial sector via `runSubMosaic()`:

```
cd <outDir>; mosaic3d -writeBlank -center <baseFlags>
    -date1 MM-DD-YYYY -date2 MM-DD-YYYY
    [-landSat <lsCulledFile>]
    [-timeOverlap]
    [-shelfMask <mosaicMask>]
    pieces/<inputFile.XXX.YYY> <dem> pieces/<mosaic-XXX.YYY>
```

For Annual products the date window passed to `mosaic3d` is expanded by −61/+121 days to capture phase data from the flanking seasons.

### `intfloat`

Hole interpolation (`interpMosiacs()`). Called per sector per component (`.vx`, `.vy`, `.dT`):

```
cd <outDir>; intfloat -thresh 50 -wdist -nr <sx> -na <sy> -ratThresh 1.
    masked/<pieceName><component> > interp/<pieceName><component>
```

### `gdalbuildvrt`

Assembles per-sector tiffs into a single VRT per component (`makeQVRTSs()`):

```
cd <outDir>/tiff; gdalbuildvrt -te xmin ymin xmax ymax -tr dx dy
    -vrtnodata <nodata> <outDir><suffix>.vrt *<suffix>.tif
```

### `gdal_translate`

Converts sector VRTs to final versioned Cloud-Optimised GeoTIFFs (`processFinalVRTs()`):

```
gdal_translate -of COG -co COMPRESS=DEFLATE -co RESAMPLING=AVERAGE
    -co OVERVIEWS=IGNORE_EXISTING -co GEOTIFF_VERSION=1.1
    -co BIGTIFF=NO -stats <vrt> ../release/<finalName>
```

### `gdal_translate` (browse JPG)

Downsamples browse COG tiff to 500 m JPEG quicklook (`writeJPG()`):

```
gdal_translate -co "QUALITY=99" -scale -of JPEG -r average
    -tr 500 500 <browseTiff> <jpgFile>
```

### `makeimageshapefile.py`

Creates a shapefile of SAR frame footprints (`makeShapeOutputs()`):

```
makeimageshapefile.py -velocity
    -firstdate YYYY:MM:DD -lastdate YYYY:MM:DD
    <inputFile> <outDir>/meta/SAR.<firstdate>.<lastdate>
```

### `pathrowshp.py`

Creates a shapefile of Landsat path/row footprints (`makeShapeOutputs()`):

```
pathrowshp.py -firstdate=YYYY:MM:DD -lastdate=YYYY:MM:DD
    <lsFile> <outDir>/meta/LS.<firstdate>.<lastdate>
```

---

## Processing flow

```
processQArgs()
  └── reads template YAML, resolves regionDefs

mkQDirs()           — create output subdirectories, copy output mask shapefile

processLsList()     — filter Landsat list by date → lsData.<d1>.<d2>

processInputFileQ() — parse SAR input file, filter by date, apply TSX rules

processSectors()
  └── writeInputFileQ() per sector  — write filtered per-sector input files
  └── runSubMosaic()   per sector  — thread: run mosaic3d

(wait for all mosaic3d threads)

maskMosaics()       — apply output shapefile mask to pieces (threaded)
interpMosiacs()     — intfloat hole fill on masked pieces (threaded)
writeQTiffs()       — write per-piece COG tiffs (threaded)
makeQVRTSs()        — gdalbuildvrt per component across all pieces

makeShapeOutputs()
  └── makeimageshapefile.py  — SAR frame shapefile
  └── pathrowshp.py          — Landsat frame shapefile

makeFinalOutput()
  └── processFinalVRTs()     — gdal_translate → release/ COG tiffs
  └── processFinalPreview()  — HSV speed browse COG tiff + JPEG
  └── makeFinalShapes()      — copy/rename meta shapefiles into release/

SatTypes()          — detect which platforms (S1A, TSX, LS8…) are present

populatePremet()    — build premet metadata records
makePremets()       — write .premet files in release/
makeSpatial()       — write .spo corner-point spatial files in release/
makeShapePremet()   — write .premet files for shapefiles
```

---

## Output file naming

Final release files follow the convention:

```
<regionID>_vel_mosaic_<Type>_<d1>_<d2>_<component>_v<VV>.<sv>.tif
```

Example: `GRE_vel_mosaic_Monthly_01Jan25_31Jan25_vv_v05.0.tif`

Components: `vx`, `vy`, `vv` (speed), `ex`, `ey`, `dT` (if `-timeOverlap`), `browse`.

---

## Key internal functions

| Function | Description |
|----------|-------------|
| `processQArgs()` | Parse args, load template, resolve region defs. |
| `processTemplate(args)` | Merge template YAML with command-line overrides; compute grid corners. |
| `processFlags(args, template)` | Parse dates, apply `-timeOverlap` to baseFlags. |
| `processInputFileQ(filename, noTSX)` | Read SAR input file; filter TSX if requested; return list of `[date, line]`. |
| `writeInputFileQ(...)` | Filter data takes by date for one sector; write sector input file. |
| `changeToNoCull(line)` | Replace culled offset filenames with `.interp.da`/`.interp.dr` unculled versions. |
| `computeTSXSeasonCount(...)` | Count TSX acquisitions per glacier per season for weighting. |
| `setLineWeight(dataTake, seasonCount)` | Set TSX weight field to `1/count`. |
| `prodType(d1, d2)` | Infer product type from duration. |
| `processSectors(...)` | Divide mosaic into sectors; set up threads; return thread list and geodat list. |
| `runSubMosaic(...)` | Thread target: run `mosaic3d` for one sector. |
| `maskMosaics(...)` | Apply shapefile mask to all sector velocity and error images (threaded). |
| `interpMosiacs(...)` | Run `intfloat` on all sectors and components (threaded). |
| `writeQTiffs(...)` | Write COG tiffs for all sectors (threaded). |
| `makeQVRTSs(...)` | Build VRTs spanning all sector tiffs per component. |
| `makeFinalOutput(...)` | Assemble release: final COGs, browse, shapes. |
| `processFinalVRTs(...)` | `gdal_translate` VRTs → versioned release COGs (threaded). |
| `processFinalPreview(...)` | Generate HSV-coloured browse COG and JPEG quicklook. |
| `populatePremet(...)` | Build list of sarfunc `premet` objects with platform/date/version metadata. |
| `makePremets(...)` | Write `.premet` files to release/. |
| `makeSpatial(...)` | Write `.spo` corner-coordinate spatial files. |
| `finalName(vrt, flags, pType)` | Construct versioned release filename from interval and component. |
