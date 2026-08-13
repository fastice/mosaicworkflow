# initRAReference.py

Builds the vr/va reference basemap used by `velocityStats --doRA` to initialise its outlier-rejection reference. The basemap is a synthetic range/azimuth velocity field derived from geometry only (zero perpendicular baseline), stored under `<frameRangeDir>/sims/mosaicOffsets.vr` and `mosaicOffsets.va`.

Run from a `velocityStats/00n0-00n9/` directory, or supply explicit paths.

---

## Why a per-frame-range basemap is needed

`velocityStats` in XY mode uses a single global vx/vy reference map. Range/azimuth (RA) components cannot share a global reference because the look direction varies between satellite tracks. Each `00n0-00n9` group of frames shares a common track geometry, so a single basemap per group is sufficient.

---

## Usage

```
initRAReference [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--frameRangeDir DIR` | `.` | Frame-range directory (e.g. `0040-0049`). |
| `--trackDir DIR` | `../..` | Track directory containing `*_00??` frame subdirs. |
| `--region NAME` | auto | `greenland` or `antarctica`. Auto-detected from `project.yaml` if not given. |
| `--dem FILE` | region default | DEM file path. |
| `--deltaT N` | `12` | Repeat-pass interval in days used for the zero-baseline secondary geodat. |
| `--ompThreads N` | `4` | OpenMP threads passed to `simoffsets`. |
| `--initializeReference` | off | Rebuild even if `sims/mosaicOffsets.vr` already exists. |

---

## Algorithm

### 1. Collect frame geodats

Scans `trackDir` for `*_00n?` subdirectories whose four-digit frame number falls within the requested range. Picks the primary geodat geojson from each (`geodat*.geojson` excluding any file with `secondary` in its name).

### 2. Build super-geodat (`buildSuperGeodat`)

Creates a single synthetic NISAR-format geodat that covers the full azimuth extent of all frames, plus ±10 km padding, and extends the range swath by ±5 km on each side.

- Selects the frame with the largest `MLAzimuthSize` as the structural base.
- Computes the combined azimuth time span: earliest `NominalTime` → latest `NominalTime + MLAzimuthSize × pixelPeriod`.
- Adds `10 000 m / mlPixA × pixPeriod` seconds of padding at each end.
- `MLAzimuthSize` = pixels needed to span the padded time range.
- `MLNearRange` −= 5 000 m; `MLFarRange` += 5 000 m; `MLRangeSize` adjusted accordingly.
- `MLCenterRange` = (`MLNearRange` + `MLFarRange`) / 2.
- Shifts `NominalTime` and `CorrectedTime` to the new padded start.
- Merges state vectors from all frame geodats, sorts by time-of-day, deduplicates.
- Writes to `sims/super.geojson`.

### 3. Build secondary geodat (`buildSecondaryGeodat`)

Deep-copies the super-geodat and shifts all times by `deltaT` days:
- `NominalTime`, `CorrectedTime`, `TimeOfFirstStateVector` all advance by `deltaT × 86 400` s.
- State vector positions are propagated with linear extrapolation: `pos_new = pos + vel × Δt`.

Because primary and secondary share the same orbital geometry, the perpendicular baseline is zero and all offset signal is pure velocity.

Writes to `sims/super.secondary.geojson`.

### 4. Run simoffsets (`runSimoffsets`)

```
simoffsets -region <region> -syncDat \
    -geodatFile super.geojson \
    -secondGeodatFile super.secondary.geojson \
    -offsetsDat offsets.dat \
    -azOffsets offsets.da \
    -dem <dem> --ompThreads <N>
```

Produces: `offsets.da`, `offsets.dr`, `offsets.dat`, `offsets.lat`, `offsets.lon`, `offsets.mask`.

### 5. Write dummy param files (`createDummyParams`)

Because the baseline is zero, all polynomial correction terms are zero:

| File | Contents |
|------|----------|
| `az.est` | Header with InSAR parameters (H, Re, RNear, Rc, nlr, nla, slpR, wl) then `0.0 0.0 0.0 0.0` |
| `rBaseline` | Same header then `0.0 0.0 0.0 0.0` |
| `baseline.26x16` | InSAR parameter header + 2 lines of `0.0 0.0 0.0 0.0` (first and second flattening baseline) |

Header values are read from the super-geodat properties.

### 6. Write mosaic3d inputFile (`createInputFile`)

```
0 0 0 0 0.2 0.2
;
1
;
nophase <abs>/super.geojson <abs>/baseline.26x16 <deltaT> 1. \
    <abs>/offsets.da <abs>/az.est <abs>/offsets.dr <abs>/rBaseline
```

All paths are absolute so mosaic3d can be run from inside `sims/`.

### 7. Run mosaic3d (`runMosaic3d`)

```
mosaic3d -noVh -no3d -rOffsets -SVConst -outputRA -center \
    inputFile <dem> mosaicOffsets
```

Produces: `sims/mosaicOffsets.vr`, `sims/mosaicOffsets.va`, `sims/mosaicOffsets.vr.geodat`.

---

## Output files

All outputs are placed in `<frameRangeDir>/sims/`:

| File | Description |
|------|-------------|
| `super.geojson` | Synthetic primary geodat spanning all frames + padding |
| `super.secondary.geojson` | Secondary geodat (super + deltaT days) |
| `offsets.da`, `offsets.dr` | Synthetic azimuth/range offsets from simoffsets |
| `offsets.dat`, `offsets.lat`, `offsets.lon`, `offsets.mask` | Grid geometry and mask |
| `az.est`, `rBaseline`, `baseline.26x16` | Zero-polynomial InSAR param files |
| `inputFile` | mosaic3d input list |
| `mosaicOffsets.vr` | Range velocity reference basemap (big-endian float32) |
| `mosaicOffsets.va` | Azimuth velocity reference basemap (big-endian float32) |
| `mosaicOffsets.vr.geodat` | Geodat for the basemap grid |

---

## Integration with velocityStats

`velocityStats --doRA` calls `initRAReference` automatically when `--initializeReference` is also given. Without `--initializeReference`, it expects `sims/mosaicOffsets.vr` to already exist and aborts with a clear error if it does not.

The basemap is loaded once per frame-range directory and used as the initial reference velocity field. Subsequent passes use the statistical mean of accepted maps, exactly as in XY mode.

---

## Key internal functions

| Function | Description |
|----------|-------------|
| `buildSuperGeodat(geodatFiles, simsDir)` | Merge frame geodats into one super-geodat spanning the full azimuth + range extent. |
| `buildSecondaryGeodat(superGeojson, simsDir, deltaT)` | Clone super-geodat and shift times by `deltaT` days; propagate SV positions. |
| `createDummyParams(simsDir, superGeojson)` | Write zero-polynomial `az.est`, `rBaseline`, `baseline.26x16`. |
| `runSimoffsets(simsDir, region, dem, ompThreads)` | Call `simoffsets` to produce synthetic offset files. |
| `createInputFile(simsDir, deltaT)` | Write `sims/inputFile` for mosaic3d. |
| `runMosaic3d(simsDir, dem)` | Call `mosaic3d -outputRA` to produce `mosaicOffsets.vr/va`. |
| `initRAReference(frameRangeDir, trackDir, region, dem, deltaT, ompThreads)` | Top-level orchestrator: calls all of the above in sequence. |
| `_parseNominalTime(timeStr)` | Parse `'HH:MM:SS.ffffff'` → datetime (day fixed to 2000-01-01). |
| `_formatNominalTime(dt)` | Format datetime → `'HH:MM:SS.ffffff'`. |
| `_parseCorrectedTime(timeStr)` | Parse `'HH MM SS.SSSSSSS'` → datetime. |
| `_formatCorrectedTime(dt)` | Format datetime → `'HH MM SS.SSSSSSS'`. |
