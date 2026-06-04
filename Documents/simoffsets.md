# simoffsets.py

Uses a velocity map and a DEM to simulate azimuth and range offset fields for use as the initial guess in SAR feature tracking (strack/strackw). Run in the pair directory for the desired image pair.

---

## Usage

```
simoffsets.py [options]
```

Must be run in the pair processing directory. The second image is located automatically by following the `*.slc` symlink that points to `../`.

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `-azOffsets FILE` | `offsets.da` | Root name for output azimuth offset file. |
| `-offsetsDat FILE` | `offsets.dat` | `.dat` or `.vrt` file defining output grid geometry. |
| `-syncDat` | False | Force re-run of `siminsar`; use `offsetsDat` to define output grid size. |
| `-noVel` | False | Compute offsets assuming zero velocity everywhere (geometry only). |
| `-fastMask` | False | Use the fast-ice mask instead of the standard ice mask. |
| `-LSB` | False | Write output in LSB byte order (default MSB). |
| `-velMap FILE` | region default | Velocity map root name (`.vx`/`.vy` suffixes appended) or `.vrt`/`.tif` path. |
| `-dem FILE` | region default | DEM file. |
| `-region NAME` | auto | Region name (`greenland` or `antarctica`). Auto-detected from geodat latitude if not given. |
| `-regionFile FILE` | None | YAML file with region-specific paths (velMap, DEM, mask etc.). |
| `-secondDir DIR` | `../secondslcdir` | Directory of the secondary image. Auto-detected from `*.slc` symlink. |
| `-geodatFile FILE` | `geodat10x2.in` | GrIMP `.in` or `.geojson` geodat file defining reference image geometry. |
| `-secondGeodatFile FILE` | `secondDir/geodatFile` | Geodat file for the secondary image. |
| `-maskInputFile FILE` | region default | Flat byte mask file (ice/fast-ice). |

---

## Algorithm

### 1. Grid setup (`runSim` → `siminsar`)

If any of `offsets.lat`, `offsets.lon`, `offsets.mask`, or `offsets.dat` are missing, or if `-syncDat` is set, `siminsar` is run to compute the lat/lon/mask grid:

- Without `-syncDat`: an `offsets.dat` file is written with a standard grid derived from the single-look image size (`r0=180`, `a0=180`, `dr=24`, `da=18`).
- With `-syncDat`: the existing `offsetsDat` file defines the grid.

```
siminsar [-LSB] -center -toLL <offsetsDat> -mask -xyDEM <dem>
    <maskFile> <geodatFile> <offsetsRoot>
```

### 2. Static (geometry) offsets (`computeStaticOffsets`)

For every output pixel:

1. Get lat/lon from the reference grid.
2. Project lat/lon into range/azimuth of the **secondary** image via `lltora` (calls the `lltora` binary).
3. Compute `dr0 = r_secondary − r_reference`, `da0 = a_secondary − a_reference`.
4. Points that fail to project (`r < 0`) are set to fill value `−2×10⁹`.
5. Fit a linear plane (`const + r·c1 + a·c2`) to each offset component.
6. Replace gross outliers (|residual| > 2 pixels) with the plane value (`fixBad`).

The polynomial coefficients are written to `<root>.poly`.

### 3. Velocity offsets (`computeVelocityRA`)

1. Read the velocity map (vx, vy in polar-stereographic metres/year).
2. Compute the image heading and the angle from the map x-axis to the look direction.
3. Bilinearly interpolate vx, vy to each output pixel location.
4. Rotate the horizontal velocity vector into the radar (range, azimuth) frame using the heading angle.
5. Flag fast-moving areas (speed > 200 m/yr within the ice mask) with bit 8 in the mask.
6. For left-looking geometry, negate the range component (`vr *= -1`).

### 4. Combining static and velocity offsets

Range:
```
dr = dr0 + vr × groundToSlant × deltaT/365
```

Azimuth:
```
da = da0 + va × deltaT/365 / slpA
```

where:
- `groundToSlant` — conversion factor from ground-range to slant-range pixels, computed from look angle geometry (`groundToSlantRangeResolution`).
- `deltaT` — days between reference and secondary acquisition dates (from geodat).
- `slpA` — single-look azimuth pixel size (m).

Missing pixels (`dr0 < −2×10⁸`) are set to `−2×10⁹` in both components.

---

## Freestanding programs called

### `siminsar`

Computes lat, lon, mask, and optionally a simulated phase or offsets grid for a SAR image pair. Called by `runSim()`:

```
siminsar [-LSB] -center -toLL <offsetsDat> -mask -xyDEM <dem>
    <maskFile> <geodatFile> <offsetsRoot>
```

Produces: `<offsetsRoot>.lat`, `<offsetsRoot>.lon`, `<offsetsRoot>.mask`, `<offsetsRoot>.simdat`.

### `lltora`

Projects lat/lon arrays into range/azimuth pixel coordinates for a given geodat and DEM. Called by `LLtoRA()`:

```
lltora <geodatFile> <dem> temp.<pid>.ll temp.<pid>.ra
```

Input/output via temporary files `temp.<pid>.ll` and `temp.<pid>.ra` (removed on completion).

---

## Output files

All outputs use root name derived from `-azOffsets` (default `offsets`):

| File | Description |
|------|-------------|
| `offsets.da` | Azimuth offset field (binary, MSB float) |
| `offsets.dr` | Range offset field (binary, MSB float) |
| `offsets.da.dat` | Geodat for azimuth offset array |
| `offsets.dr.dat` | Geodat for range offset array |
| `offsets.vrt` | Two-band VRT combining `.dr` (band 1) and `.da` (band 2), with `deltaT` metadata |
| `offsets.poly` | Linear polynomial coefficients for range and azimuth static offsets |
| `offsets.lat` | Latitude array (from siminsar, if regenerated) |
| `offsets.lon` | Longitude array (from siminsar, if regenerated) |
| `offsets.mask` | Ice mask array (from siminsar, if regenerated) |
| `offsets.simdat` | Simulation metadata (from siminsar, if regenerated) |
| `fail.simoffsets` | Created at start, removed on successful completion; presence indicates failure |

---

## Failure flag

`fail.simoffsets` is created at startup and deleted on clean exit. Its presence indicates the program did not complete successfully.

---

## Key internal functions

| Function | Description |
|----------|-------------|
| `simOffsetsProcessArgs1(fp)` | Parse args with `argparse`; resolve region and defaults; check input file existence. |
| `resolveRegion(args)` | Determine region from geodat latitude or `-region`; apply velMap/DEM/mask overrides. |
| `checkFiles(myRegion, secondDir, args, fp)` | Verify velMap, DEM, secondDir, geodat, and mask files all exist. |
| `runSim(geodatFile, offsetsDat, dem, maskInputFile, syncDat, byteOrder)` | Optionally write `offsets.dat`; run `siminsar` to produce lat/lon/mask grid. |
| `computeStaticOffsets(offsets1, secondGeoDat, dem)` | Compute geometry-only offsets via `lltora`; fit and clean outliers. |
| `computePoly(dr0, da0, r1, a1)` | Fit a linear plane to range and azimuth offset fields. |
| `computePlane(x, y, z)` | Least-squares fit of `const + x·c1 + y·c2` to valid points. |
| `fixBad(d, coeff, r1, a1)` | Replace outliers (residual > 2 px) with plane value. |
| `computeVelocityRA(velMap, offsets1, srsInfo)` | Interpolate velocity map; rotate to radar frame; flag fast areas. |
| `groundToSlantRangeResolution(offsets)` | Convert ground-range to slant-range pixels using look angle geometry. |
| `LLtoRA(lat, lon, geodatFile, dem)` | Call `lltora` binary to project lat/lon → range/azimuth; return 2D arrays. |
| `writeOffPoly(fileName, coeffR, coeffA)` | Write polynomial coefficients to text file. |
