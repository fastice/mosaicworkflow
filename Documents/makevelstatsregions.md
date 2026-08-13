# makevelstatsregions.py

Run from a `track-N/` directory after `velThumbs` has produced per-segment velocity mosaics.  For each `velocityStats/F1-F2` frame-range sub-directory, finds all `*_*/velocity/` outputs whose frame number falls within `[F1, F2]`, unions their geographic bounding boxes, and writes (or updates) a `tiepoints/vel_thumb_header_F1dashF2` file with the resolved `resolution` line.

Supports both mosaic3d output formats:
- **Tiff mode** — reads extent from `mosaicOffsets.vrt` (GDAL geotransform, metres → km)
- **Binary mode** — reads extent from `mosaicOffsets.vx.geodat` (fallback)

The VRT is tried first; if absent the geodat is used.

---

## Usage

```
makevelstatsregions [--dem FILE] [--TSX]
```

Must be run from the `track-N/` directory.  `project.yaml` is read from the parent directory.

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--dem FILE` | auto | DEM path written into new header files. Auto-detected from `regionFile`/`region` in `project.yaml` if not given. |
| `--TSX` | off | TSX mode: set `baseDir` to `../` (one level up) and disable the 10 km geographic padding added to each axis. |

---

## Algorithm

1. Read `../project.yaml` (or legacy `../sensor.yaml`) to get `framePattern` and the DEM path from the region file.
2. Glob `*_*/velocity` to find all per-segment velocity directories.
3. Glob `velocityStats/*-*` to enumerate frame-range directories.
4. For each frame-range `F1-F2`:
   - Parse the frame number from each velocity directory name using `framePattern`.
   - Skip any directory that contains an `Exclude` file.
   - For all velocity directories with frame ∈ `[F1, F2]`, read `mosaicOffsets.vx.geodat` and call `geo.boundsInKm()` to get the geographic bounding box.
   - Union the bounding boxes across all matching directories.
   - Add 10 km padding on every side (disabled by `--TSX`).
   - Derive pixel size from the last geodat read.
   - Call `writeVelThumb()` to write or update the header file.

---

## Output: `tiepoints/vel_thumb_header_F1dashF2`

Each header file drives a `velThumbs` run for the corresponding frame range.

**New file** (no prior header exists):

```
DEM       = <dem path>
track_root= <abs path to track dir>
tie_dir   = <abs path to track dir>/tiepoints
vel_dir = <abs path to track dir>/velocity
extra_flags ="-no3d"
resolution = "<x0> <y0> <xs> <ys> <dx> <dy>"
```

**Existing file** (header already present): all lines except `resolution` are preserved verbatim; the `resolution` line is replaced with the newly computed values.

### Resolution fields

| Field | Meaning |
|-------|---------|
| `x0` | West edge of bounding box (km), rounded, minus padding |
| `y0` | South edge of bounding box (km), rounded, minus padding |
| `xs` | Width in pixels |
| `ys` | Height in pixels |
| `dx` | Pixel width in km |
| `dy` | Pixel height in km |

---

## Integration with setupNISARTracks

`setupNISARTracks --runVelstatsregions` calls `makevelstatsregions` (no arguments) in each `track-N/` directory via `subprocess.run`.

## Notes

- An `Exclude` file at the top level of a segment directory (`*_*/Exclude`) causes that segment to be skipped for all frame ranges.
- `framePattern` from `project.yaml` (e.g. `'00??'`) controls which characters of the directory suffix are the frame number.  The prefix characters before the first `?` are stripped before converting to an integer.
