# makeframetie.py

Splits a per-year tie-point plan file into per-frame-range sub-plans, runs `tie_script` on each, then runs `vel_thumbs` to generate velocity thumbnails. Run from the `tiepoints/` directory of a track.

---

## Usage

```
makeframetie.py [options] tie_plan
```

---

## Arguments and options

| Argument / Option | Description |
|-------------------|-------------|
| `tie_plan` | Tie-plan file to process (e.g. `tie_plan2024` or `tie_plan2024-26` for TSX multi-track). |
| `--overWrite` | Pass `--overWrite` to `vel_thumbs` to rerun existing products. |
| `--keepVz` | Pass `--keepVz` to `vel_thumbs` to retain `.vz` and `.vz.geodat` files. |

---

## Behaviour

1. **Detects sensor and frame directories** (`getSensorInfo`) from `project.yaml` (preferred) or `sensor.yaml` two directories up. Reads `sensor`, `framePattern`, and derives frame-range directories from `../velocityStats/*-*/`.
2. For each frame range `X-Y`:
   - Reads the input `tie_plan` and filters segment lines to include only frames in `[X, Y]`.
   - Writes `tie_planYEAR.XdashY` containing only those segments plus a matching `all*dash*` tiefile block.
   - Writes `vel_thumb_plan_all*dash*` from the corresponding `vel_thumb_header_XdashY` file.
   - Runs `tie_script tie_planYEAR.XdashY` to generate tiepoints for that frame range.
   - Runs `vel_thumbs [--overWrite] [--keepVz] vel_thumb_plan_...` to generate velocity thumbnails.

For TSX multi-track, filenames include the track number as a suffix and the header file is `vel_thumb_header_XdashY-TRACK`.

---

## Programs called

| Program | Purpose |
|---------|---------|
| `tie_script` | Runs tiepoint estimation for the per-frame sub-plan |
| `vel_thumbs` | Generates velocity thumbnail maps for each frame range |

---

## Key internal functions

| Function | Description |
|----------|-------------|
| `getSensorInfo(tiePlanFile)` | Detects sensor, track directory, frame ranges, and `framePattern` from `project.yaml`/`sensor.yaml` or path. |
