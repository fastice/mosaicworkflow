# cleanoffmerge.py

Run inside an individual offset directory (e.g. `track-N/YYYYMMDD_YYYYMMDD_00??/`).  Merges fast and slow offsets when available, applies bad-pixel lists produced by `autoclean`, then cleans scene edges and removes isolated islands/pin-holes via `intfloat`.

---

## Usage

```
cleanoffmerge [--sensor SENSOR]
```

Must be run from the offset pair directory.

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--sensor SENSOR` | `S1` | Sensor name used to look up interpolation parameters: `S1`, `TSX`, or `CSK`. |

---

## Algorithm

### 1. Fast/slow merge

If `fast/azimuth.offsets.noclean.fast` exists, `mergeoff.py` is called to merge the fast offsets with the slow offsets (requires `azimuth.offsets.slow`).

### 2. Load offsets

Reads `azimuth.offsets` (and the companion range file).  Uses a VRT file (`offsets.range-azimuth.vrt`) if present; otherwise uses `.dat` sidecar files.

### 3. Hand-cleaning masks

Applies `badoffsets_poly.idl` and `badoffsets_auto.idl` (IDL save files) if present.

### 4. Autocull bad-pixel lists

Detects mode from which list files exist:

**RA mode** — triggered when `badoffsets_auto.list.dr` and/or `badoffsets_auto.list.da` are present (written by `autoclean --doRA`):

- `badoffsets_auto.list.dr` → zeroed in the **range** offset array only (`rgOff`).
- `badoffsets_auto.list.da` → zeroed in the **azimuth** offset array only (`azOff`).
- Each component is removed independently so a bad range pixel does not discard its azimuth measurement and vice versa.

**XY mode** — triggered when `badoffsets_auto.list` is present (written by `autoclean` in XY mode):

- Applies the index list to both range and azimuth offset arrays simultaneously via `off.removeList`.

### 5. Post-removal cleaning (both modes)

After bad-pixel removal the following steps are applied in order:

1. **Edge clipping** (`fixEdgeEffects`) — clips the first and last 4 azimuth lines of data within each range column when the data stripe starts within 120 lines of the image edge; prevents false matches at scene boundaries.
2. **Ragged-edge removal** (`cleanEdges`) — morphological open (erosion then dilation) on the valid-data mask; pixels that survive erosion but not the open are fringe artefacts and are zeroed.
3. **Island/pin-hole filter** (`islandThreshAndPinHole`) — runs `intfloat -wdist` on a temporary copy of the offsets to fill small holes and remove isolated islands, using sensor-specific thresholds from `sarfunc.sensorDefinitions`.

### 6. Write results

Writes cleaned `azimuth.offsets` and `range.offsets` (and sigma files) back to the current directory.

---

## Key files

| File | Description |
|------|-------------|
| `badoffsets_auto.list` | Flat binary int32 indices (XY mode); applied to both components |
| `badoffsets_auto.list.dr` | Flat binary int32 indices (RA mode); applied to range only |
| `badoffsets_auto.list.da` | Flat binary int32 indices (RA mode); applied to azimuth only |
| `badoffsets_poly.idl` | IDL save file from hand polygon cleaning |
| `badoffsets_auto.idl` | IDL save file from automated IDL cleaning |
| `fast/azimuth.offsets.noclean.fast` | Presence triggers fast/slow merge |
| `offsets.range-azimuth.vrt` | VRT file; used instead of `.dat` sidecars when present |

---

## Notes

- RA vs XY mode is selected automatically from which list files are present; no flag is needed.
- `cleanoffmerge` is normally invoked by `cleanoff` (the per-pair shell script), not directly.
- The `--sensor` flag must match the sensor used to produce the offsets so that the correct interpolation thresholds are applied.

---

## Key internal functions

| Function | Description |
|----------|-------------|
| `cleanoffmergeArgs()` | Parse `--sensor` argument. |
| `islandThreshAndPinHole(off, sensorInfo)` | Fill pin-holes and remove islands via `intfloat`; returns cleaned `off` object. |
| `cleanEdges(off)` | Remove ragged fringe pixels via morphological open on the valid-data mask. |
| `fixEdgeEffects(off)` | Clip leading/trailing scene-edge pixels to suppress false matches at data/no-data transitions. |
