# mergeoff

Merges regular (slow) speckle-tracked offsets with fast-tracked offsets for a SAR image pair. Where both exist, takes a weighted average using sensor-specific weights (`regW`, `fastW` from the sensor definition). Where only fast offsets exist, fills in from the fast data. Also applies an optional solid-earth correction (`offsets.SECorrection`) and handles IDL hand-cleaning flags (`fast/badoffsets.idl`). Usually called as part of the `cleanoffmerge` pipeline.

Run from within the pair processing directory.

---

## Usage

```
mergeoff.py sensor
```

| Argument | Description |
|----------|-------------|
| `sensor` | Sensor name: `S1`, `CSK`, or `TSX` |

---

## Expected inputs

| File | Description |
|------|-------------|
| `azimuth.offsets.slow` / `offsets.slow.vrt` | Regular (slow) tracked offsets |
| `fast/azimuth.offsets.noclean.fast` / `fast/offsets.noclean.fast.vrt` | Unfiltered fast offsets |
| `fast/azimuth.offsets.fast` / `fast/offsets.fast.vrt` | Fast offsets (post-clean) |
| `offsets.SECorrection` + `.vrt` | Optional solid-earth range correction |
| `fast/badoffsets.idl` | Optional IDL hand-edit bad-pixel list |

---

## Output

`azimuth.offsets` — merged azimuth offset file (and associated range/sigma files).

---

## Part of the mosaicworkflow package.
