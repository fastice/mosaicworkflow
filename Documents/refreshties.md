# refreshties.py

Top-level driver that runs `maketies.py` and `makeframetie.py` in parallel across multiple track directories. Typically run from the main project directory (one level above the individual track directories).

---

## Usage

```
refreshties.py [years ...] [options]
```

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `years` | current year | One or more years to process. |
| `-toRun TRACKS` | `['track-26','track-74','track-90','track-112','track-141','track-170']` | Python list of track directories to process, e.g. `"['track-26','track-74']"`. |
| `-tieFiles FILE` | None | YAML file with time-varying tiepoints; passed to `maketies.py`. |
| `-tiesOnly` | False | Run `maketies.py` only; skip `makeframetie.py`. |
| `-velthumbsOnly` | False | Run `makeframetie.py` only; skip `maketies.py`. |
| `-phase` | False | Phase tiepoints mode — sets `-phase` flag on `maketies.py` and implies `-tiesOnly`. |
| `-winter` | False | Use winter tiepoints (passes `-winter` to `maketies.py`). |
| `-noPrompt` | False | Run without confirmation prompt. |
| `-noQuadFit` | False | Pass `--noQuadFit` to `tie_script` (via `maketies.py`) to skip the `-deltaBQ` quadratic baseline correction estimate (`rBaseline.quad`). |
| `--overWrite` | False | Pass `--overWrite` to `makeframetie.py` to rerun existing products. |
| `--keepVz` | False | Pass `--keepVz` to `makeframetie.py` to retain `.vz` and `.vz.geodat` files. |
| `--useSquint` | False | Pass `--useSquint` to `makeframetie.py` (squint heading correction in `mosaic3d` and `tiepoints -motion`). |
| `--tiff` | False | Pass `--tiff` to `makeframetie.py`/`vel_thumbs` so velocity products are written as GeoTIFF (`mosaic3d -GTiff`) regardless of `project.yaml` `velThumbOutput`. Forces tiff on only — it never selects binary. `setupS1Tracks` passes this unconditionally. |

---

## Behaviour

For each track directory in `toRun`, a thread is spawned that:

1. Detects the sensor from `project.yaml` in the project root (preferred), falling back to `sensor.yaml` with a warning, then falling back to path-based detection (`TSX`, `CSK`, `Sentinel1`, `ALOS2`, `NISAR`).
2. Unless `-velthumbsOnly`: runs `maketies.py -run [flags] [years]` inside the track directory.
3. Unless `-tiesOnly`: changes into `tiepoints/` and runs `makeframetie.py [--overWrite] [--keepVz] [--useSquint] [--tiff] tie_planYEAR[suffix]` for each year.

All threads run concurrently (up to 24 at a time via `u.runMyThreads`). stdout and stderr for each track are written to `<trackdir>/stdout` and `<trackdir>/stderr`.

---

## Programs called

| Program | When |
|---------|------|
| `maketies.py` | unless `-velthumbsOnly` |
| `makeframetie.py` | unless `-tiesOnly` or `-phase` |
