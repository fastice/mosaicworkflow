# maketies.py

Builds per-year tie-point plan files (`tie_planYEAR`) for a single track directory and writes a `refreshTies` csh script that chains them together. Run from inside a track directory (e.g. `track-26/`).

---

## Usage

```
maketies.py [years ...] [options]
```

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `years` | current year through present | One or more years to process. |
| `-run` | False | Execute the generated `refreshTies` script immediately after building it. |
| `-winter` | False | Use winter tiepoints (passed to `setuptopstie.py` as `-winter`). |
| `-phase` | False | Generate phase tiepoint plans (passed to `setuptopstie.py` as `-phase`; only active for Sentinel1/ALOS2). |
| `-tieFiles FILE` | None | YAML file with time-varying tiepoints; passed to `setuptopstie.py`. |

---

## Behaviour

1. **Detects sensor and paths** (`getSensorTrackInfo`) from `../sensor.yaml` or the working directory path (`TSX`, `CSK`, `Sentinel1`, `ALOS2`, `NISARTest`). Derives `tieDir` from the path accordingly.
2. Creates `tieDir/old/` if it does not exist.
3. For each year:
   - Calls `setuptopstie.py [flags] YEAR` via `check_output` to generate `tie_planYEAR[-suffix]` in the current directory.
   - Moves any existing `tie_planYEAR` to `tieDir/old/` and renames the new one into `tieDir/`.
   - If `setuptopstie.py` returned a non-zero image count, appends `tie_script tie_planYEAR` to `refreshTies`.
4. Calls `makeAll` to produce `tie_planAll[-suffix]` covering all years and appends it to `refreshTies`.
5. If a `tie_planSpecial[-suffix]` file exists, appends it to `refreshTies`.
6. If `-run` was given, executes `refreshTies` via `csh` from `tieDir`.

For multi-track glaciers (e.g. Jak, Helheim) a track-number suffix is appended to filenames to keep plans separate.

---

## Programs called

| Program | Purpose |
|---------|---------|
| `setuptopstie.py` | Generates the per-year `tie_plan` from SAR pair data |
| `tie_script` | Called from the generated `refreshTies` csh script (if `-run`) |
