# setuptopstie.py

Generates a `tie_planYEAR` file for a single year by scanning SAR pair directories, grouping pairs by temporal baseline (nDays), and writing tiepoint plan entries in the format expected by `tie_script`. Run from inside a track directory (e.g. `track-26/`).

---

## Usage

```
setuptopstie.py [options] YEAR
```

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `year` | required | Year to generate the tie plan for (2008–2035). |
| `-winter` | False | Replace the standard `extraties` path with the winter tiepoints file for `YEAR`. |
| `-phase` | False | Generate a phase-mode tie plan (see [Phase mode](#phase-mode) below). |
| `-tieFiles FILE` | None | YAML file with time-varying tiepoints; used to generate multiple `extraties` blocks with date-bounded intervals. |
| `-tie_plan_suffix SUFFIX` | `` | Suffix appended to the output filename (`tie_planYEAR-SUFFIX`). Used for multi-track glaciers. |

---

## Behaviour

1. **Detects sensor and paths** from `../project.yaml` (preferred) or `../sensor.yaml`. Falls back to path-based detection (`TSX`, `CSK`, `ALOS2`, `NISARTest`, `Sentinel`).
2. **Reads the tie-plan header** (`tiepoints/tie_plan_header` or a year-specific variant) for `extraties`, `DEM`, and other global settings.
3. **Reads any per-track flags** from a `flagsNNNN[-YEAR]` file in the tiepoints directory.
4. **Scans pair directories** (`*_FRAME` pattern) for the year. For each directory, reads pair metadata from a `.pairinfo` file if present, otherwise parses a `strackin.*` file and writes the `.pairinfo` for future runs.
5. **Groups pairs by nDays** and sorts within each group by date.
6. **Writes `tie_planYEAR`**:
   - Copies the header block.
   - Writes `extraties` / bedrock tiefile blocks (one per tiepoint interval if `-tieFiles` given).
   - For each nDays group, writes `use`, `tiefile`, and segment lines (see [Phase mode](#phase-mode) for prefix details).
   - Writes a final `all YEAR` tiefile that collects all orbit groups.
7. Prints the total segment count to stdout (used by `maketies.py`).

---

## Phase mode

When `-phase` is given, `tieCodes()` returns:

| Entry type | Prefix | Suffix |
|------------|--------|--------|
| `use` line (segment processed but not tiepointed) | `pps` | `uw` |
| `tiefile` segment line | `pn` | `uw` |

In normal (speckle) mode both prefixes are blank and suffixes are empty.

`pps` tells `tie_script` to run `process_phase` (baseline estimation + flat-earth removal) for the segment but **not** add it to `tie_segs` for the mosaic run. `pn` adds the segment with phase available but speckle tracking disabled.

---

## Key internal functions

| Function | Description |
|----------|-------------|
| `getMySensorAndTrack()` | Reads `project.yaml`/`sensor.yaml` or inspects the path for sensor, track, region, and root directory. |
| `getSensorInfoAndHeader(sensor, year, sensor_yaml)` | Returns the header file path, `maxDays`, and `framePattern` for the detected sensor. |
| `getPairs(year, phaseTies, winterTies, region, framePattern)` | Scans pair directories; returns a dict keyed by nDays with lists of pair metadata dicts. |
| `tieCodes(phaseTies, track)` | Returns prefix/suffix codes controlling the `use` and `tiefile` line types in the plan. |
| `setupExtraTies(fout, tieFilesData, track, tracknum, year, extraTies)` | Writes bedrock tiefile blocks; returns `tieIntervalCodes` dict if multiple tie-file intervals are in use. |
| `processNDays(fout, pairs, track, tracknum, nDays, year, flags, phaseTies, tieIntervalCodes)` | Writes `subdir`, `use`, and segment lines for one nDays group; returns list of orbit group names. |
| `processAll(fout, allOrbGroups, year, tracknum)` | Writes the final `all YEAR` tiefile that copies all orbit groups. |
| `checkOffsets(imageDir)` | Verifies that `azimuth.offsets` and `range.offsets` exist and match the expected size from the `.dat` file. Returns an error-comment prefix if invalid. |
| `getPrefix(imageDir, phaseTies)` | Returns `'## Exclude '`, `'## Special '`, an offset-error string, or `''` depending on directory state. |
