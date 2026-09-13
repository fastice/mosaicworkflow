# makevelnoclean

Reproduces velocity directories without culling, creating `velocity_nocull` (or `velocity_SVConst`) siblings of existing `velocity` directories. Runs `mosaic3d` via threaded `runOff` scripts. Can also redo existing culled velocity directories with updated flags.

Run from the directory above the track directories (e.g., the main project directory).

---

## Usage

```
makevelnoclean.py [options]
```

No positional arguments. Scans every `track-*` directory below the current one (or the current directory itself when run from inside a track), unless `-toRun` names an explicit list.

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `-reset` | False | Reprocess even if `velocity_nocull/mosaicOffsets.vx` already exists |
| `-redoculled` | False | Rerun existing culled `velocity/` directories (adds `-vzFlag 3` to `runOff`). **Mutually exclusive with the nocull path** — when set, `velocity_nocull/` is not touched and `-reset` has no effect. Rebuilding both sets takes two invocations: one with `-redoculled`, one with `-reset` |
| `-tiff` | False | Force GeoTIFF output (`mosaic3d -GTiff`) regardless of `project.yaml` `velThumbOutput`. Also *removes* `-GTiff` from an inherited `runOff` when tiff is not selected, so the format always matches what the project asks for |
| `-threads=N` | 12 | Number of parallel threads (max 30) |
| `-noprompt` | False | Run without interactive prompt |
| `-useHeader` | False | Read output resolution from `tiepoints/vel_thumb_header_*` files |
| `-SVAlongTrack` | — | Use azimuth-varying SV baseline/offset corrections |
| `-SVConst` | — | Use constant range/azimuth SV corrections |
| `-firstdate=YYYY:MM:DD` | 1990:12:31 | Only process velocities whose `mosaicOffsets.vx` date is on or after this date |
| `-lastdate=YYYY:MM:DD` | 2100:12:31 | Only process velocities on or before this date |
| `-reSize` | False | Reset `inputFile` bounding box to natural size (for `makevelstatsregions.py`) |
| `-resolution="X0 Y0 XS YS DX DY"` | None | Override resolution string in `inputFile` |
| `-toRun="['track-X',...]"` | every `track-*` here | Python list of track directories to process. The default is discovered from the working directory and numerically sorted (`track-2` before `track-10`); running from inside a track dir uses `['.']`. Errors if no `track-*` is found and none was given |
| `-help` | — | Print usage and exit |

---

## What it does

For each `velocity` directory matching the date filter:

1. Creates a `velocity_nocull` (or `velocity_SVConst`) sibling directory
2. Builds a new `inputFile` substituting uncleaned offset file paths (`*.interp.da`, `*.interp.dr`) and optionally appending fast-offset lines
3. Copies `runOff` from the original velocity directory and injects `-SVConst` (and `-vzFlag 3` when `-redoculled`, and `-GTiff` when tiff output is selected) into the `mosaic3d` command line
4. Runs each `runOff` in a thread via `csh`

---

## Part of the mosaicworkflow package.
