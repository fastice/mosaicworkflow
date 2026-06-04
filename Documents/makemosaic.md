# makemosaic.py

Front-end driver that loops over a date range and calls `setupquarters.py` once per time interval to produce velocity mosaics.

---

## Usage

```
makemosaic.py [options]
```

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--firstdate YYYY-MM-DD` | `2015-01-01` | Start of date range. Adjusted to align with interval boundaries. |
| `--lastdate YYYY-MM-DD` | today | End of date range. |
| `--inputFile FILE` | None | Input file for radar data (overrides template default). |
| `--interval` | `monthly` | Time interval: `s1cycle`, `s1-12day`, `monthly`, `quarterly`, `annual`, `multiYear`, `quarterlyJFM`. |
| `--noReprocess` | False | Skip reprocessing; only reformat existing data. |
| `--noLandsat` | False | Do not include Landsat data. |
| `--landsatPath PATH` | None | Path to directory containing Landsat list files. |
| `--LSFitType` | `custom.mask` | Landsat fit type: `custom.mask`, `custom`, `static`, `static.mask`. |
| `--mosaicMask FILE` | None | Mask applied by mosaic3d during mosaicking. |
| `--outputMask FILE` | None | Shapefile for masking final output. |
| `--baseFlags STR` | None | Override default flags passed to the mosaicker. |
| `--keepFast` | False | Suppress fast-area culling (retains mélange). |
| `--check` | False | Set up command and check masks but do not run. |
| `--template FILE` | `mosaic.template.yaml` | YAML template that defines the mosaic geometry and parameters. |
| `--noLabel` | False | Suppress processing source label in output. |
| `--mosaicsSetupFile FILE` | `mosaicsSetup.yaml` | YAML file with seasonal mask paths and interval configuration. |

---

## Interval types

| Interval | Period | Notes |
|----------|--------|-------|
| `s1cycle` | 6 or 12 days | Follows Sentinel-1 repeat cycle; 12-day before 2016-09-20 and after 2021-12-24, 6-day between. TSX excluded. |
| `s1-12day` | 12 days | Fixed 12-day cadence. TSX excluded. |
| `monthly` | 1 month | Snaps `firstdate` to first of month. |
| `quarterly` | ~3 months | Start month from `seasonData.quarterlyDates` and `quarterlyMonth` in `mosaicsSetup.yaml`. |
| `annual` | ~1 year | Start month from `seasonData.annualFirstMonth`. |
| `multiYear` | entire range | Single product spanning `firstdate`→`lastdate`; loop runs once. |

---

## Configuration files

### `mosaicsSetup.yaml`

Read by `mosaicfunc.readYaml`. Expected keys:

| Key | Description |
|-----|-------------|
| `landsatPath` | Path to Landsat list files (overridden by `--landsatPath`). |
| `mosaicMask` | Explicit mask file (overridden by `--mosaicMask`). |
| `intervalMasks` | Dict mapping interval name → mask path template (with `YY` and `SS` placeholders). |
| `whichMask` | Dict mapping interval name → list of 12 season codes, one per calendar month. |
| `seasonData` | Dict with `quarterlyDates`, `quarterlyMonth`, `annualFirstMonth` for date alignment. |

### `mosaic.template.yaml`

Passed directly to `setupquarters.py` via `--template`. Defines mosaic grid, platform lists, and other mosaicker parameters.

---

## Landsat list files

Located under `landsatPath/listfiles/` with names of the form:

```
Listfile.<YEAR>.<fitType>
```

Files spanning the requested date range are merged into a temporary `mergedList.<year1>-<year2>` file and passed to `setupquarters.py` via `--lsFile`.

---

## Seasonal mask selection

When `intervalMasks` is present in `mosaicsSetup.yaml` and no explicit `--mosaicMask` is given, the mask is selected as follows:

1. Compute the central date of the interval.
2. Look up the season code for that month from `whichMask[interval]`.
3. Substitute `YY` (2-digit year) and `SS` (season code) into the mask path template.

---

## Freestanding programs called

### `setupquarters.py`

The only external program called (via `subprocess.call`). Invoked once per time interval with arguments assembled by `makeCommand()`:

```
setupquarters.py
    [--template FILE]
    --firstdate YYYY-MM-DD
    --lastdate  YYYY-MM-DD
    [--noReprocess]
    [--noCull]          # when --keepFast
    [--noTSX]           # when interval is s1cycle or s1-12day
    [--noLabel]
    [--baseFlags "..."]
    [--outputMask FILE]
    [--mosaicMask FILE]
    [--inputFile FILE]
    [--lsFile mergedList]
```

---

## Processing flow

```
makemosaicArgs()
  └── reads mosaicsSetup.yaml
  └── adjustFirstDate() — snaps firstdate to interval boundary

getLists()            — collect Landsat list files
createMergedList()    — merge into single list file

loop over intervals:
  getMosaicMask()     — select seasonal or explicit mask
  makeCommand()       — assemble setupquarters.py command string
  call(command)       — run setupquarters.py
  incrementDate()     — advance to next interval start
```

---

## Key internal functions

| Function | Description |
|----------|-------------|
| `adjustFirstDate(firstdate, interval, seasonData)` | Snaps `firstdate` to the correct boundary for the chosen interval. |
| `incrementDate(myDate, interval)` | Advances a date by one interval period. |
| `findLastDate(firstDate, myArgs)` | Returns `(endDate, year)` for the interval starting at `firstDate`. |
| `getLists(landsatPath, fitType)` | Globs Landsat list files matching `fitType`. |
| `createMergedList(listFiles, firstDate, lastDate)` | Concatenates list files spanning the date range into one merged file. |
| `getNewMask(interval, firstDate, lastDate, mosaicsSetup)` | Resolves seasonal mask path from template. |
| `getMosaicMask(mosaicsSetup, interval, firstDate, lastDate, keepFast)` | Returns the mask file to use, or None. |
| `makeCommand(firstDate, lastDate, mergedList, mosaicMaskFile, myArgs)` | Builds the `setupquarters.py` command string. |
