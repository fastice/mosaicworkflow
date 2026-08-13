# grepdate

Reports the processing status of all SAR image pairs for a given year. For each product directory, shows processing state, temporal baseline to the next image, orbit/frame, date, SLC presence, and day-of-year. Autodetects sensor from the current path (Sentinel, TSX, CSK, NISAR).

Supersedes the legacy csh scripts `greptops` and `greptopsuw` (which are kept as thin wrappers for backwards compatibility).

Run from within a track directory.

---

## Usage

```
grepdate [options] year
```

| Argument | Description |
|----------|-------------|
| `year` | Four-digit year to report on |

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--sensor NAME` | auto | Sensor name (`S1`, `TSX`, `CSK`) — only needed if not detectable from path |
| `--firstdate YYYY-MM-DD` | Jan 1 of year | Override start of date range |
| `--lastdate YYYY-MM-DD` | Dec 31 of year | Override end of date range |
| `--frame FRAME` | `*` (all) | Single frame number, or `frame1,frame2` for a range |

---

## Output columns

```
<status> <dt>d <orbit>_<frame>  <date>  <slc>  <doy>
```

| Field | Meaning |
|-------|---------|
| `status` | `xf-`/`xf v`/`x.-`/`x-v` = processed (fast/slow × velocity); `.--` = runboth only; `o--` = no runboth |
| `dt` | Days to next image of same frame |
| `slc` | `+` if `.slc` file present, `-` if not |
| `doy` | Day of year |

---

## Examples

All frames for 2024:
```
grepdate 2024
```

Single frame, custom date range:
```
grepdate 2024 --frame 450 --firstdate 2024-06-01 --lastdate 2024-09-30
```

---

## Part of the mosaicworkflow package.
