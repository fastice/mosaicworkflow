# velocityStats.py

Compute per-frame-range mean velocity and outlier-rejection thresholds from a stack of per-pair velocity mosaics. Run from a `velocityStats/X-Y/` directory (or from the parent `velocityStats/` directory — it loops over all `X-Y` subdirectories it finds). The outputs are used by `mosaic3d` as error weights for culling bad pairs during the final velocity mosaic.

---

## Usage

```
velocityStats.py [-nocull] [-sumVz] [-velmap=velmap] [-region=region] [-help]
```

---

## Options

| Option | Description |
|--------|-------------|
| `-nocull` | Use `velocity_nocull/` subdirectories instead of `velocity/` |
| `-sumVz` | Also accumulate and write the vertical (`.vz`) channel |
| `-velmap=FILE` | Base-map velocity file used to fill pixels with few observations; default is the region default from `sarfunc.defaultRegionDefs` |
| `-region=NAME` | Override region detection; valid values are `greenland` and `antarctica`. When omitted the region is read from `../../region` or inferred from the geodat EPSG |
| `-help` | Print usage and exit |

---

## Behaviour

1. **Discovers frame-range directories** — globs for `X-Y` subdirectories in the current directory (e.g. `0000-0009`, `0010-0019`).
2. **Collects velocity files** — lists `../*_*/velocity` (or `velocity_nocull`) directories relative to the `velocityStats/` parent.
3. **Two-pass outlier rejection** per frame range:
   - **Pass 1** — flags pixels where the per-pair speed exceeds 3× the base-map speed (for pixels faster than 100 m/yr). Accumulates mean and variance.
   - **Pass 2** — re-reads the pass-1 mean/sigma and rejects pixels more than 3σ from the mean.
4. **Writes outputs** to `X-Y/velocity` (or `velocity_nocull`):
   - `.vx`, `.vy` — mean velocity components (m/yr, big-endian float32)
   - `.ex`, `.ey` — per-pixel error thresholds (m/yr)
   - `.navg` — count of accepted observations per pixel
   - `.raw.vx`, `.raw.vy` — pass-1 mean before base-map blending
   - Matching `.geodat` sidecars for all binary files
   - `.vz` (if `-sumVz`) — mean vertical velocity channel

---

## Sigma (error threshold) computation

The `.ex`/`.ey` thresholds fed to `mosaic3d` are computed adaptively per pixel based on the number of observations and the base-map speed:

| Speed regime | Threshold logic |
|---|---|
| ≥ 7 observations | Sample standard deviation of the stack |
| 1–6 observations | 15 + 0.2 × \|base velocity\| m/yr, tightened for mid-speed pixels |
| No observations | Filled from base map; error set to large value |
| Slow (< 80 m/yr) | Capped at 30 m/yr |
| Mid (80–300 m/yr) | Lower-bounded at 7.5% of base speed + 7.5 m/yr |
| Fast (> 300 m/yr) | Lower-bounded at 75 m/yr |
| Extreme (> 1500 m/yr) | Lower-bounded at 75 + 5% of base speed |

For pixels inside a high-variability shapefile region (`sigmaShape` from `sarfunc.defaultRegionDefs`), the threshold is taken directly from the data variance regardless of count.

---

## Directory layout expected

```
velocityStats/
    0000-0009/        ← frame-range output dirs
    0010-0019/
    ...
    ../               ← track dir
        <orbit>_<frame>/
            velocity/
                mosaicOffsets.vx
                mosaicOffsets.vy
                mosaicOffsets.vx.geodat
                ...
```

---

## Key internal functions

| Function | Description |
|----------|-------------|
| `setupFirstVel` | Allocates accumulators and interpolates the base-map onto the output grid |
| `findKeepers` | Flags outlier pixels using either the base-map speed threshold (pass 1) or the running mean/sigma (pass 2) |
| `sumStats` | Accumulates sum and sum-of-squares for mean/variance |
| `normStats` | Normalises accumulators and builds adaptive error thresholds |
| `makeSigMask` | Rasterises the high-variability shapefile into a per-pixel scale factor |
| `getRegion` | Detects region from `../../region` file or geodat EPSG |

---

## Programs called

None — `velocityStats` is self-contained; it reads binary flat files and geodats via `utilities` and writes results directly.
