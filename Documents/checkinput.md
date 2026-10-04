# checkinput

Checks every file a `mosaic3d` master input file refers to, and checks that each
baseline and offset-parameter fit has a usable solution. It is the usual check
after [`makeMaster`](makeMaster.md) builds a master, and `makeMaster` prints the
command as its suggested next step. It works on any master: Sentinel-1 or TSX
offsets, the S1 unwrapped-phase project, or a mix of these.

```
checkinput.py <inputFile> [--phaseTieThresh N] [--phaseCoverage F] [--sigmaTieThresh S]
checkinput    <inputFile> ...          # same program
```

Without the two exclude options, `checkinput` only reads files: it reports
problems and changes nothing. A full Sentinel-1 offsets master (about 20,000
lines) takes about 40 minutes over NFS. The 3,956-line phase master takes
about 40 seconds.

This program moved here from `/Users/ian/insarScripts/bin/checkinput.py` on
2026-10-04. The old copy was renamed `checkinput.py.movedtomosaicworkflow`.

## Options

| Option | Default | Effect |
|---|---|---|
| `--sigmaTieThresh S` | off | Warn when a phase baseline's sigma is greater than S. |
| `--phaseTieThresh N` | off | **Writes `Exclude.pending`** in the frame dir of every phase product with N or fewer phase tiepoints. |
| `--phaseCoverage F` | off | For a phase product with 10 or fewer tiepoints, computes the fraction of valid phase pixels. If it is below F, **writes `Exclude.pending`**. The old `createMaster.ipynb` used 0.1. |

These values are fixed in the code: `azSigma = 1.0` (az.est sigma warning),
`rgSigma = 1.0` (rBaseline sigma warning) and `svConst = True` (see below).

## What it checks, per line

Each data line has 9 fields:
`phase geodat baseline nDays weight azOff azParams rgOff rParams`.
A field containing `None` is skipped.

1. **Offsets** (`azimuth.offsets`, `range.offsets`) must exist and hold at
   least `nr × na × 4` bytes, using the size from the `.vrt` or `.dat`. If
   only the TIFF+VRT version exists, `<name>.vrt` is used and the size check
   is skipped.
2. **geodat** is read to get the acquisition date and size.
3. **Azimuth parameters** (`az.est`) and **range baseline** (`rBaseline`): both
   need a solution and a sigma. A warning is printed for sigma = 0, for sigma above
   `azSigma`/`rgSigma`, or for 10 or fewer tiepoints.
4. **Phase baseline** (skipped for `nophase`): the same tests, with sigma checked
   against `--sigmaTieThresh`. Then the exclude options above are applied.

At the end it prints, for each project directory, the date of the most recent
TSX and Sentinel-1 acquisition. The project directory is the fourth path
component from the end of the geodat path, for example `Sentinel1` or
`SentinelGreenlandUnwrapNew`.

### Parameter file formats

| Form | How it is read |
|---|---|
| Legacy text (`az.est`, `rBaseline`, `baseline.10x2`) | `sigma = …` and `Ntiepoints = used/given` lines, plus a numeric parameter line. With `svConst`, `.const` is appended to `az.est` and `.deltabp` to `rBaseline`. |
| YAML (`az.est.yaml`, `rBaseline.deltabp.yaml`, `baseline.30x6.yaml`) | The first `sigma:` key and the first `nTiepointsUsed:` or `nTiepoints:` key. Sigma ≤ 0 is the no-solution marker. Parsed line by line because `rBaseline.deltabp.yaml` keeps a legacy `;` header. `.const`/`.deltabp` are never appended to a `.yaml` name. |
| Legacy name, YAML on disk | If `az.est.const` (or `rBaseline.deltabp`) is missing but `<name>.yaml` exists, the YAML file is checked instead. |

For a phase-baseline YAML, the sigma read is the one from the solution block
(`noSquint:`), which is the fit `tiepoints` kept, with or without the
ionosphere screen (`usingIon`).

### Phase files

The coverage test reads legacy flat files as big-endian float32 (`nr × na`
from the geodat). It reads `.vrt` phase files (the S1 phase-only products'
`phase.uw.interp.vrt`) through GDAL. `Exclude.pending` goes in the frame dir:
the part before `/motion` for legacy `…/motion/phase.reflat.uw`, or the phase
file's own directory for `…/<frame>/phase.uw.interp.vrt`.

## Output lines

| Message | Meaning |
|---|---|
| `<n> Missing <file>` | A referenced file does not exist, even after the `.vrt`/`.yaml` fallbacks. |
| `<n> no solution for <file>` | A parameter fit has no solution: no sigma, sigma ≤ 0, or 0 tiepoints. |
| `<n> only <k> tiepoints for <file>` | The fit used 10 or fewer tiepoints. |
| `<n> warning sigma = <s> for <file>` | Sigma is 0 or above the relevant threshold. |
| `Excluding …` | An `Exclude.pending` file was written (only with the exclude options). |

`<n>` is the 1-based line number among the data lines.

## Typical results (2026-10-03/04)

- **S1 phase master** (`/Volumes/insar3/ian/SentinelGreenlandUnwrapNew/masterInput/track-all/inputFile`):
  no missing files, no "no solution" entries, and 98 frames with 10 or fewer
  phase tiepoints. The most recent acquisition is 2026-04-18.
- **S1 offsets master** (`/Volumes/insar9/ian/Sentinel1/Release/masterInput/track-all/inputFile`):
  runs to completion. Before the 2026-10-04 changes it crashed on the first
  VRT-migrated frame. Real findings are 54 YAML parameter files with no
  solution and 934 sigma warnings. 1,082 false "Missing" lines from legacy
  `az.est`/`rBaseline` names (now YAML only on disk) were removed by the
  YAML fallback.
