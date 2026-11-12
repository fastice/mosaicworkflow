'''
Software versions written into the released products, so every product can be traced to the
commits of the code that made it (beneath the product version in the file name).

Each repo's version is bumped on every commit by the pre-commit hook
(mosaicSource/scripts/preCommitVersionBump.sh): pyproject.toml for the Python packages, VERSION
for mosaicSource. The commit is the short git hash, marked -dirty.<first 8 of the hash of
git diff HEAD> when run from uncommitted changes -- the same convention mosaicSource compiles
into its binaries (scripts/makeVersionHeader.sh).
'''
import glob
import os
import re
import subprocess
from osgeo import gdal


def gitCommit(repoDir):
    ''' Short commit of repoDir, with -dirty.<diff hash> for uncommitted changes '''
    def git(*args, **kwargs):
        return subprocess.run(['git', '-C', repoDir, *args], capture_output=True, **kwargs)
    result = git('rev-parse', '--short', 'HEAD', text=True)
    if result.returncode != 0:
        return 'unknown'
    commit = result.stdout.strip()
    if git('status', '--porcelain', '--untracked-files=no', text=True).stdout.strip():
        diff = git('diff', 'HEAD').stdout
        diffHash = subprocess.run(['git', 'hash-object', '--stdin'], input=diff,
                                  capture_output=True).stdout.decode().strip()[:8]
        commit = f'{commit}-dirty.{diffHash}'
    return commit


def packageVersion(module):
    '''
    (version, commit) of an installed package. The version is read from pyproject.toml beside
    the package (editable installs, so always current) and falls back to the installed metadata.
    '''
    packageDir = os.path.dirname(os.path.abspath(module.__file__))
    repoDir = os.path.dirname(packageDir)
    version = 'unknown'
    pyproject = os.path.join(repoDir, 'pyproject.toml')
    if os.path.exists(pyproject):
        with open(pyproject) as fp:
            match = re.search(r'^version\s*=\s*"([^"]+)"', fp.read(), re.MULTILINE)
        if match:
            version = match.group(1)
    else:
        try:
            from importlib.metadata import version as installedVersion
            version = installedVersion(module.__name__)
        except Exception:
            pass
    return version, gitCommit(repoDir)


def pythonVersions():
    ''' Metadata keys for the Python packages that make the mosaics '''
    import mosaicworkflow
    import mosaicfunc
    meta = {}
    for module in [mosaicworkflow, mosaicfunc]:
        version, commit = packageVersion(module)
        meta[f'{module.__name__}Version'] = version
        meta[f'{module.__name__}Commit'] = commit
    return meta


def cVersionsFromPieces(piecesDir):
    '''
    mosaicSource* metadata keys read back from the mosaic3d pieces (GTiff, or the VRT sidecars
    of flat-binary pieces), so the release records the binary that actually ran. Pieces made by
    different builds are joined with " | ". Pieces from before versioning give no keys.
    '''
    values = {}
    for pieceFile in sorted(glob.glob(f'{piecesDir}/*.tif') + glob.glob(f'{piecesDir}/*.vrt')):
        ds = gdal.Open(pieceFile)
        if ds is None:
            continue
        for key, value in (ds.GetMetadata() or {}).items():
            if key.startswith('mosaicSource'):
                values.setdefault(key, [])
                if value not in values[key]:
                    values[key].append(value)
        ds = None
    return {key: ' | '.join(v) for key, v in values.items()}
