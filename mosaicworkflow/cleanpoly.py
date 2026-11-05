#!/usr/bin/env python3
"""
cleanpoly.py

Python replacement for clean_list + cleantool_velocity.pro + make_velpoly_bad.pro.

For each velocity directory (track-N/<frame>/velocity):
  1. Display the per-frame speed (velocity/mosaicOffsets) at one screen pixel per
     image pixel, with the existing polygons from mosaicOffsetsbadVerts.txt, and
     let the user add polygons/circles or toggle existing ones on/off.
  2. Convert the active polygons (lat/lon/z) to slant range/azimuth with lltora,
     fill them on the offset grid and write the flat indices (uint32, a*nr + r)
     to <frame>/badoffsets_poly.list. cleanoffmerge uses this in preference to
     the IDL badoffsets_poly.idl.
  3. Run the frame's cleanoff script.

mosaicOffsetsbadVerts.txt keeps the IDL format (rows of 4f15.6): a header row
(999, polyId, nVerts, 999) followed by nVerts rows (lat, lon, on/off, z), so
polygons drawn with either tool can be edited with the other.

Keys: p add polygon, c add circle, e toggle polygon on/off, v edit vertices,
Esc cancel, q quit.
Mouse while drawing: left mark, shift-left (or middle, or Backspace) erase
last, right accept.
Vertex mode: left drag moves a vertex, shift-left (or middle) deletes one,
right inserts one on the nearest edge; v or Esc leaves the mode.
"""
import argparse
import datetime
import glob
import os
import sys
import tempfile
from subprocess import call, run
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import scipy.io as io
from osgeo import gdal, osr
from matplotlib.path import Path

vertsName = 'mosaicOffsetsbadVerts.txt'
listName = 'badoffsets_poly.list'
demNorth = '/Volumes/insar7/ian/gimp/gimp2/270m/dem.gimp2.270mex.tif'
demSouth = '/Volumes/insarb/ian/fromJonathan/newAntarctic/dem1.tif'
# display backdrop used by rsfig_lite (the south RAMP mosaic no longer exists)
sarNorth = '/Volumes/insar7/gmap/Mosaics/GreenlandComplete/mosaic100extend/' \
    'mosaic100m'
noData = -2.e9


def cleanPolyArgs():
    parser = argparse.ArgumentParser(
        description='\033[1mDisplay per-frame velocities, edit bad-data '
        'polygons, write badoffsets_poly.list and run cleanoff (Python '
        'version of clean_list)\033[0m',
        epilog='Directories may be velocity dirs or the frame dirs above them. '
        'A quoted wildcard is expanded.')
    parser.add_argument('dirs', nargs='*', help='velocity (or frame) dirs')
    parser.add_argument('--listFile', default=None,
                        help='file with one velocity (or frame) dir per line')
    parser.add_argument('--firstDate', default='2025-01-01',
                        help='only frames whose first image date is on or '
                        'after this (YYYY-MM-DD); 1900-01-01 for all '
                        '[2025-01-01]')
    parser.add_argument('--nCleanoff', type=int, default=4,
                        help='cleanoffs run at once in the background [4]')
    parser.add_argument('--south', action='store_true',
                        help='Antarctic (EPSG:3031) [north, EPSG:3413]')
    parser.add_argument('--noReplace', action='store_true',
                        help=f'Skip dirs that already have {vertsName}')
    parser.add_argument('--indicesOnly', action='store_true',
                        help=f'No display: rebuild {listName} from existing '
                        'polygons (e.g. after offsets are re-measured)')
    parser.add_argument('--noCleanoff', action='store_true',
                        help='Do not run cleanoff afterwards (by default it '
                        f'runs, in the background, when {listName} changed)')
    parser.add_argument('--velocityFile', default='mosaicOffsets',
                        help='velocity file root [mosaicOffsets]')
    parser.add_argument('--dem', default=None,
                        help='DEM for vertex heights [gimp2 270m or dem1]')
    parser.add_argument('--maxv', type=float, default=3000.,
                        help='display clip speed m/yr [3000]')
    parser.add_argument('--exponent', type=float, default=0.1,
                        help='SAR backdrop stretch image**exponent [0.1]')
    args = parser.parse_args()
    if args.listFile is not None:
        args.dirs += [x.strip() for x in open(args.listFile)
                      if len(x.strip()) > 0 and not x.startswith('#')]
    if len(args.dirs) == 0:
        parser.error('no directories given')
    args.firstDate = datetime.datetime.strptime(args.firstDate,
                                                '%Y-%m-%d').date()
    if args.dem is None:
        args.dem = demSouth if args.south else demNorth
    return args


def getDirs(dirArgs):
    ''' Expand globs and return (velocityDir, frameDir) pairs '''
    dirs = []
    for myArg in dirArgs:
        for myDir in sorted(glob.glob(myArg)) if any(
                c in myArg for c in '*?[') else [myArg]:
            myDir = os.path.abspath(myDir.rstrip('/'))
            if os.path.isdir(f'{myDir}/velocity'):
                myDir = f'{myDir}/velocity'
            if not os.path.isdir(myDir):
                print(f'skipping {myDir}: not a directory')
                continue
            dirs.append((myDir, os.path.dirname(myDir)))
    return dirs


def firstImageDate(velDir, velFile):
    ''' First image date from velocity/<velFile>.meta (None if unknown) '''
    metaFile = f'{velDir}/{velFile}.meta'
    if not os.path.exists(metaFile):
        return None
    for line in open(metaFile):
        if line.startswith('First Image Date'):
            try:
                return datetime.datetime.strptime(
                    line.split('=')[1].strip(), '%b:%d:%Y').date()
            except ValueError:
                return None
    return None

# ---------------------------------------------------------------------------
# Polygon file (IDL mosaicOffsetsbadVerts.txt format)
# ---------------------------------------------------------------------------


def readVerts(vertsFile):
    ''' Return list of polygons {id, lat, lon, on, z} '''
    polys = []
    if not os.path.exists(vertsFile):
        return polys
    data = np.loadtxt(vertsFile, ndmin=2)
    i = 0
    while i < len(data):
        row = data[i]
        if row[0] == 999 and row[3] == 999:
            nVerts = int(row[2])
            verts = data[i + 1:i + 1 + nVerts]
            if len(verts) == nVerts and nVerts > 0:
                polys.append({'id': int(row[1]), 'lat': verts[:, 0].copy(),
                              'lon': verts[:, 1].copy(),
                              'on': verts[0, 2] == 1, 'z': verts[:, 3].copy()})
            i += 1 + nVerts
        else:
            i += 1
    return polys


def writeVerts(vertsFile, polys):
    ''' Write polygons in the IDL format, keeping the old file as ...old '''
    if os.path.exists(vertsFile):
        call(['cp', vertsFile, vertsFile + 'old'])
    rows = []
    for poly in polys:
        rows.append([999, poly['id'], len(poly['lat']), 999])
        for lat, lon, z in zip(poly['lat'], poly['lon'], poly['z']):
            rows.append([lat, lon, 1 if poly['on'] else 0, z])
    if len(rows) == 0:
        rows = [[0, 0, 0, 0]]
    np.savetxt(vertsFile, np.array(rows), fmt='%15.6f', delimiter='')

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


class Projection:
    ''' Polar stereographic (m) <-> lat/lon '''

    def __init__(self, south):
        ps = osr.SpatialReference()
        ps.ImportFromEPSG(3031 if south else 3413)
        ll = osr.SpatialReference()
        ll.ImportFromEPSG(4326)
        for srs in [ps, ll]:
            srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        self.toLL = osr.CoordinateTransformation(ps, ll)
        self.toXY = osr.CoordinateTransformation(ll, ps)

    def xyToLL(self, x, y):
        pts = np.array(self.toLL.TransformPoints(
            np.column_stack([x, y]).astype(float)))
        return pts[:, 1], pts[:, 0]

    def llToXY(self, lat, lon):
        pts = np.array(self.toXY.TransformPoints(
            np.column_stack([lon, lat]).astype(float)))
        return pts[:, 0], pts[:, 1]


def demHeights(demFile, lat, lon):
    ''' Bilinear DEM heights at lat/lon '''
    ds = gdal.Open(demFile)
    srs = osr.SpatialReference(wkt=ds.GetProjection())
    ll = osr.SpatialReference()
    ll.ImportFromEPSG(4326)
    for s in [srs, ll]:
        s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    pts = np.array(osr.CoordinateTransformation(ll, srs).TransformPoints(
        np.column_stack([lon, lat]).astype(float)))
    gt = ds.GetGeoTransform()
    band = ds.GetRasterBand(1)
    z = np.zeros(len(lat))
    for k, (x, y) in enumerate(pts[:, 0:2]):
        # pixel-centre coordinates
        c = (x - gt[0]) / gt[1] - 0.5
        r = (y - gt[3]) / gt[5] - 0.5
        c0 = int(np.clip(np.floor(c), 0, ds.RasterXSize - 2))
        r0 = int(np.clip(np.floor(r), 0, ds.RasterYSize - 2))
        w = band.ReadAsArray(c0, r0, 2, 2).astype(float)
        fc, fr = np.clip(c - c0, 0, 1), np.clip(r - r0, 0, 1)
        z[k] = (w[0, 0] * (1 - fc) * (1 - fr) + w[0, 1] * fc * (1 - fr) +
                w[1, 0] * (1 - fc) * fr + w[1, 1] * fc * fr)
    return z


def findGeodat(frameDir):
    ''' Same choice as make_velpoly_bad: first of geodat?x* geodat10x* '''
    files = sorted(glob.glob(f'{frameDir}/geodat?x*')) + \
        sorted(glob.glob(f'{frameDir}/geodat10x*'))
    if len(files) == 0:
        return None
    return files[0]


def runLLtoRA(geodat, lat, lon, z, workDir):
    ''' Run lltora, return slant range, azimuth (single-look pixels) '''
    with tempfile.NamedTemporaryFile('w', dir=workDir, prefix='tmpCleanPoly',
                                     suffix='.lle', delete=False) as fp:
        for vals in zip(lat, lon, z):
            print('%15.6f%15.6f%15.6f' % vals, file=fp)
        print('&', file=fp)
        llFile = fp.name
    try:
        result = run(['lltora', geodat, llFile], capture_output=True,
                     text=True)
    finally:
        os.remove(llFile)
    # data rows follow the '# n' line and end at '&'
    rows, inData = [], False
    for line in result.stdout.splitlines():
        if line.startswith('#'):
            inData = True
            continue
        if line.strip() == '&':
            break
        if inData:
            pieces = line.split()
            if len(pieces) == 5:
                rows.append([float(x) for x in pieces])
    if len(rows) != len(lat):
        sys.exit(f'lltora returned {len(rows)} points for {len(lat)}:\n'
                 f'{result.stdout}\n{result.stderr}')
    rows = np.array(rows)
    return rows[:, 0], rows[:, 1]


def readOffsetGrid(frameDir):
    ''' R0, A0, nr, na, rLooks, aLooks from azimuth.offsets.dat/offsets.dat '''
    for datFile in ['azimuth.offsets.dat', 'offsets.dat']:
        if os.path.exists(f'{frameDir}/{datFile}'):
            vals = open(f'{frameDir}/{datFile}').read().split()[0:6]
            R0, A0, nr, na, rLooks, aLooks = [float(x) for x in vals]
            return R0, A0, int(nr), int(na), rLooks, aLooks
    return None


def polygonIndices(x, y, nr, na):
    ''' Flat indices (y*nr + x) of offset pixels inside polygon (x, y) in
    offset-pixel units. As IDL polyfillv: vertices are truncated to integers,
    then pixel (i, j) is in when its centre (i + 0.5, j + 0.5) is. '''
    x, y = np.trunc(x), np.trunc(y)
    x0, x1 = max(int(np.floor(x.min())), 0), min(int(np.ceil(x.max())), nr-1)
    y0, y1 = max(int(np.floor(y.min())), 0), min(int(np.ceil(y.max())), na-1)
    if x1 < x0 or y1 < y0:
        return np.array([], dtype='u4')
    xg, yg = np.meshgrid(np.arange(x0, x1 + 1), np.arange(y0, y1 + 1))
    inside = Path(np.column_stack([x, y])).contains_points(
        np.column_stack([xg.ravel() + 0.5, yg.ravel() + 0.5]))
    return (yg.ravel()[inside] * nr + xg.ravel()[inside]).astype('u4')


def polyIndicesForFrame(frameDir, polys):
    ''' Union of offset indices for the active polygons (None if no grid) '''
    grid = readOffsetGrid(frameDir)
    geodat = findGeodat(frameDir)
    if grid is None or geodat is None:
        print(f'{frameDir}: missing offsets.dat or geodat, no {listName}')
        return None
    R0, A0, nr, na, rLooks, aLooks = grid
    bad = [np.array([], dtype='u4')]
    for poly in polys:
        if not poly['on']:
            continue
        r, a = runLLtoRA(geodat, poly['lat'], poly['lon'], poly['z'],
                         frameDir)
        bad.append(polygonIndices((r - R0) / rLooks, (a - A0) / aLooks,
                                  nr, na))
    return np.unique(np.concatenate(bad)).astype('u4')


def makePolyList(frameDir, velDir):
    ''' Write badoffsets_poly.list from velocity/mosaicOffsetsbadVerts.txt '''
    vertsFile = f'{velDir}/{vertsName}'
    if not os.path.exists(vertsFile):
        print(f'{frameDir}: no {vertsName}, {listName} not written')
        return False
    bad = polyIndicesForFrame(frameDir, readVerts(vertsFile))
    if bad is None:
        return False
    # compare with what cleanoff applies now (.list, else the IDL .idl)
    listFile = f'{frameDir}/{listName}'
    idlFile = f'{frameDir}/badoffsets_poly.idl'
    if os.path.exists(listFile):
        old = np.fromfile(listFile, dtype='u4')
    elif os.path.exists(idlFile):
        old = np.unique(io.readsav(idlFile)['idel_poly']).astype('u4')
    else:
        old = np.array([], dtype='u4')
    if np.array_equal(old, bad):
        print(f'{frameDir}: polygon offsets unchanged ({len(bad)})')
        return False
    bad.tofile(listFile)
    print(f'{frameDir}: {len(bad)} offsets in {listName} (was {len(old)})')
    return True


def runCleanoff(frameDir):
    ''' Run the frame's cleanoff, output to cleanoff.cleanpoly.log '''
    log = f'{frameDir}/cleanoff.cleanpoly.log'
    status = call(f'cd {frameDir} ; setenv OPENBLAS_NUM_THREADS 1 ; '
                  f'cleanoff >& {log}', shell=True, executable='/bin/csh')
    print(f'\ncleanoff {"done" if status == 0 else "FAILED"} for {frameDir}'
          f'{"" if status == 0 else f" (see {log})"}')
    return frameDir, status

# ---------------------------------------------------------------------------
# Velocity image
# ---------------------------------------------------------------------------


def readSpeed(velDir, velFile):
    ''' Return speed (row 0 = north), x, y of the first pixel centre and the
    pixel size (m) '''
    root = f'{velDir}/{velFile}'
    if os.path.exists(root + '.vrt'):
        gdal.PushErrorHandler('CPLQuietErrorHandler')
        ds = gdal.Open(root + '.vrt')
        gdal.PopErrorHandler()
        vx = ds.GetRasterBand(1).ReadAsArray().astype(float)
        vy = ds.GetRasterBand(2).ReadAsArray().astype(float)
        gt = ds.GetGeoTransform()
        xc, yc = gt[0] + 0.5 * gt[1], gt[3] + 0.5 * gt[5]
        dx, dy = gt[1], -gt[5]
    elif os.path.exists(root + '.vx'):
        vals = []
        for line in open(root + '.vx.geodat'):
            if line.startswith((';', '#', '&')):
                continue
            vals += [float(v) for v in line.split()]
        nx, ny, dx, dy, x0, y0 = vals[0:6]
        nx, ny = int(nx), int(ny)
        # GrIMP raw xy images: row 0 at the southern origin (pixel centre, km)
        vx = np.flipud(np.fromfile(root + '.vx', '>f4').reshape(ny, nx))
        vy = np.flipud(np.fromfile(root + '.vy', '>f4').reshape(ny, nx))
        vx, vy = vx.astype(float), vy.astype(float)
        xc, yc = x0 * 1000., y0 * 1000. + (ny - 1) * dy
    else:
        return None
    speed = np.hypot(vx, vy)
    speed[np.logical_or(vx <= noData * 0.1, ~np.isfinite(speed))] = np.nan
    return speed, xc, yc, dx, dy


def sarBackdrop(shape, xc, yc, dx, dy, south, exponent):
    ''' SAR mosaic brightness (0-1) sampled (nearest) at the velocity pixel
    centres, stretched as rsfig_lite: bytscl(image^exponent) - 0.8 min '''
    if south or not os.path.exists(sarNorth):
        return np.ones(shape)
    ny, nx = shape
    vals = []
    for line in open(sarNorth + '.geodat'):
        if not line.startswith((';', '#', '&')):
            vals += [float(v) for v in line.split()]
    mx, my, mdx, mdy, mx0, my0 = vals[0:6]
    mx, my = int(mx), int(my)
    sar = np.memmap(sarNorth, dtype='>f4', mode='r', shape=(my, mx))
    # mosaic row 0 is its southern edge, origin is a pixel centre (km)
    cols = np.rint((xc + np.arange(nx) * dx - mx0 * 1000.) / mdx).astype(int)
    rows = np.rint((yc - np.arange(ny) * dy - my0 * 1000.) / mdy).astype(int)
    okC, okR = (cols >= 0) & (cols < mx), (rows >= 0) & (rows < my)
    image = np.zeros(shape)
    if okC.any() and okR.any():
        image[np.ix_(okR, okC)] = sar[rows[okR]][:, cols[okC]]
    image[image <= 0] = 2
    image = image ** exponent
    image = np.floor(255.999 * (image - image.min()) /
                     max(image.max() - image.min(), 1e-9))
    image = image - 0.8 * image.min()
    return image / max(image.max(), 1e-9)


def speedImage(speed, backdrop, maxv):
    ''' rsfig_lite colours: hue = log speed, saturation rising with speed,
    value = SAR backdrop; no-data pixels show the backdrop in gray '''
    from matplotlib.colors import hsv_to_rgb
    v = np.nan_to_num(speed, nan=0.)
    v = np.minimum(v, maxv)
    v[(v > 0) & (v < 1.5)] = 1.5
    hue = np.maximum(np.log(np.maximum(v, 1e-9)), 0) / np.log(maxv)
    sat = np.minimum((0.5 + v / 125.) / 1.5, 1.)
    hue[v == 0], sat[v == 0] = 0., 0.
    rgb = hsv_to_rgb(np.dstack([np.clip(hue, 0, 1), sat, backdrop]))
    return (255 * rgb).astype(np.uint8)


class PolyEditor:
    ''' One-to-one pixel display of a frame's speed with polygon editing '''

    def __init__(self, velDir, polys, speed, xc, yc, dx, dy, proj, args):
        import tkinter as tk
        from PIL import Image, ImageTk
        self.tk = tk
        self.polys, self.proj, self.args = polys, proj, args
        self.xc, self.yc, self.dx, self.dy = xc, yc, dx, dy
        self.speed = speed
        self.modified, self.save = False, False
        self.mode, self.clicks, self.tmpItems = None, [], []
        self.drag = None
        ny, nx = speed.shape
        #
        self.root = tk.Tk()
        self.root.title(f'cleanpoly: {os.path.dirname(velDir)}')
        maxW = self.root.winfo_screenwidth() - 80
        maxH = self.root.winfo_screenheight() - 160
        frame = tk.Frame(self.root)
        frame.pack(fill='both', expand=True)
        self.canvas = tk.Canvas(frame, width=min(nx, maxW),
                                height=min(ny, maxH), bg='black',
                                highlightthickness=0,
                                scrollregion=(0, 0, nx, ny))
        if nx > maxW:
            xs = tk.Scrollbar(frame, orient='horizontal',
                              command=self.canvas.xview)
            xs.pack(side='bottom', fill='x')
            self.canvas.configure(xscrollcommand=xs.set)
        if ny > maxH:
            ys = tk.Scrollbar(frame, orient='vertical',
                              command=self.canvas.yview)
            ys.pack(side='right', fill='y')
            self.canvas.configure(yscrollcommand=ys.set)
        self.canvas.pack(side='left', fill='both', expand=True)
        backdrop = sarBackdrop(speed.shape, xc, yc, dx, dy, args.south,
                               args.exponent)
        img = speedImage(speed, backdrop, args.maxv)
        self.photo = ImageTk.PhotoImage(Image.fromarray(img, mode='RGB'))
        self.canvas.create_image(0, 0, image=self.photo, anchor='nw')
        #
        self.help = tk.StringVar()
        self.status = tk.StringVar()
        tk.Label(self.root, textvariable=self.help, anchor='w',
                 font=('Helvetica', 12, 'bold')).pack(fill='x')
        tk.Label(self.root, textvariable=self.status, anchor='w',
                 font=('Courier', 11)).pack(fill='x')
        self.setHelp()
        self.drawPolys()
        self.root.bind('<Key>', self.onKey)
        self.canvas.bind('<Motion>', self.onMotion)
        self.canvas.bind('<Button-1>', lambda e: self.onClick(e, 1))
        self.canvas.bind('<Button-2>', lambda e: self.onClick(e, 2))
        self.canvas.bind('<Button-3>', lambda e: self.onClick(e, 3))
        # no middle button needed: shift-left acts as middle
        self.canvas.bind('<Shift-Button-1>', lambda e: self.onClick(e, 2))
        self.canvas.bind('<B1-Motion>', self.onDrag)
        self.canvas.bind('<ButtonRelease-1>', self.onRelease)
        self.root.protocol('WM_DELETE_WINDOW', self.quit)
        self.root.focus_force()

    def run(self):
        self.root.mainloop()
        return self.save

    def setHelp(self, text=None):
        self.help.set(text or '"q" quit (accept), "p" add poly, '
                      '"c" add circle, "e" toggle poly on/off, '
                      '"v" edit vertices')

    # coordinate conversions (canvas pixel edges: column c spans c..c+1)
    def canvasToXY(self, cx, cy):
        return (self.xc + (cx - 0.5) * self.dx,
                self.yc - (cy - 0.5) * self.dy)

    def xyToCanvas(self, x, y):
        return ((x - self.xc) / self.dx + 0.5,
                (self.yc - y) / self.dy + 0.5)

    def eventCanvas(self, event):
        return self.canvas.canvasx(event.x), self.canvas.canvasy(event.y)

    def polyCanvas(self, poly):
        x, y = self.proj.llToXY(poly['lat'], poly['lon'])
        return self.xyToCanvas(x, y)

    def drawPolys(self):
        self.canvas.delete('poly')
        for poly in self.polys:
            cx, cy = self.polyCanvas(poly)
            pts = np.column_stack([np.append(cx, cx[0]),
                                   np.append(cy, cy[0])]).ravel().tolist()
            self.canvas.create_line(*pts, width=3, tags='poly',
                                    fill='#3080ff' if poly['on']
                                    else '#909090')
            if self.mode == 'vertex':
                for x, y in zip(cx, cy):
                    self.canvas.create_rectangle(x - 3, y - 3, x + 3, y + 3,
                                                 outline='yellow',
                                                 tags='poly')

    def onMotion(self, event):
        cx, cy = self.eventCanvas(event)
        col, row = int(cx), int(cy)
        ny, nx = self.speed.shape
        if 0 <= col < nx and 0 <= row < ny:
            x, y = self.canvasToXY(col + 0.5, row + 0.5)
            lat, lon = self.proj.xyToLL([x], [y])
            sp = self.speed[row, col]
            spText = 'no data' if np.isnan(sp) else f'{sp:8.1f} m/yr'
            self.status.set(f'x {x/1000.:9.2f} y {y/1000.:9.2f} km  '
                            f'lat {lat[0]:9.5f} lon {lon[0]:10.5f}  '
                            f'speed {spText}')

    def onKey(self, event):
        key = event.char.lower()
        if event.keysym == 'Escape' or (self.mode == 'vertex' and key == 'v'):
            self.clearTmp()
            self.mode, self.clicks, self.drag = None, [], None
            self.drawPolys()
            self.setHelp()
            return
        if event.keysym == 'BackSpace' and self.mode in ['poly', 'circle']:
            if len(self.clicks) > 0:
                self.clicks.pop()
            self.drawTmp()
            return
        if self.mode is not None and key != 'q':
            return
        if key in 'qx' and key != '':
            self.quit()
        elif key == 'p':
            self.mode, self.clicks = 'poly', []
            self.setHelp('POLYGON: left mark, shift-left/Backspace erase '
                         'last, right accept')
        elif key == 'c':
            self.mode, self.clicks = 'circle', []
            self.setHelp('CIRCLE: left centre then radius point, '
                         'shift-left/Backspace erase, right accept')
        elif key == 'e':
            self.mode = 'edit'
            self.setHelp('EDIT: left click near a polygon to toggle it')
        elif key == 'v':
            self.mode = 'vertex'
            self.drawPolys()
            self.setHelp('VERTICES: left drag move, shift-left delete, '
                         'right insert on edge; "v"/Esc done')

    def clearTmp(self):
        for item in self.tmpItems:
            self.canvas.delete(item)
        self.tmpItems = []

    def circlePoints(self):
        (cx0, cy0), (cx1, cy1) = self.clicks[0:2]
        rad = np.hypot(cx1 - cx0, cy1 - cy0)
        theta = np.arange(20) * 2 * np.pi / 20
        return cx0 + rad * np.cos(theta), cy0 + rad * np.sin(theta)

    def drawTmp(self):
        self.clearTmp()
        if self.mode == 'poly' and len(self.clicks) > 1:
            self.tmpItems.append(self.canvas.create_line(
                *np.array(self.clicks).ravel().tolist(), width=2,
                fill='yellow'))
        elif self.mode == 'circle' and len(self.clicks) == 2:
            cx, cy = self.circlePoints()
            pts = np.column_stack([np.append(cx, cx[0]),
                                   np.append(cy, cy[0])]).ravel().tolist()
            self.tmpItems.append(self.canvas.create_line(*pts, width=2,
                                                         fill='yellow'))
        for (cx, cy) in self.clicks:
            self.tmpItems.append(self.canvas.create_oval(
                cx - 2, cy - 2, cx + 2, cy + 2, outline='yellow'))

    def onClick(self, event, button):
        cx, cy = self.eventCanvas(event)
        if self.mode == 'edit' and button == 1:
            self.toggleNearest(cx, cy)
            self.mode = None
            self.setHelp()
            return
        if self.mode == 'vertex':
            self.vertexClick(cx, cy, button)
            return
        if self.mode not in ['poly', 'circle']:
            return
        if button == 1:
            if self.mode == 'circle' and len(self.clicks) == 2:
                self.clicks[1] = (cx, cy)
            else:
                self.clicks.append((cx, cy))
        elif button == 2 and len(self.clicks) > 0:
            self.clicks.pop()
        elif button == 3:
            self.acceptShape()
            return
        self.drawTmp()

    def acceptShape(self):
        if self.mode == 'poly' and len(self.clicks) >= 3:
            cx, cy = np.array(self.clicks).T
        elif self.mode == 'circle' and len(self.clicks) == 2:
            cx, cy = self.circlePoints()
        else:
            cx = None
        if cx is not None:
            x, y = self.canvasToXY(cx, cy)
            lat, lon = self.proj.xyToLL(x, y)
            z = demHeights(self.args.dem, lat, lon)
            newId = max([p['id'] for p in self.polys], default=-1) + 1
            self.polys.append({'id': newId, 'lat': lat, 'lon': lon,
                               'on': True, 'z': z})
            self.modified = True
            self.drawPolys()
        self.clearTmp()
        self.mode, self.clicks = None, []
        self.setHelp()

    def nearestVertex(self, cx, cy, maxDist=8):
        ''' (poly, index) of the vertex within maxDist canvas pixels '''
        best, dMin = None, maxDist
        for poly in self.polys:
            px, py = self.polyCanvas(poly)
            d = np.hypot(px - cx, py - cy)
            if d.min() <= dMin:
                best, dMin = (poly, int(np.argmin(d))), d.min()
        return best

    def nearestEdge(self, cx, cy, maxDist=8):
        ''' (poly, index after which to insert, cx, cy) for the closest edge
        point within maxDist canvas pixels '''
        best, dMin = None, maxDist
        for poly in self.polys:
            px, py = self.polyCanvas(poly)
            qx, qy = np.roll(px, -1), np.roll(py, -1)
            ex, ey = qx - px, qy - py
            t = np.clip(((cx - px) * ex + (cy - py) * ey) /
                        np.maximum(ex**2 + ey**2, 1e-12), 0, 1)
            fx, fy = px + t * ex, py + t * ey
            d = np.hypot(fx - cx, fy - cy)
            i = int(np.argmin(d))
            if d[i] <= dMin:
                best, dMin = (poly, i, fx[i], fy[i]), d[i]
        return best

    def canvasToLL(self, cx, cy):
        x, y = self.canvasToXY(np.atleast_1d(cx), np.atleast_1d(cy))
        lat, lon = self.proj.xyToLL(x, y)
        return lat[0], lon[0]

    def vertexClick(self, cx, cy, button):
        if button == 1:
            self.drag = self.nearestVertex(cx, cy)
        elif button == 2:
            hit = self.nearestVertex(cx, cy)
            if hit is not None and len(hit[0]['lat']) > 3:
                poly, i = hit
                for key in ['lat', 'lon', 'z']:
                    poly[key] = np.delete(poly[key], i)
                self.modified = True
                self.drawPolys()
        elif button == 3:
            hit = self.nearestEdge(cx, cy)
            if hit is not None:
                poly, i, fx, fy = hit
                lat, lon = self.canvasToLL(fx, fy)
                z = demHeights(self.args.dem, [lat], [lon])[0]
                for key, val in zip(['lat', 'lon', 'z'], [lat, lon, z]):
                    poly[key] = np.insert(poly[key], i + 1, val)
                self.modified = True
                self.drawPolys()

    def onDrag(self, event):
        if self.mode != 'vertex' or self.drag is None:
            return
        poly, i = self.drag
        poly['lat'][i], poly['lon'][i] = self.canvasToLL(
            *self.eventCanvas(event))
        self.modified = True
        self.drawPolys()

    def onRelease(self, event):
        if self.mode != 'vertex' or self.drag is None:
            return
        poly, i = self.drag
        poly['z'][i] = demHeights(self.args.dem, poly['lat'][i:i + 1],
                                  poly['lon'][i:i + 1])[0]
        self.drag = None

    def toggleNearest(self, cx, cy):
        best, dMin = None, np.inf
        for poly in self.polys:
            px, py = self.polyCanvas(poly)
            d = np.min(np.hypot(px - cx, py - cy))
            if d < dMin:
                best, dMin = poly, d
        if best is not None:
            best['on'] = not best['on']
            self.modified = True
            self.drawPolys()

    def quit(self):
        if self.modified:
            from tkinter import messagebox
            answer = messagebox.askyesnocancel(
                'cleanpoly', 'Save edited velocity polygon file?')
            if answer is None:
                return
            self.save = answer
        self.root.destroy()


def main():
    args = cleanPolyArgs()
    proj = Projection(args.south)
    dirs = []
    for velDir, frameDir in getDirs(args.dirs):
        date1 = firstImageDate(velDir, args.velocityFile)
        if date1 is None or date1 < args.firstDate:
            print(f'skipping {frameDir}: first date {date1} before '
                  f'{args.firstDate}')
        elif args.noReplace and os.path.exists(f'{velDir}/{vertsName}'):
            print(f'skipping {velDir}: {vertsName} exists')
        else:
            dirs.append((velDir, frameDir))
    pool = ThreadPoolExecutor(max_workers=max(args.nCleanoff, 1))
    jobs = []
    for k, (velDir, frameDir) in enumerate(dirs):
        vertsFile = f'{velDir}/{vertsName}'
        print(f'\033[1m[{k + 1}/{len(dirs)}] {frameDir}\033[0m')
        if not args.indicesOnly:
            velData = readSpeed(velDir, args.velocityFile)
            if velData is None:
                print(f'{velDir}: no {args.velocityFile} velocity, skipping')
                continue
            polys = readVerts(vertsFile)
            editor = PolyEditor(velDir, polys, *velData, proj, args)
            if editor.run():
                writeVerts(vertsFile, polys)
                print(f'saved {vertsFile}')
            else:
                print('no edits saved')
        changed = makePolyList(frameDir, velDir)
        if changed and not args.noCleanoff:
            if os.path.exists(f'{frameDir}/cleanoff'):
                jobs.append(pool.submit(runCleanoff, frameDir))
            else:
                print(f'{frameDir}: no cleanoff script, not run')
    if len(jobs) > 0:
        print(f'waiting for {sum(not j.done() for j in jobs)} cleanoff(s)')
    pool.shutdown(wait=True)
    failed = [j.result()[0] for j in jobs if j.result()[1] != 0]
    print(f'cleanoff ran for {len(jobs)} frame(s), {len(failed)} failed')
    for frameDir in failed:
        print(f'  FAILED: {frameDir}')


if __name__ == '__main__':
    main()
