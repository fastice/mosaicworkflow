#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Build the frozen geomosaic -fillMask raster for an image-mosaic template.

The template's trackOrder tracks are tier 0 and each list in fillTrackOrder is
a fill tier (1, 2, ...). For every cell the fraction of nDays windows in which
each tier acquired it is computed from the date databases and the products'
geodat footprints. The mask value is the highest tier geomosaic may use there:

    0  off ice, where a tier-0 track covers >= threshold of the windows, or
       where no fill tier does (locked: only the trackOrder tracks are used)
    k  otherwise, the first fill tier k covering >= threshold

Each fill tier's region is then grown by -grow km into the locked tier-0
region, but never within -coastBuffer km of the ice margin. Footprints from
geodat corners overstate the valid data at swath edges by ~1 km, so without
this a thin sliver at the hole edge would be locked to tracks that do not
fill it. Tier 0 still wins wherever it has data.

Run once and keep the file fixed, so that every mosaic made with it uses the
same geometry at each pixel.
"""
import argparse
import json
import os
from collections import defaultdict
from datetime import datetime, timedelta
import numpy as np
from osgeo import gdal, ogr, osr
from scipy.ndimage import binary_dilation, distance_transform_edt
import mosaicfunc as mosf
import utilities as u

GIMPICEMASK = '/Volumes/insar7/ian/gimp/icerock/GimpIceMask_90m_2015_v1.2.tif'


def parseArgs():
    parser = argparse.ArgumentParser(
        description='Make the geomosaic -fillMask (max tier) raster for a '
        'template with trackOrder and fillTrackOrder')
    parser.add_argument('template', type=str, help='mosaic template yaml')
    parser.add_argument('maskFile', type=str, help='output GeoTIFF')
    parser.add_argument('-firstDate', '--firstDate', type=str, default=None,
                        help='YYYY-MM-DD [lastDate - 182 days]')
    parser.add_argument('-lastDate', '--lastDate', type=str, default=None,
                        help='YYYY-MM-DD [last date in the DBs]')
    parser.add_argument('-nDays', '--nDays', type=int, default=6,
                        help='mosaic window length in days [6]')
    parser.add_argument('-threshold', '--threshold', type=float, default=0.5,
                        help='fraction of windows that counts as covered [0.5]')
    parser.add_argument('-resolution', '--resolution', type=float,
                        default=500., help='mask cell size in m [500]')
    parser.add_argument('-iceMask', '--iceMask', type=str,
                        default=GIMPICEMASK, help='ice mask (1 = ice)')
    parser.add_argument('-grow', '--grow', type=float, default=2.,
                        help='grow fill regions into tier 0 by km [2]')
    parser.add_argument('-coastBuffer', '--coastBuffer', type=float,
                        default=30., help='never grow within km of the ice '
                        'margin [30]')
    args = parser.parse_args()
    return args


def readDates(dbPaths):
    ''' {date: [(dbPath, track, product), ...]} over all DBs '''
    myDB = mosf.sarDB()
    myDB.readDB(dbPaths)
    byDate = defaultdict(list)
    for dbPath, entries in myDB.byDateAndPath.items():
        for myDate, prods in entries.items():
            for prod in prods:
                byDate[myDate].append((dbPath, prod['track'],
                                       prod['product']))
    return byDate


def footprintWkt(dbPath, track, product, toPS):
    ''' Footprint polygon (WKT, template CRS) from the product geodat '''
    for geodat in ['geodat10x2.geojson', 'geodat30x6.geojson']:
        geoFile = f'{dbPath}/track-{track}/{product}/{geodat}'
        if os.path.exists(geoFile):
            with open(geoFile) as fp:
                corners = json.load(fp)['geometry']['coordinates'][0]
            # geodat corners are [lat, lon]
            xy = [toPS.TransformPoint(lat, lon)[0:2] for lat, lon in corners]
            return 'POLYGON ((' + ','.join(f'{x} {y}' for x, y in xy) + '))'
    return None


def rasterizeWkts(wkts, gridInfo):
    ''' 0/1 array of cells inside any of the polygons '''
    nx, ny, gt, srs = gridInfo
    mem = gdal.GetDriverByName('MEM').Create('', nx, ny, 1, gdal.GDT_Byte)
    mem.SetGeoTransform(gt)
    mem.SetProjection(srs.ExportToWkt())
    vds = ogr.GetDriverByName('Memory').CreateDataSource('')
    layer = vds.CreateLayer('fp', srs, ogr.wkbPolygon)
    for wkt in wkts:
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(ogr.CreateGeometryFromWkt(wkt))
        layer.CreateFeature(feature)
    gdal.RasterizeLayer(mem, [1], layer, burn_values=[1])
    return mem.GetRasterBand(1).ReadAsArray()


def tierCoverage(byDate, tierTracks, firstDate, lastDate, nDays, gridInfo,
                 toPS):
    '''
    For each tier, the max over its tracks of the fraction of nDays windows
    that acquired each cell.
    '''
    nWin = ((lastDate - firstDate).days + nDays) // nDays
    trackTier = {t: k for k, tracks in enumerate(tierTracks) for t in tracks}
    # (track, window) -> footprints
    windows = defaultdict(list)
    fpCache = {}
    myDate = firstDate
    while myDate <= lastDate:
        for dbPath, track, product in byDate.get(myDate.strftime('%Y-%m-%d'),
                                                 []):
            if track not in trackTier:
                continue
            key = (dbPath, track, product.split('_')[1])
            if fpCache.get(key) is None:
                fpCache[key] = footprintWkt(dbPath, track, product, toPS)
            if fpCache[key] is not None:
                windows[(track, (myDate - firstDate).days // nDays)].append(
                    fpCache[key])
        myDate += timedelta(days=1)
    nx, ny = gridInfo[0], gridInfo[1]
    counts = {t: np.zeros((ny, nx), np.int16) for t in trackTier}
    for (track, _), wkts in windows.items():
        counts[track] += rasterizeWkts(wkts, gridInfo)
    coverage = []
    for tracks in tierTracks:
        coverage.append(np.max([counts[t] for t in tracks], axis=0) / nWin)
        print(f'tier {len(coverage) - 1} {tracks}: max windows '
              f'{np.max([counts[t].max() for t in tracks])} of {nWin}')
    return coverage, nWin


def iceOnGrid(iceMask, gridInfo):
    ''' Ice where at least half the mask pixels in a cell are ice (averaged) '''
    nx, ny, gt, srs = gridInfo
    warped = gdal.Warp('', iceMask, format='MEM', dstSRS=srs.ExportToWkt(),
                       outputBounds=[gt[0], gt[3] + gt[5] * ny,
                                     gt[0] + gt[1] * nx, gt[3]],
                       width=nx, height=ny, resampleAlg='average',
                       outputType=gdal.GDT_Float32)
    return warped.GetRasterBand(1).ReadAsArray() >= 0.5


def main():
    args = parseArgs()
    gdal.UseExceptions()
    template = mosf.readTemplate(args.template)
    tierTracks = [[str(t) for t in template['trackOrder']]] + \
        [[str(t) for t in tier] for tier in template.get('fillTrackOrder', [])]
    if len(tierTracks) < 2 or '*' in tierTracks[0]:
        u.myerror('template needs an explicit trackOrder and a fillTrackOrder')
    # the same DBs setupimagemosaic reads
    dbPaths = mosf.mosaicTempArgs(args.template)[2]
    byDate = readDates(dbPaths)
    lastDate = datetime.strptime(args.lastDate or max(byDate), '%Y-%m-%d')
    firstDate = datetime.strptime(args.firstDate, '%Y-%m-%d') \
        if args.firstDate else lastDate - timedelta(days=182)
    print(f'{firstDate:%Y-%m-%d} to {lastDate:%Y-%m-%d}, {args.nDays} day '
          f'windows, threshold {args.threshold}')
    # Grid covering the template mosaic
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(template.get('epsg', 3413)))
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    x0 = template['xll'] - template['dx'] / 2
    y0 = template['yll'] - template['dy'] / 2
    res = args.resolution
    nx = int(np.ceil(template['sx'] * template['dx'] / res))
    ny = int(np.ceil(template['sy'] * template['dy'] / res))
    gt = (x0, res, 0., y0 + ny * res, 0., -res)
    gridInfo = (nx, ny, gt, srs)
    ll = osr.SpatialReference()
    ll.ImportFromEPSG(4326)
    ll.SetAxisMappingStrategy(osr.OAMS_AUTHORITY_COMPLIANT)  # lat, lon
    toPS = osr.CoordinateTransformation(ll, srs)
    #
    coverage, nWin = tierCoverage(byDate, tierTracks, firstDate, lastDate,
                                  args.nDays, gridInfo, toPS)
    ice = iceOnGrid(args.iceMask, gridInfo)
    locked = (coverage[0] >= args.threshold) | ~ice
    growable = ice & (distance_transform_edt(ice) * res / 1000. >=
                      args.coastBuffer)
    nGrow = int(round(args.grow * 1000. / res))
    disk = np.hypot(*np.mgrid[-nGrow:nGrow + 1, -nGrow:nGrow + 1]) <= nGrow
    maxTier = np.zeros((ny, nx), np.uint8)
    for k in range(len(tierTracks) - 1, 0, -1):
        fill = (coverage[k] >= args.threshold) & ~locked
        if nGrow > 0:
            fill |= binary_dilation(fill, structure=disk) & growable
        maxTier[fill] = k
    cellKm2 = res * res / 1e6
    for k in range(len(tierTracks)):
        print(f'maxTier {k}: {np.sum(ice & (maxTier == k)) * cellKm2:.0f} '
              'km2 of ice')
    # Write
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(args.maskFile, nx, ny, 1, gdal.GDT_Byte,
                       options=['COMPRESS=DEFLATE', 'TILED=YES'])
    ds.SetGeoTransform(gt)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(maxTier)
    ds.SetMetadata({'template': os.path.abspath(args.template),
                    'tiers': json.dumps(tierTracks),
                    'firstDate': f'{firstDate:%Y-%m-%d}',
                    'lastDate': f'{lastDate:%Y-%m-%d}',
                    'nDays': str(args.nDays), 'nWindows': str(nWin),
                    'threshold': str(args.threshold),
                    'grow': str(args.grow),
                    'coastBuffer': str(args.coastBuffer),
                    'iceMask': args.iceMask,
                    'created': datetime.now().strftime('%Y-%m-%d %H:%M')})
    ds = None
    print(f'wrote {args.maskFile}')


if __name__ == "__main__":
    main()
