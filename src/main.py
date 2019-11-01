#!/usr/bin/env python3
"""Analyze block sizes
"""

import argparse
import logging
import os
from os.path import join as pjoin
from logging import debug, info
import numpy as np
import igraph
import networkx as nx
import matplotlib.pyplot as plt
import cv2
import pickle as pkl

# for function calculate_real_area
import pyproj
import shapely
import shapely.ops as ops
from shapely.geometry.polygon import Polygon
from functools import partial
import plotly.graph_objects as go
import plotly
from haversine import haversine
from scipy import stats

MAX = 999999999

##########################################################
def xnet2igraph_batch(xnetdir):
    """Convert dir containing xnet graphs

    Args:
    xnetdir(str): path to the xnet graphs
    """
    # https://github.com/filipinascimento/xnet/blob/master/xnetwork/__init__.py
    import xnet
    outdir = '/tmp/'
    xnetdir = '/tmp/'
    acc = 0
    for xnetgraphpath in os.listdir(xnetdir):
        if not xnetgraphpath.endswith('.xnet'): continue
        g = xnet.xnet2igraph(pjoin(xnetdir, xnetgraphpath))
        coords = [ (x,y) for x, y in zip(g.vs['posx'], g.vs['posy']) ]
        outfilename = pjoin(outdir, os.path.splitext(xnetgraphpath)[0] + '.graphml')
        igraph.write(g, outfilename, 'graphml')
        acc += 1
    info('Sucessfully converted {} xnet files.'.format(acc))

##########################################################
def polyarea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

##########################################################
def calculate_real_area(coords, unit='km'):
    # https://gis.stackexchange.com/questions/127607/area-in-km-from-polygon-of-coordinates
    try:
        geom = Polygon(coords)
        geom_area = ops.transform(
            partial(
                pyproj.transform,
                pyproj.Proj(init='EPSG:4326'),
                pyproj.Proj(
                    proj='aea',
                    lat1=geom.bounds[1],
                    lat2=geom.bounds[3])),
            geom)
        if unit == 'km': factor = 1000000
        else: factor = 1
        return geom_area.area / factor
    except Exception as e:
        return -1 # Invalid polygon

##########################################################
def load_areas(outdir):
    """Load @areas from @outdir

    Args:
    allareas(dict): filename (without ext) as keys and areas as values
    outdir(str): output directory
    """

    allareas = {}
    for f in os.listdir(outdir):
        if not f.endswith('.csv'): continue

        k = os.path.splitext(f)[0]
        csvpath = pjoin(outdir, f)
        allareas[k] = np.loadtxt(csvpath, dtype=float, skiprows=1)

    return allareas

##########################################################
def plot_distributions(allareas, outdir):
    """Plot distributions for each array of areas

    Args:
    allareas(dict): filename as key and areas as values
    epsilon(float): value to consider as invalid area
    """

    figs = {}
    for k in ['hist', 'boxplot']:
        figs[k] = go.Figure()

    for k, areas in allareas.items():
        figs['hist'].add_trace(
                go.Histogram(x=stats.zscore(areas),
                    histnorm='probability',
                    name=k,
                    # nbinsx=500,
                    # xbins=dict(size=1000000)
                    xbins=dict(size=.1)
                    ))
        figs['boxplot'].add_trace(
                go.Box(y=areas,
                    name=k,
            ))

    figs['hist'].update_traces(opacity=0.75)
    plotly.offline.plot(figs['hist'], filename=pjoin(outdir, 'hist.html'), auto_open=False)
    plotly.offline.plot(figs['boxplot'], filename=pjoin(outdir, 'boxplot.html'), auto_open=False)

#########################################################
def plot_graph_raster(graphsdir, skeldir):
    if os.path.exists(skeldir):
        info('{} already exists. Skipping ...'.format(skeldir))
        return

    os.makedirs(skeldir)

    info('Reading graphs from {} ...'.format(graphsdir))
    allareas = {}
    figfactor = 10000
    for filepath in os.listdir(graphsdir):
        if not filepath.endswith('.graphml'): continue
        info(' *' + filepath)
        g = igraph.Graph.Read(pjoin(graphsdir, filepath))
        lonrange = [np.min(g.vs['posx']), np.max(g.vs['posx'])]
        latrange = [np.min(g.vs['posy']), np.max(g.vs['posy'])]

        lonfigsize = (lonrange[1] - lonrange[0])*figfactor
        latfigsize = (latrange[1] - latrange[0])*figfactor

        visual = dict(
            vertex_size = 0,
            vertex_color = [0, 0, 0, 0],
            edge_color = [0, 0, 0, 0],
            edge_width = 2,
            bbox = (lonfigsize, latfigsize)
        )

        outpath = pjoin(skeldir, os.path.splitext(filepath)[0] + '.png')
        layout = [ (x, -y) for x, y in zip(g.vs['posx'], g.vs['posy']) ]
        igraph.plot(g, target=outpath, layout=layout, **visual)

##########################################################
def colorize(labels):
    ulabels = np.unique(labels)
    _colors = np.random.randint(0, 255, (len(ulabels), 3))

    labeled_img = np.zeros((labels.shape[0], labels.shape[1], 3), dtype=np.uint8)
    for i, lab in enumerate(ulabels[1:]): # 0: skeleton, 2: external part
        z = np.where(labels == lab)
        labeled_img[z] = _colors[i]

    return labeled_img

##########################################################
def get_components_from_raster(rasterdir, outdir):
    outpath = pjoin(outdir, 'components.pkl')

    if os.path.exists(outpath):
        info('{} already exists. Skipping ...'.format(outdir))
        return pkl.load(open(outpath, 'rb'))

    info('Reading raster graphics from {} ...'.format(rasterdir))
    components = {}

    for filepath in os.listdir(rasterdir):
        if not filepath.endswith('.png'): continue
        info(' *' + filepath)
        imgpath = pjoin(rasterdir, filepath)
        img = cv2.imread(imgpath, 0)
        img = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)[1]

        # kernel = np.ones((3,3),np.uint8)
        # img = cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel)

        ret, labels = cv2.connectedComponents(img)
        components[os.path.splitext(filepath)[0]] = labels
    pkl.dump(components, open(outpath, 'wb'))
    return components

##########################################################
def generate_components_vis(components, compdir):
    info('Generating components visualization ...')
    if os.path.exists(compdir):
        info('{} already exists. Skipping ...'.format(compdir))
        return

    os.makedirs(compdir)

    for k, labels in components.items():
        info(' *' + k)
        labeled_img = colorize(labels)
        outpath = pjoin(compdir, k + '_labels.png')
        cv2.imwrite(outpath, labeled_img)

##########################################################
def compute_raster_real_conversion(rasterranges, lonlatranges):
    conversionfactors = {}
    for k, c in lonlatranges.items():
        coords = np.array([ [c[0], c[1]], [c[0], c[3]], [c[2], c[3]], [c[2], c[1]], ])
        real = calculate_real_area(coords)
        r = rasterranges[k]
        raster = (r[2] - r[0]) * (r[3] - r[1])
        conversionfactors[k] = real/raster
    return conversionfactors

##########################################################
def calculate_block_areas(labels, lonlatranges, outdir):
    outpath = pjoin(outdir, 'areas.pkl')
    if os.path.exists(outpath):
        info('{} already exists. Skipping ...'.format(outpath))
        return pkl.load(open(outpath, 'rb'))

    areas, rasterranges = calculate_raster_areas(labels, outdir)
    conversionfactors = compute_raster_real_conversion(rasterranges, lonlatranges)

    for k in areas.keys():
        areas[k] = areas[k] * conversionfactors[k]
    pkl.dump(areas, open(outpath, 'wb'))
    return areas

##########################################################
def calculate_raster_areas(labels, outdir):
    info('Computing block areas from components ...')

    areas = {}
    ranges = {}
    for k, label in labels.items():
        info(' *' + k)
        ulabels, area = np.unique(label, return_counts=True)
        areas[k] = np.array(area)
        skelids = np.where(label == 1) # 1: skeleton (cv2)
        ranges[k] = np.array([ np.min(skelids[0]),  np.min(skelids[1]),
            np.max(skelids[0]), np.max(skelids[1]) ])
    return areas, ranges

##########################################################
def filter_areas(areas):
    for k, v in areas.items():
        areas[k] = v[2:]
    return areas

##########################################################
def compute_graph_statistics(graphsdir):
    info('Computing graph statistics ...')
    segmean = {}
    segstd = {}
    udistmean = {}
    udiststd = {}
    wdistmean = {}
    wdiststd = {}

    for filepath in os.listdir(graphsdir):
        if not filepath.endswith('.graphml'): continue

        k = os.path.splitext(filepath)[0]
        info(' *' + filepath)
        g = igraph.Graph.Read(pjoin(graphsdir, filepath))
        d = np.array(g.shortest_paths())
        wd = np.array(g.shortest_paths(weights=g.es['weight']))
        nvertices = len(g.vs)

        # avgdist = 0.0
        # avgwdist = 0.0
        ndists = int((nvertices * (nvertices-1)) / 2)

        acc = 0
        inds = np.ndarray((ndists, 2), dtype=int)

        for i in range(nvertices):
            for j in range(i+1, nvertices):
                inds[acc] = [i, j]
                # avgdist += d[i, j] / ndists
                # avgwdist += d[i, j] / ndists
                acc += 1

        udistmean[k] = np.mean(d[inds[:, 0], inds[:, 1]])
        udiststd[k] = np.std(d[inds[:, 0], inds[:, 1]])
        wdistmean[k] = np.mean(wd[inds[:, 0], inds[:, 1]])
        wdiststd[k] = np.std(wd[inds[:, 0], inds[:, 1]])
        segmean[k] = np.mean(g.es['weight'])
        segstd[k] = np.std(g.es['weight'])
    return segmean, segstd, udistmean, udiststd, wdistmean, wdiststd

##########################################################
def compute_statistics(graphsdir, allareas, outdir):
    """Compute statistics from areas

    Args:
    areas(dict): city as keys and list of areas as values
    """
    outpath = pjoin(outdir, 'summary.csv')
    if os.path.exists(outpath):
        info('{} already exists. Skipping ...'.format(outpath))
        return

    fh = open(outpath, 'w')
    fh.write('city,nblocks,areamean,areastd,areacv,areamin,areamax,areaentropy,areaeveness,segmean,segstd,udistmean,udiststd,wdistmean,wdiststd\n')

    segmean, segstd, udistmean, udiststd, wdistmean, wdiststd = compute_graph_statistics(graphsdir)
    for k, areas in allareas.items():
        a = areas[2:] # 0: skeleton, 1: background
        arel = a / np.sum(a)
        entropy = -np.sum(arel * np.log(arel))
        evenness = np.exp(entropy) / len(arel)
        st = [k, len(a), np.mean(a), np.std(a), np.std(a)/np.mean(a),
                np.min(a), np.max(a), entropy, evenness,
                segmean[k], segstd[k],
                udistmean[k], udiststd[k], wdistmean[k], wdiststd[k]
                ]
        fh.write(','.join([ str(s) for s in st]) + '\n')

##########################################################
def generate_test_graphs(outdir):
    sz = 20 
    g = igraph.Graph.Lattice([sz, sz], circular=False) # lattice
    coords = np.array(g.layout('grid').coords)
    _range = (np.max(coords, 0) - np.min(coords, 0)) * 10
    g.vs['posx'] = coords[:, 0] / _range[0]
    g.vs['posy'] = coords[:, 1] / _range[1]
    outfilename = pjoin(outdir, 'lattice.graphml')
    igraph.write(g, outfilename, 'graphml')

    g = igraph.Graph.Erdos_Renyi(60, p=1) # lattice
    coords = np.array(g.layout('random').coords)
    _range = (np.max(coords, 0) - np.min(coords, 0)) * 10
    g.vs['posx'] = coords[:, 0] / _range[0]
    g.vs['posy'] = coords[:, 1] / _range[1]
    outfilename = pjoin(outdir, 'erdos.graphml')
    igraph.write(g, outfilename, 'graphml')

##########################################################
def add_weights_to_edges(graphsdir, weightdir):
    if os.path.exists(weightdir):
        info('{} already exists. Skipping ...'.format(weightdir))
        return
    os.makedirs(weightdir)

    for filepath in os.listdir(graphsdir):
        if not filepath.endswith('.graphml'): continue
        outpath = pjoin(weightdir, filepath)
        info(' *' + filepath)
        g = igraph.Graph.Read(pjoin(graphsdir, filepath))
        for e in g.es:
            coordu = np.array([ g.vs[e.source]['posy'], g.vs[e.source]['posx'] ])
            coordv = np.array([ g.vs[e.target]['posy'], g.vs[e.target]['posx'] ])
            e['weight'] = haversine(coordu, coordv, unit='km') # in meters
        igraph.write(g, outpath, 'graphml')

##########################################################
def get_maps_ranges(graphsdir):
    info('Getting map ranges ...')
    ranges = {}
    for filepath in os.listdir(graphsdir):
        k = os.path.splitext(filepath)[0]
        if not filepath.endswith('.graphml'): continue
        info(' *' + filepath)
        g = igraph.Graph.Read(pjoin(graphsdir, filepath))
        lon = g.vs['posx']
        lat = g.vs['posy']
        ranges[k] = np.array([np.min(lon), np.min(lat), np.max(lon), np.max(lat)])
    return ranges

##########################################################
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('graphsdir', help='Graphs directory')
    parser.add_argument('--outdir', default='/tmp/', help='Output directory')
    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s] %(message)s',
    datefmt='%Y%m%d %H:%M', level=logging.INFO)

    skeldir = pjoin(args.outdir, 'skel')
    compdir = pjoin(args.outdir, 'comp')
    weightdir = pjoin(args.outdir, 'weighted')

    # generate_test_graphs(args.graphsdir)

    add_weights_to_edges(args.graphsdir, weightdir)
    lonlatranges = get_maps_ranges(args.graphsdir)
    plot_graph_raster(args.graphsdir, skeldir)
    components = get_components_from_raster(skeldir, args.outdir)
    generate_components_vis(components, compdir)
    allareas = calculate_block_areas(components, lonlatranges, args.outdir)
    compute_statistics(weightdir, allareas, args.outdir)
    filteredareas = filter_areas(allareas) # 0: skeleton, 1: background
    plot_distributions(filteredareas, args.outdir)

if __name__ == "__main__":
    main()

