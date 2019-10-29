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

MAX = 999999999

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

def polyarea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def calculate_real_area(coords):
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
        return geom_area.area
    except Exception as e:
        return -1 # Invalid polygon
############################################################
def compute_block_areas(graphsdir, outdir, epsilon=0.00000001):
    """Calculate block areas from each graph in @graphsdir

    Args:
    graphsdir(str): path to the dir containing the graphml files

    Returns:
    dict: game of original file as key and areas as values
    """

    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Polygon
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)

    info('Reading directory {} ...'.format(graphsdir))
    allareas = {}
    for filepath in os.listdir(graphsdir):
        info('Reading file {} ...'.format(filepath))
        g = nx.read_graphml(pjoin(graphsdir, filepath))
        g = g.to_undirected()
        cycles = nx.cycle_basis(g)
        ncycles = len(cycles)
        info('Cycles found: {}'.format(ncycles))
        areas = np.ndarray(ncycles, dtype=float)
        # patches= []

        # coordsmin = np.array([MAX, MAX], dtype=float)
        # coordsmax = np.array([-MAX, -MAX], dtype=float)

        # for cycle in cycles:
            # for nodeid in cycle:
                # aux = np.array([g.node[nodeid]['posx'], g.node[nodeid]['posy']])
                # if aux[0] < coordsmin[0]: coordsmin[0] = aux[0]
                # if aux[1] < coordsmin[1]: coordsmin[1] = aux[1]
                # if aux[0] > coordsmax[0]: coordsmax[0] = aux[0]
                # if aux[1] > coordsmax[1]: coordsmax[1] = aux[1]

        for j, cycle in enumerate(cycles):
            n = len(cycle)
            coords = np.ndarray((n, 2), dtype=float)
            for i, nodeid in enumerate(cycle):
                coords[i] = np.array([g.node[nodeid]['posx'], g.node[nodeid]['posy']])

            areas[j] = calculate_real_area(coords)

            # coords = (coords - coordsmin) / (coordsmax-coordsmin)

            # if areas[j] > epsilon:
                # patches.append(Polygon(coords))
                # f = coords

        # p = PatchCollection(patches, alpha=0.1)
        # colors = 100*np.random.rand(len(patches))
        # p.set_array(np.array(colors))
        # ax.add_collection(p)
        # filename = os.path.splitext(filepath)[0]
        # ax.set_title(filename)
        # xticks = np.array(ax.get_xticks())
        # yticks = np.array(ax.get_yticks())
        # ticksrange = np.array([xticks[-1] - xticks[0], yticks[-1] - yticks[0]],
                              # dtype=float)
        # coordsrange = coordsmax - coordsmin
        # factor = coordsrange / ticksrange
        # xticks_new = xticks
        # ax.set_xticklabels((xticks-xticks[0])*factor[0] + coordsmin[0])
        # ax.set_yticklabels((xticks-xticks[0])*factor[1] + coordsmin[1])
        # outvispath = pjoin(outdir,  filename + '.pdf')
        # plt.savefig(outvispath)

        errorsind = areas < epsilon
        validind = areas > epsilon
        info('Number of valid {}, invalid {}'.format(np.sum(validind), np.sum(errorsind)))
        info('Block areas mean {:4f} ± {:4f}'.format(np.mean(areas[validind]),
                                                             np.std(areas[validind])))
        allareas[os.path.splitext(filepath)[0]] = areas
    return allareas
##########################################################
def dump_areas(allareas, outdir):
    """Dump @areas in outdir

    Args:
    areas(dict): filename (without ext) as keys and areas as values
    outdir(str): output directory
    """

    for filename, areas in allareas.items():
        outpath = pjoin(outdir, filename + '.csv')
        np.savetxt(outpath, areas, delimiter=',',
                   header='aream2', comments='')

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
def plot_distributions(allareas, epsilon, outdir):
    """Plot distributions for each array of areas

    Args:
    allareas(dict): filename as key and areas as values
    epsilon(float): value to consider as invalid area
    """

    import plotly.graph_objects as go
    import plotly
    from scipy import stats
    figs = {}
    for k in ['hist', 'boxplot']:
        figs[k] = go.Figure()

    for k, areas in allareas.items():
        errorsind = areas < epsilon
        validind = areas > epsilon
        figs['hist'].add_trace(
                go.Histogram(x=stats.zscore(areas[validind]),
                    histnorm='probability',
                    name=k,
                    # nbinsx=500,
                    # xbins=dict(size=1000000)
                    xbins=dict(size=.1)
                    ))
        figs['boxplot'].add_trace(
                go.Box(y=areas[validind],
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
        info(filepath)
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
        # print(lab)
        z = np.where(labels == lab)
        labeled_img[z] = _colors[i]

    return labeled_img

##########################################################

def get_components_from_raster(rasterdir, compdir):
    if os.path.exists(compdir):
        info('{} already exists. Skipping ...'.format(compdir))
        return pkl.load(open(pjoin(compdir, 'components.pkl'), 'rb'))

    os.makedirs(compdir)

    info('Reading raster graphics from {} ...'.format(rasterdir))
    components = {}

    for filepath in os.listdir(rasterdir):
        if not filepath.endswith('.png'): continue
        info(filepath)
        imgpath = pjoin(rasterdir, filepath)
        img = cv2.imread(imgpath, 0)
        img = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)[1]

        # kernel = np.ones((3,3),np.uint8)
        # img = cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel)

        ret, labels = cv2.connectedComponents(img)
        components[os.path.splitext(filepath)[0]] = labels
    pkl.dump(components, open(pjoin(compdir, 'components.pkl'), 'wb'))
    return components

##########################################################
def generate_components_vis(components, compdir):
    info('Generating components visualization ...')

    for k, labels in components.items():
        labeled_img = colorize(labels)
        outpath = pjoin(compdir, k + '_labels.png')
        cv2.imwrite(outpath, labeled_img)

##########################################################
def calculate_block_areas(labels, outdir):
    outpath = pjoin(outdir, 'areas.pkl')

    if os.path.exists(outpath):
        info('{} already exists. Skipping ...'.format(outpath))
        return pkl.load(open(outpath, 'rb'))

    info('Computing block areas from components ...')

    areas = {}
    for k, label in labels.items():
        info(k)

        ulabels = np.unique(label)
        areas[k] = np.ndarray(len(ulabels), dtype=int)

        for j, ulabel in enumerate(ulabels):
            areas[k][j] = np.sum(np.where(label == ulabel))

    pkl.dump(areas, open(outpath, 'wb'))
    return areas

def filter_areas(areas):
    for k, v in areas.items():
        areas[k] = v[2:]
    return areas
##########################################################

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('graphsdir', help='Graphs directory')
    parser.add_argument('--outdir', default='/tmp/', help='Output directory')
    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s] %(message)s',
    datefmt='%Y%m%d %H:%M', level=logging.INFO)
    epsilon = 0.01 # In m2

    # if not os.path.exists(args.outdir):
        # os.makedirs(args.outdir, exist_ok=True)
        # allareas = compute_block_areas(args.graphsdir, args.outdir, epsilon)
        # dump_areas(allareas, args.outdir)

    # allareas = load_areas(args.outdir)
    # plot_distributions(allareas, epsilon, args.outdir)

    skeldir = pjoin(args.outdir, 'skel')
    compdir = pjoin(args.outdir, 'comp')

    plot_graph_raster(args.graphsdir, skeldir)
    components = get_components_from_raster(skeldir, compdir)
    # generate_components_vis(components, compdir)
    areas = calculate_block_areas(components, args.outdir)
    areas = filter_areas(areas) # 0: skeleton, 1: background
    plot_distributions(areas, epsilon, args.outdir)

if __name__ == "__main__":
    main()

