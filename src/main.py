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
        patches= []

        coordsmin = np.array([MAX, MAX], dtype=float)
        coordsmax = np.array([-MAX, -MAX], dtype=float)

        for cycle in cycles:
            for nodeid in cycle:
                aux = np.array([g.node[nodeid]['posx'], g.node[nodeid]['posy']])
                if aux[0] < coordsmin[0]: coordsmin[0] = aux[0]
                if aux[1] < coordsmin[1]: coordsmin[1] = aux[1]
                if aux[0] > coordsmax[0]: coordsmax[0] = aux[0]
                if aux[1] > coordsmax[1]: coordsmax[1] = aux[1]

        for j, cycle in enumerate(cycles):
            n = len(cycle)
            coords = np.ndarray((n, 2), dtype=float)
            for i, nodeid in enumerate(cycle):
                coords[i] = np.array([g.node[nodeid]['posx'], g.node[nodeid]['posy']])

            areas[j] = calculate_real_area(coords)

            coords = (coords - coordsmin) / (coordsmax-coordsmin)

            if areas[j] > epsilon:
                patches.append(Polygon(coords))
                f = coords

        p = PatchCollection(patches, alpha=0.1)
        colors = 100*np.random.rand(len(patches))
        p.set_array(np.array(colors))
        ax.add_collection(p)
        filename = os.path.splitext(filepath)[0]
        ax.set_title(filename)
        xticks = np.array(ax.get_xticks())
        yticks = np.array(ax.get_yticks())
        ticksrange = np.array([xticks[-1] - xticks[0], yticks[-1] - yticks[0]],
                              dtype=float)
        coordsrange = coordsmax - coordsmin
        factor = coordsrange / ticksrange
        xticks_new = xticks
        ax.set_xticklabels((xticks-xticks[0])*factor[0] + coordsmin[0])
        ax.set_yticklabels((xticks-xticks[0])*factor[1] + coordsmin[1])
        outvispath = pjoin(outdir,  filename + '.pdf')
        plt.savefig(outvispath)

        errorsind = areas < epsilon
        validind = areas > epsilon
        info('Number of valid {}, invalid {}'.format(np.sum(validind), np.sum(errorsind)))
        info('Block areas mean {:4f} ± {:4f}'.format(np.mean(areas[validind]),
                                                             np.std(areas[validind])))
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
    figs = {}
    for k in ['hist', 'hist500k', 'boxplot']:
        figs[k] = go.Figure()

    # patches = []
    # fig = go.Figure(data=[go.Histogram(x=areas[validind], nbinsx=500)])

    for k, areas in allareas.items():
        errorsind = areas < epsilon
        validind = areas > epsilon

        figs['hist'].add_trace(go.Histogram(x=areas[validind], histnorm='probability',
                                   # nbinsx=500,
                                   xbins=dict(size=1000000)
                                   ))
        figs['hist500k'].add_trace(go.Histogram(x=areas[areas<500000], histnorm='probability',
                                   # nbinsx=500,
                                   # xbins=dict(size=1000000)
                                   ))
        figs['boxplot'].add_trace(go.Box(y=areas[validind]))

    # figs['hist'].update_layout(barmode='overlay')
    figs['hist'].update_traces(opacity=0.75)
    plotly.offline.plot(figs['hist'], filename=pjoin(outdir, 'hist.html'), auto_open=False)
    plotly.offline.plot(figs['hist500k'], filename=pjoin(outdir, 'hist500k.html'), auto_open=False)
    plotly.offline.plot(figs['boxplot'], filename=pjoin(outdir, 'boxplot.html'), auto_open=False)
    # py.plot(data, filename='basic_line', auto_open=False)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('graphsdir', help='Graphs directory')
    parser.add_argument('--outdir', default='/tmp/', help='Output directory')
    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s] %(message)s',
    datefmt='%Y%m%d %H:%M', level=logging.INFO)
    # epsilon = 0.00000001 # In m2
    epsilon = 0.01 # In m2

    # if not os.path.exists(args.outdir):
    if True:
        os.makedirs(args.outdir, exist_ok=True)
        allareas = compute_block_areas(args.graphsdir, args.outdir, epsilon)
        dump_areas(allareas, args.outdir)

    return
    allareas = load_areas(args.outdir)
    plot_distributions(allareas, epsilon, args.outdir)

if __name__ == "__main__":
    main()

