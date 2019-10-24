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
def compute_block_areas(graphsdir, epsilon=0.00000001):
    """Calculate block areas from each graph in @graphsdir

    Args:
    graphsdir(str): path to the dir containing the graphml files

    Returns:
    dict: game of original file as key and areas as values
    """

    info('Reading directory {} ...'.format(graphsdir))
    allareas = {}
    for filepath in os.listdir(graphsdir):
        info('Reading file {} ...'.format(filepath))
        g = nx.read_graphml(pjoin(graphsdir, filepath))
        cycles = nx.cycle_basis(g)
        ncycles = len(cycles)
        info('Cycles found: {}'.format(ncycles))
        areas = np.ndarray(ncycles, dtype=float)
        for j, cycle in enumerate(cycles):
            n = len(cycle)
            coords = np.ndarray((n, 2), dtype=float)
            for i, nodeid in enumerate(cycle):
                coords[i] = np.array([g.node[nodeid]['posx'], g.node[nodeid]['posy']])

            areas[j] = calculate_real_area(coords)
            # res = polyarea(coords[:, 0], coords[:, 1])


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
def plot_distributions(allareas, epsilon):
    """Plot distributions for each array of areas

    Args:
    allareas(dict): filename as key and areas as values
    epsilon(float): value to consider as invalid area
    """

    import plotly.graph_objects as go
    fig = go.Figure()


    # patches = []
    # fig = go.Figure(data=[go.Histogram(x=areas[validind], nbinsx=500)])
    for k, areas in allareas.items():
        errorsind = areas < epsilon
        validind = areas > epsilon
        fig.add_trace(go.Histogram(x=areas[validind], nbinsx=500, histnorm='probability'))
        # break

        # fig.add_trace(go.Box(y=areas[validind]))

        # Plot scatter alltogether
        # fig = go.Figure(data=[go.Box(y=areas[validind],
                                     # boxpoints='all', # can also be outliers, or suspectedoutliers, or False
                                     # jitter=0.3, # add some jitter for a better separation between points
                                     # pointpos=-1.8 # relative position of points wrt box
                                     # )])



        # fig = go.Figure()
        # fig.add_trace(go.Scatter(y=areas[validind],
                                 # mode='lines',
                                 # name='lines'))
    fig.show()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('graphsdir', help='Graphs directory')
    parser.add_argument('--outdir', default='/tmp/', help='Output directory')
    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s] %(message)s',
    datefmt='%Y%m%d %H:%M', level=logging.INFO)
    epsilon = 0.00000001 # In m2

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        allareas = compute_block_areas(args.graphsdir, epsilon)
        dump_areas(allareas, args.outdir)

    allareas = load_areas(args.outdir)
    plot_distributions(allareas, epsilon)
if __name__ == "__main__":
    main()

