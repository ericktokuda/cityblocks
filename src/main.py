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

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('graphsdir', help='Graphs directory')
    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s] %(message)s',
    datefmt='%Y%m%d %H:%M', level=logging.DEBUG)

    for filepath in os.listdir(args.graphsdir):
        g = nx.read_graphml(pjoin(args.graphsdir, filepath))
        cycles = nx.cycle_basis(g)
        for cycle in cycles:
            n = len(cycle)
            coords = np.ndarray((n, 2), dtype=float)
            for i, nodeid in enumerate(cycle):
                coords[i, 0] = g.node[nodeid]['posx']
                coords[i, 1] = g.node[nodeid]['posy']
            # area = polyarea(coords[:, 0], coords[:, 1])
            # print(area)

            import pyproj
            import shapely
            import shapely.ops as ops
            from shapely.geometry.polygon import Polygon
            from functools import partial

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

            print(geom_area.area)
    

if __name__ == "__main__":
    main()

