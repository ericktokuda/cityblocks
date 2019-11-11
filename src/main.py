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
# import networkx as nx
# import matplotlib.pyplot as plt
import cv2
import pickle as pkl
import pandas as pd
# import graph_tool
# import graph_tool.centrality
# import graph_tool.topology

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
from scipy.stats.stats import pearsonr

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
        g.vs['x'] = g.vs['posx']
        g.vs['y'] = g.vs['posy']
        del g.vs['posx']
        del g.vs['posy']
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
def plot_distributions2(outdir):
    """Plot distributions for each array of areas

    Args:
    allareas(dict): filename as key and areas as values
    epsilon(float): value to consider as invalid area
    """

    figs = {}
    for k in ['entropylog', 'entropycvlog', 'entropyhist']:
        figs[k] = go.Figure()

    ########################################################## Area plots
    areaspath = pjoin(outdir, 'areas.pkl')
    allareas = pkl.load(open(areaspath, 'rb'))
    allareas = filter_areas(allareas) # 0: skeleton, 1: background

    df = pd.read_csv(pjoin(outdir, 'summary.csv'))


    ########################################################## Entropy plots
    for _, row in df.iterrows():
         figs['entropylog'].add_trace(go.Scatter(
            x=[row.areaentropy],
            y=[row.wdistmean/row.segmean],
            mode='markers',
            marker_size=20,
            name=row.city,
            ))
    entropylogpearson = pearsonr(df.areaentropy, np.log(df.wdistmean/df.segmean))
    figs['entropylog'].update_layout(
            title="Entropy of areas vs. average distance (pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Entropy of block area",
            yaxis_title="Average distance",
            yaxis_type="log"
            )

    figs['entropyhist'].add_trace(
            go.Histogram(x=df.areaentropy,
                histnorm='probability',
                ))
    figs['entropyhist'].update_layout(
            title="Entropy of block area",
            xaxis_title="Entropy",
            yaxis_title="Relative frequency",
            )

    ########################################################## Evenness cv plots
    for _, row in df.iterrows():
        figs['entropycvlog'].add_trace(go.Scatter(
            x=[row.areaentropy],
            y=[row.wdiststd/row.wdistmean/row.segmean],
            mode='markers',
            marker_size=20,
            name=row.city,
            ))
    entropycvlogpearson = pearsonr(df.areaentropy, np.log(df.wdiststd/df.wdistmean/df.segmean))
    figs['entropycvlog'].update_layout(
            title="Entropy of areas vs. coefficient of variance of the distance (pearson: {:.2f})".\
                    format(entropycvlogpearson[0]),
            xaxis_title="Entropy of block area",
            yaxis_title="Coefficient of variance of the distance",
            yaxis_type='log',
            )
    ########################################################## Save to file
    plotly.offline.plot(figs['entropylog'], filename=pjoin(outdir, 'entropylog.html'), auto_open=False)
    plotly.offline.plot(figs['entropyhist'], filename=pjoin(outdir, 'entropyhist.html'), auto_open=False)
    plotly.offline.plot(figs['entropycvlog'], filename=pjoin(outdir, 'entropycvlog.html'), auto_open=False)

    figs['entropylog'].write_image("/tmp/foo.png")

##########################################################
def plot_distributions(outdir):
    """Plot distributions for each array of areas

    Args:
    allareas(dict): filename as key and areas as values
    epsilon(float): value to consider as invalid area
    """

    figs = {}
    for k in ['hist', 'boxplot', 'entropy', 'evenness', 'entropycv', 'entropyz', 'entropylog', 'entropyhist',
            'coefvar', 'coefvarexp', 'betweenness']:
        figs[k] = go.Figure()

    ########################################################## Area plots
    areaspath = pjoin(outdir, 'areas.pkl')
    allareas = pkl.load(open(areaspath, 'rb'))
    allareas = filter_areas(allareas) # 0: skeleton, 1: background


    for k, areas in allareas.items():
        figs['hist'].add_trace(
                # go.Histogram(x=stats.zscore(areas),
                go.Histogram(x=areas,
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
    figs['hist'].update_layout(
            title="Area sizes",
            xaxis_title="Area (km^2)",
            yaxis_title="Relative frequency",
            )
    df = pd.read_csv(pjoin(outdir, 'summary.csv'))

    ########################################################## Entropy plots
    for _, row in df.iterrows():
        figs['entropy'].add_trace(go.Scatter(
            x=[row.areaentropy],
            y=[row.wdistmean/row.segmean],
            mode='markers',
            marker_size=20,
            name=row.city,
            ))
        figs['entropylog'].add_trace(go.Scatter(
            x=[row.areaentropy],
            y=[row.wdistmean/row.segmean],
            mode='markers',
            marker_size=20,
            name=row.city,
            ))
        figs['entropyz'].add_trace(go.Scatter(
            x=[row.areaentropy],
            y=[(row.wdistmean-np.mean(df.wdistmean)) / np.std(df.wdiststd)],
            mode='markers',
            marker_size=20,
            name=row.city,
            ))
    figs['entropy'].update_layout(
            title="Entropy of areas X mean distance (normalized)",
            xaxis_title="Entropy of block area",
            yaxis_title="Shortest path length (normalized by segment length mean)",
            )
    entropylogpearson = pearsonr(df.areaentropy, np.log(df.wdistmean/df.segmean))
    figs['entropylog'].update_layout(
            title="Entropy of areas X mean distance (normalized) pearson: " + str(entropylogpearson),
            xaxis_title="Entropy of block area",
            yaxis_title="Shortest path length (normalized by segment length mean) - LOG",
            yaxis_type="log"
            )
    figs['entropyz'].update_layout(
            title="Entropy of areas X mean distance (z-score)",
            xaxis_title="Entropy of block area",
            yaxis_title="Shortest path length (standardized among all measures)",
            )

    figs['entropyhist'].add_trace(
            go.Histogram(x=df.areaentropy,
                histnorm='probability',
                # name=k,
                # xbins=dict(size=.1)
                ))

    ########################################################## Evenness plots
    for _, row in df.iterrows():
        figs['evenness'].add_trace(go.Scatter(
            x=[row.areaeveness],
            y=[row.wdistmean/row.segmean],
            mode='markers',
            marker_size=40,
            name=row.city,
            ))
    figs['evenness'].update_layout(
            title="Evenness of areas X mean distance (normalized)",
            xaxis_title="Evenness of block area",
            yaxis_title="Shortest path length (normalized by segment length mean)",
            )

    ########################################################## Evenness cv plots
    for _, row in df.iterrows():
        figs['entropycv'].add_trace(go.Scatter(
            x=[row.areaentropy],
            y=[row.wdiststd/row.wdistmean/row.segmean],
            mode='markers',
            marker_size=40,
            name=row.city,
            ))
    figs['entropycv'].update_layout(
            title="Entropy of areas X CV of mean distance (normalized)",
            xaxis_title="Entropy of block area",
            yaxis_title="Coefficient of dispersion of the shortest path length (normalized by segment length mean)",
            )
    ########################################################## Coefficient of variation
    for _, row in df.iterrows():
        figs['coefvar'].add_trace(go.Scatter(
            x=[row.areamean/row.areasum],
            y=[row.wdistmean/row.segmean],
            mode='markers',
            marker_size=40,
            name=row.city,
            ))
    figs['coefvar'].update_layout(
            title="Length of shortest paths X Rel. size of blocks",
            xaxis_title="Average relative size of blocks",
            yaxis_title="Length of shortest paths normalized by average segment length",
            )
    ########################################################## Coefficient of variation exp
    for _, row in df.iterrows():
        figs['coefvarexp'].add_trace(go.Scatter(
            x=[row.areamean/row.areasum],
            y=[row.wdistmean/row.segmean],
            mode='markers',
            marker_size=40,
            name=row.city,
            ))
    figs['coefvarexp'].update_layout(
            title="Length of shortest paths X Rel. size of blocks",
            xaxis_title="Average relative size of blocks",
            yaxis_title="Length of shortest paths normalized by average segment length",
            )
    ########################################################## betweenness
    for _, row in df.iterrows():
        figs['betweenness'].add_trace(go.Scatter(
            x=[row.areamean/row.areasum],
            y=[row.betwvmean],
            mode='markers',
            marker_size=40,
            name=row.city,
            ))
    figs['betweenness'].update_layout(
            title="Betweenness x coefficience of variantion",
            xaxis_title="Coefficient of variation",
            yaxis_title="Average betweenness",
            )
    ########################################################## Save to file
    plotly.offline.plot(figs['hist'], filename=pjoin(outdir, 'hist.html'), auto_open=False)
    plotly.offline.plot(figs['boxplot'], filename=pjoin(outdir, 'boxplot.html'), auto_open=False)
    plotly.offline.plot(figs['entropy'], filename=pjoin(outdir, 'entropy.html'), auto_open=False)
    plotly.offline.plot(figs['entropylog'], filename=pjoin(outdir, 'entropylog.html'), auto_open=False)
    plotly.offline.plot(figs['entropyz'], filename=pjoin(outdir, 'entropyz.html'), auto_open=False)
    plotly.offline.plot(figs['entropyhist'], filename=pjoin(outdir, 'entropyhist.html'), auto_open=False)
    plotly.offline.plot(figs['evenness'], filename=pjoin(outdir, 'evenness.html'), auto_open=False)
    plotly.offline.plot(figs['entropycv'], filename=pjoin(outdir, 'entropycv.html'), auto_open=False)
    plotly.offline.plot(figs['coefvar'], filename=pjoin(outdir, 'coefvar.html'), auto_open=False)
    plotly.offline.plot(figs['coefvarexp'], filename=pjoin(outdir, 'coefvarexp.html'), auto_open=False)
    plotly.offline.plot(figs['betweenness'], filename=pjoin(outdir, 'betweenness.html'), auto_open=False)
    figs['entropylog'].write_image("/tmp/foo.pdf")

#########################################################
def plot_graph_raster(graphsdir, skeldir):
    if os.path.exists(skeldir):
        info('{} already exists. Skipping ...'.format(skeldir))
        return

    os.makedirs(skeldir)

    info('Reading graphs from {} ...'.format(graphsdir))
    allareas = {}
    figfactor = 20000
    for filepath in os.listdir(graphsdir):
        if not filepath.endswith('.graphml'): continue
        info(' *' + filepath)
        g = igraph.Graph.Read(pjoin(graphsdir, filepath))
        for ed in g.es():
            if 'geometry' not in ed.attributes(): continue # split lines
            if ed['geometry'] == None or ed['geometry'] == '': continue

            y = ed['geometry'].split('(')[1][:-1].split(', ')
            artnodes = []
            nartnodes = len(y)
            nnodes = len(g.vs)

            for j in range(nartnodes):
                aux = y[j].split(' ')
                g.add_vertex(nnodes + j, x=float(aux[0]), y=float(aux[1]))

            g.add_edge(g.vs[ed.source], g.vs[nnodes]) #source->first
            for j in range(1, nartnodes): #internal nodes
                g.add_edge(g.vs[nnodes+j-1], g.vs[nnodes+j])
            g.add_edge(g.vs[nnodes+nartnodes-1], g.vs[ed.target]) #internal->target
            g.delete_edges((ed.source, ed.target))

        g.to_undirected()
        lonrange = [np.min(g.vs['x']), np.max(g.vs['x'])]
        latrange = [np.min(g.vs['y']), np.max(g.vs['y'])]

        lonfigsize = (lonrange[1] - lonrange[0])*figfactor
        latfigsize = (latrange[1] - latrange[0])*figfactor

        COL='#000000'
        # COL='#4e7f80'
        visual = dict(
            vertex_size = 0,
            # vertex_size = 5,
            vertex_frame_width = 0,
            vertex_color = COL,
            edge_color = COL,
            edge_width = 2,
            bbox = (lonfigsize, latfigsize)
        )

        outpath = pjoin(skeldir, os.path.splitext(filepath)[0] + '.png')
        layout = [ (x, -y) for x, y in zip(g.vs['x'], g.vs['y']) ]
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

        ret, labels = cv2.connectedComponents(img, connectivity=4)
        components[os.path.splitext(filepath)[0]] = labels
    pkl.dump(components, open(outpath, 'wb'))
    return components

##########################################################
def generate_components_vis(components, compdir):
    if os.path.exists(compdir):
        info('{} already exists. Skipping ...'.format(compdir))
        return

    info('Generating components visualization ...')
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
def compute_statistics(graphsdir, allareas, outdir):
    """Compute statistics from areas

    Args:
    areas(dict): city as keys and list of areas as values
    """
    outpath = pjoin(outdir, 'summary.csv')
    if os.path.exists(outpath):
        info('{} already exists. Skipping ...'.format(outpath))
        return

    info('Computing graph statistics...')
    fh = open(outpath, 'w', buffering=1)
    fh.write('city,nvertices,nedges,nblocks,areasum,areamean,areastd,areacv,areamin,areamax,' \
            'areaentropy,areaeveness,segmean,segstd,udistmean,udiststd,' \
            'wdistmean,wdiststd,betwvmean,betwvstd\n')

    segmean = {}
    segstd = {}
    udistmean = {}
    udiststd = {}
    wdistmean = {}
    wdiststd = {}
    nvertices = {}
    betwvmean = {}
    betwvstd = {}
    betwemean = {}
    betwestd = {}
    nedges = {}
    errors = []

    for k in allareas.keys():
        info(' *' + k)
        try:
            filepath = pjoin(graphsdir, k + '.graphml')
            g = igraph.Graph.Read(filepath)
            # gt = graph_tool.load_graph(filepath)

            n = len(g.vs)
            ndists = int((n * (n-1)) / 2)

            udsum = 0.0
            ud2sum = 0.0
            wdsum = 0.0
            wd2sum = 0.0

            ########################################################## slower
            # dists = graph_tool.topology.shortest_distance(gt)
            # wdists = graph_tool.topology.shortest_distance(gt, weights=gt.edge_properties['weight'])
            # for i in range(n):
                # for j in range(i+1, n):
                    # udsum += dists[i][j]
                    # ud2sum += dists[i][j]*dists[i][j]

                    # wdsum += wdists[i][j]
                    # wd2sum += wdists[i][j]*wdists[i][j]
            ##########################################################

            # Using the function for all vertice at once crashes
            for i in range(n):
                aux = np.array(g.shortest_paths(source=i, mode='ALL'))[0][i+1:]
                udsum += np.sum(aux)
                ud2sum += np.sum(np.square(aux))

                aux = np.array(g.shortest_paths(source=i, mode='ALL',
                    weights=g.es['weight']))[0][i+1:]
                wdsum += np.sum(aux)
                wd2sum += np.sum(np.square(aux))

            segmean[k] = np.mean(g.es['weight'])
            segstd[k] = np.std(g.es['weight'])
            
            # bv, be = graph_tool.centrality.betweenness(gt)
            bv = g.betweenness()
            betwvmean[k] = np.mean(bv)
            betwvstd[k] = np.std(bv)

            udistmean[k] = udsum / ndists
            udiststd[k] = ( np.sum(ud2sum) - ((np.sum(udsum)**2)/ndists)) / ndists
            wdistmean[k] = wdsum / ndists
            wdiststd[k] = ( np.sum(wd2sum) - ((np.sum(wdsum)**2)/ndists)) / ndists
            nvertices[k] = len(g.vs)
            nedges[k] = len(g.es)
            ##########################################################

            a = allareas[k][2:] # 0: skeleton, 1: background
            arel = a / np.sum(a)
            entropy = -np.sum(arel * np.log(arel))
            evenness = np.exp(entropy) / len(arel)
            st = [k, nvertices[k], nedges[k], len(a), np.sum(a),
                    np.mean(a), np.std(a), np.std(a)/np.mean(a),
                    np.min(a), np.max(a), entropy, evenness,
                    segmean[k], segstd[k],
                    udistmean[k], udiststd[k], wdistmean[k], wdiststd[k],
                    betwvmean[k], betwvstd[k]
                    ]
            fh.write(','.join([ str(s) for s in st]) + '\n')
        except:
            errors.append(k)
    info('Errors:')
    info(errors)
    fh.close()

##########################################################
def generate_test_graphs(outdir):
    sz = 20 
    g = igraph.Graph.Lattice([sz, sz], circular=False) # lattice
    coords = np.array(g.layout('grid').coords)
    _range = (np.max(coords, 0) - np.min(coords, 0)) * 10
    g.vs['x'] = coords[:, 0] / _range[0]
    g.vs['y'] = coords[:, 1] / _range[1]
    outfilename = pjoin(outdir, 'lattice.graphml')
    igraph.write(g, outfilename, 'graphml')

    g = igraph.Graph.Erdos_Renyi(60, p=1) # lattice
    coords = np.array(g.layout('random').coords)
    _range = (np.max(coords, 0) - np.min(coords, 0)) * 10
    g.vs['x'] = coords[:, 0] / _range[0]
    g.vs['y'] = coords[:, 1] / _range[1]
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
        g = g.components(mode='weak').giant()
        if 'x' in g.vs.attributes():
            g.vs['x'] = np.array(g.vs['x']).astype(float)
            g.vs['y'] = np.array(g.vs['y']).astype(float)
        else:
            g.vs['x'] = np.array(g.vs['posx']).astype(float)
            g.vs['y'] = np.array(g.vs['posy']).astype(float)

        del g.vs['id'] # Avoid future warnings
        for e in g.es:
            try:
                e['weight'] = float(e['length'])
            except:
                coordu = np.array([ g.vs[e.source]['y'], g.vs[e.source]['x'] ])
                coordv = np.array([ g.vs[e.target]['y'], g.vs[e.target]['x'] ])
                e['weight'] = haversine(coordu, coordv, unit='km') # in meters

        igraph.write(g, outpath, 'graphml')

##########################################################
def get_maps_ranges(graphsdir):
    info('Getting map ranges ...')
    ranges = {}
    for filepath in os.listdir(graphsdir):
        k = os.path.splitext(filepath)[0]
        if not filepath.endswith('.graphml'): continue
        g = igraph.Graph.Read(pjoin(graphsdir, filepath))
        lon = g.vs['x']
        lat = g.vs['y']
        ranges[k] = np.array([np.min(lon), np.min(lat), np.max(lon), np.max(lat)])
    return ranges

##########################################################
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('graphsdir', help='Graphs directory')
    parser.add_argument('outdir', default='/tmp/', help='Output directory')
    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s] %(message)s',
    datefmt='%Y%m%d %H:%M', level=logging.INFO)

    skeldir = pjoin(args.outdir, 'skel')
    compdir = pjoin(args.outdir, 'comp')
    weightdir = pjoin(args.outdir, 'weighted')

    # generate_test_graphs(args.graphsdir)

    plot_distributions2('/tmp/')
    return

    add_weights_to_edges(args.graphsdir, weightdir)
    lonlatranges = get_maps_ranges(weightdir)
    plot_graph_raster(weightdir, skeldir)
    components = get_components_from_raster(skeldir, args.outdir)
    generate_components_vis(components, compdir)
    allareas = calculate_block_areas(components, lonlatranges, args.outdir)
    compute_statistics(weightdir, allareas, args.outdir)
    plot_distributions(args.outdir)

if __name__ == "__main__":
    main()

