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
import cv2
import pickle as pkl
import pandas as pd

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
from itertools import groupby

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
def colorize(labels):
    ulabels = np.unique(labels)
    _colors = np.random.randint(0, 255, (len(ulabels), 3))

    labeled_img = np.zeros((labels.shape[0], labels.shape[1], 3), dtype=np.uint8)
    for i, lab in enumerate(ulabels[1:]): # 0: skeleton, 2: external part
        z = np.where(labels == lab)
        labeled_img[z] = _colors[i]

    return labeled_img

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
def compute_areas_entropy(outdir):
    areaentropy = {}
    allareas = pkl.load(open(pjoin(outdir, 'areas.pkl'), 'rb'))
    del allareas['0662602_Rolling_Hills']

    binsize = [0.001, 0.01, 0.1, 1, 10]
    for b in binsize:
        binstr = 'areasentropy' + str(b).replace('.', '')
        areaentropy[binstr] = {}
        for idx, areas in allareas.items():
            areamax = np.max(areas)

            nticks = int(areamax / b) + 1
            _ticks = np.ndarray(nticks)

            _ticks[0] = 0
            for i in range(1, nticks):
                _ticks[i] = _ticks[i-1] + b
            hist,_ = np.histogram(areas, bins=_ticks)

            hist = hist[hist != 0] # Removing values=0
            relfreq = hist / np.sum(hist)

            areaentropy[binstr][idx] = -np.sum(relfreq * np.log(relfreq))
    return areaentropy

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
                e['weight'] = float(e['length']) / 1000
            except:
                coordu = np.array([ g.vs[e.source]['y'], g.vs[e.source]['x'] ])
                coordv = np.array([ g.vs[e.target]['y'], g.vs[e.target]['x'] ])
                e['weight'] = haversine(coordu, coordv) # in meters

        igraph.write(g, outpath, 'graphml')

##########################################################
def get_maps_ranges(graphsdir, outdir):
    outpath = pjoin(outdir, 'lonlatranges.pkl')

    if os.path.exists(outpath):
        info('{} already exists. Skipping ...'.format(outpath))
        return pkl.load(open(outpath, 'rb'))

    info('Getting map ranges ...')
    ranges = {}
    for filepath in os.listdir(graphsdir):
        k = os.path.splitext(filepath)[0]
        if not filepath.endswith('.graphml'): continue
        g = igraph.Graph.Read(pjoin(graphsdir, filepath))
        lon = g.vs['x']
        lat = g.vs['y']
        ranges[k] = np.array([np.min(lon), np.min(lat), np.max(lon), np.max(lat)])

    pkl.dump(ranges, open(outpath, 'wb'))
    return ranges

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
def get_components_from_raster(rasterdir, outdir):
    outpath = pjoin(outdir, 'components.pkl')

    if os.path.exists(outpath):
        info('{} already exists. Skipping ...'.format(outpath))
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
    fh.write('city,nvertices,nedges,degreeentropy,diameter,nblocks,areasum,' \
             'areamean,areastd,areacv,areamin,areamax,areasentropy0001,' \
             'areasentropy001,areasentropy01,areasentropy1,areasentropy10,' \
             'areadiventropy,areaeveness,segmean,segstd,udistmean,udiststd,' \
             'wdistmean,wdiststd,wdistfilteredmean,wdistfilteredstd,betwvmean,betwvstd\n')

    errors = []

    for k in sorted(allareas.keys()):
        info(' *' + k)
        # if True:
        try:
            filepath = pjoin(graphsdir, k + '.graphml')
            g = igraph.Graph.Read(filepath)

            n = len(g.vs)
            ndists = int((n * (n-1)) / 2)

            udsum = 0.0
            ud2sum = 0.0
            wdsum = 0.0
            wd2sum = 0.0
            wdsumfiltered = 0.0
            wd2sumfiltered = 0.0
            ndistsfiltered = 0

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

                wdsumfiltered += np.sum(aux[aux > 1])
                wd2sumfiltered += np.sum(np.square(aux[aux>1]))
                ndistsfiltered += len(aux[aux>1])

            segmean = np.mean(g.es['weight'])
            segstd = np.std(g.es['weight'])
            
            bv = g.betweenness()
            betwvmean = np.mean(bv)
            betwvstd = np.std(bv)

            udistmean = udsum / ndists
            udiststd = ( np.sum(ud2sum) - ((np.sum(udsum)**2)/ndists)) / ndists
            wdistmean = wdsum / ndists
            wdiststd = ( np.sum(wd2sum) - ((np.sum(wdsum)**2)/ndists)) / ndists
            wdistfilteredmean = wdsumfiltered / ndistsfiltered
            wdistfilteredstd = ( np.sum(wd2sumfiltered) - ((np.sum(wdsumfiltered)**2)/ndistsfiltered)) / ndistsfiltered
            nvertices = len(g.vs)
            nedges = len(g.es)
            diameter = g.diameter(weights='weight')
            ##########################################################

            a = allareas[k][2:] # 0: skeleton, 1: background
            arel = a / np.sum(a)
            diventropy = -np.sum(arel * np.log(arel))
            evenness = np.exp(diventropy) / len(arel)

            ##########################################################
            degrees = g.degree()
            N = np.sum(degrees)
            freq = np.array([len(list(group)) for key, group in groupby(degrees)]) / N
            degreeentropy = -np.sum(freq * np.log(freq))
            ##########################################################
            binsize = [0.001, 0.01, 0.1, 1, 10]
            areasentropy = {}
            for b in binsize:
                areas = allareas[k]
                areamax = np.max(areas)
                nticks = int(areamax / b) + 1
                _ticks = np.ndarray(nticks)

                _ticks[0] = 0
                for i in range(1, nticks):
                    _ticks[i] = _ticks[i-1] + b
                hist,_ = np.histogram(areas, bins=_ticks)

                hist = hist[hist != 0] # Removing values=0
                relfreq = hist / np.sum(hist)
                areasentropy[b] = -np.sum(relfreq * np.log(relfreq))

            ##########################################################
            st = [k, nvertices, nedges, degreeentropy, diameter,
                  len(a), np.sum(a), np.mean(a), np.std(a), np.std(a)/np.mean(a),
                  np.min(a), np.max(a),
                  areasentropy[0.001], areasentropy[0.01], areasentropy[0.1],
                  areasentropy[1],areasentropy[10], diventropy, evenness,
                  segmean, segstd, udistmean, udiststd, wdistmean, wdiststd,
                  wdistfilteredmean, wdistfilteredstd, betwvmean, betwvstd,
                  ]
            fh.write(','.join([ str(s) for s in st]) + '\n')
        except:
            errors.append(k)

    info('Errors:')
    info(errors)
    fh.close()

##########################################################
def plot_distributions(outdir, lonlatranges):
    """Plot distributions for each array of areas

    Args:
    allareas(dict): filename as key and areas as values
    epsilon(float): value to consider as invalid area
    """

    info('Plotting distributions of data from {} ...'.format(pjoin(outdir, 'summary.csv')))
    figs = {}
    for k in ['meanpathvsdiventr', 'meanpathvsareasentr001',  'meanpathvsareasentr01',
              'meanpathvsareasentr1', 'meanpathvsdiventrfiltered',
              'degrentrvsareaentr0001', 'degrentrvsareaentr001',
              'degrentrvsareaentr01', 'degrentrvsareaentr1', 'degrentrvsareaentr10',
              'degrentrvsareadiventr',
              'entropydiagonallog', 'meanpathcvvsareasdiventr',
              'meanpathvsdegrentr',  'meanpathvsareaevenness', 'meanpathcvvsareaevenness',
              'entropyhist', 'areaevennesshist']:
        figs[k] = go.Figure()

    ########################################################## Area plots
    areaspath = pjoin(outdir, 'areas.pkl')
    allareas = pkl.load(open(areaspath, 'rb'))
    allareas = filter_areas(allareas) # 0: skeleton, 1: background

    df = pd.read_csv(pjoin(outdir, 'summary.csv'))

    ########################################################## Entropy plots
    for _, row in df.iterrows():
         figs['meanpathvsdiventr'].add_trace(go.Scatter(
            x=[row.areadiventropy],
            y=[np.log(row.wdistmean/row.segmean)],
            mode='markers', marker_size=10, name=row.city,))
    entropylogpearson = pearsonr(df.areadiventropy, np.log(df.wdistmean/df.segmean))
    figs['meanpathvsdiventr'].update_layout(
            title="Entropy of areas vs. logarithm of the normalized average path length (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Entropy of block area",
            yaxis_title="Logarithm of the normalized average path length",
            # yaxis_type="log"
            )

    for _, row in df.iterrows():
         figs['meanpathvsareasentr001'].add_trace(go.Scatter(
            x=[row.areasentropy001],
            y=[row.wdistmean/row.segmean/row.diameter],
            mode='markers', marker_size=10, name=row.city,))
    entropylogpearson = pearsonr(df.areasentropy001, df.wdistmean/df.segmean/df.diameter)
    figs['meanpathvsareasentr001'].update_layout(
            title="Entropy of areas vs. the normalized average path length (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Entropy of block area",
            yaxis_title="Normalized average path length (/meanseglenght/diameter)",
            # yaxis_type="log"
            )

    for _, row in df.iterrows():
         figs['meanpathvsareasentr01'].add_trace(go.Scatter(
            x=[row.areasentropy01],
            y=[row.wdistmean/row.segmean/row.diameter],
            mode='markers', marker_size=10, name=row.city,))
    entropylogpearson = pearsonr(df.areasentropy01, df.wdistmean/df.segmean/df.diameter)
    figs['meanpathvsareasentr01'].update_layout(
            title="Entropy of areas vs. the normalized average path length (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Entropy of block area",
            yaxis_title="Normalized average path length (/meanseglenght/diameter)",
            # yaxis_type="log"
            )

    for _, row in df.iterrows():
         figs['meanpathvsareasentr1'].add_trace(go.Scatter(
            x=[row.areasentropy1],
            y=[row.wdistmean/row.segmean/row.diameter],
            mode='markers', marker_size=10, name=row.city,))
    entropylogpearson = pearsonr(df.areasentropy1, df.wdistmean/df.segmean/df.diameter)
    figs['meanpathvsareasentr1'].update_layout(
            title="Entropy of areas vs. the normalized average path length (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Entropy of block area",
            yaxis_title="Normalized average path length (/meanseglenght/diameter)",
            # yaxis_type="log"
            )

    df2 = df[ (df.areasum > 5) & (df.areasum < 7) ].copy()
    if df2.shape[0] == 0: info('Empty dataframe after filtering')
    df2.city.to_csv(pjoin(outdir, 'cities_130.csv'), header=False, index=False)
    for _, row in df2.iterrows():
         figs['meanpathvsdiventrfiltered'].add_trace(go.Scatter(
            x=[row.areadiventropy],
            y=[row.wdistmean],
            mode='markers',
            marker_size=10,
            name=row.city,
            ))
    entropylogpearson = pearsonr(df2.areadiventropy, np.log(df2.wdistmean))
    figs['meanpathvsdiventrfiltered'].update_layout(
            title="Divisional entropy of areas vs. average path length non-normalized (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Divisional entropy of block area",
            yaxis_title="Average path length",
            yaxis_type="log"
            )

    for _, row in df.iterrows():
         figs['degrentrvsareaentr0001'].add_trace(go.Scatter(
            x=[row['areasentropy0001']],
            y=[row.degreeentropy],
            mode='markers', marker_size=10, name=row.city,))
    entropylogpearson = pearsonr(df['areasentropy0001'], df.degreeentropy)
    figs['degrentrvsareaentr0001'].update_layout(
            title="Entropy of areas (bin=0.001) vs. Entropy of degree (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Entropy of block area",
            yaxis_title="Entropy of degrees",
            )

    for _, row in df.iterrows():
         figs['degrentrvsareaentr001'].add_trace(go.Scatter(
            x=[row['areasentropy001']],
            y=[row.degreeentropy],
            mode='markers', marker_size=10, name=row.city,))
    entropylogpearson = pearsonr(df['areasentropy001'], df.degreeentropy)
    figs['degrentrvsareaentr001'].update_layout(
            title="Entropy of areas (bin=0.01) vs. Entropy of degree (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Entropy of block area",
            yaxis_title="Entropy of degrees",
            )

    for _, row in df.iterrows():
         figs['degrentrvsareaentr01'].add_trace(go.Scatter(
            x=[row['areasentropy01']],
            y=[row.degreeentropy],
            mode='markers', marker_size=10, name=row.city,))
    entropylogpearson = pearsonr(df['areasentropy01'], df.degreeentropy)
    figs['degrentrvsareaentr01'].update_layout(
            title="Entropy of areas (bin=0.1) vs. Entropy of degree (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Entropy of block area",
            yaxis_title="Entropy of degrees",
            )

    for _, row in df.iterrows():
         figs['degrentrvsareaentr1'].add_trace(go.Scatter(
            x=[row['areasentropy1']],
            y=[row.degreeentropy],
            mode='markers', marker_size=10, name=row.city,))
    entropylogpearson = pearsonr(df['areasentropy1'], df.degreeentropy)
    figs['degrentrvsareaentr1'].update_layout(
            title="Entropy of areas (bin=1) vs. Entropy of degree (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Entropy of block area",
            yaxis_title="Entropy of degrees",
            )


    for _, row in df.iterrows():
         figs['degrentrvsareaentr10'].add_trace(go.Scatter(
            x=[row['areasentropy10']],
            y=[row.degreeentropy],
            mode='markers', marker_size=10, name=row.city,))
    entropylogpearson = pearsonr(df['areasentropy10'], df.degreeentropy)
    figs['degrentrvsareaentr10'].update_layout(
            title="Entropy of areas (bin=10) vs. Entropy of degree (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Entropy of block area",
            yaxis_title="Entropy of degrees",
            )

    for _, row in df.iterrows():
         figs['degrentrvsareadiventr'].add_trace(go.Scatter(
            x=[row.areadiventropy],
            y=[row.degreeentropy],
            mode='markers', marker_size=10, name=row.city,))
    entropylogpearson = pearsonr(df.areadiventropy, df.degreeentropy)
    figs['degrentrvsareadiventr'].update_layout(
            title="Divisional entropy of areas vs. Entropy of degree (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Divisional entropy of block area",
            yaxis_title="Entropy of degrees",
            )

    for _, row in df.iterrows():
        figs['meanpathcvvsareasdiventr'].add_trace(go.Scatter(
            x=[row.areadiventropy],
            y=[np.log(row.wdiststd/row.wdistmean/row.segmean)],
            mode='markers', marker_size=10, name=row.city,))
    entropycvlogpearson = pearsonr(df.areadiventropy, np.log(df.wdiststd/df.wdistmean/df.segmean))
    figs['meanpathcvvsareasdiventr'].update_layout(
            title="Entropy of areas vs. logarithm of the coefficient of variation of the normalized shortest path length (Pearson: {:.2f})".\
                    format(entropycvlogpearson[0]),
            xaxis_title="Entropy of block area",
            yaxis_title="Logarithm of the coefficient of variation of the normalized shortest path length",
            # yaxis_type='log',
            )

    for _, row in df.iterrows():
         figs['meanpathvsdegrentr'].add_trace(go.Scatter(
            x=[row.degreeentropy],
            y=[row.wdistmean/row.segmean],
            mode='markers', marker_size=10, name=row.city,))
    entropylogpearson = pearsonr(df.degreeentropy, np.log(df.wdistmean/df.segmean))
    figs['meanpathvsdegrentr'].update_layout(
            title="Entropy of degrees vs. mean path length norm. by mean segment length (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Entropy of degrees",
            yaxis_title="Mean path length (normalized by mean segment length)",
            yaxis_type="log"
            )

    
    for _, row in df.iterrows():
         figs['meanpathvsareaevenness'].add_trace(go.Scatter(
            x=[row.areaeveness],
            y=[np.log(row.wdistmean)],
            mode='markers', marker_size=10, name=row.city,))
    entropylogpearson = pearsonr(df.areaeveness, np.log(df.wdistmean))
    figs['meanpathvsareaevenness'].update_layout(
            title="Logarithm of the average path length vs. evenness of block areas (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Evenness of block areas",
            yaxis_title="Logarithm of the average path length",
            # yaxis_type="log"
            )

    # print(np.log(df.wdiststd/df.wdistmean))
    # for _, row in df.iterrows():
         # figs['meanpathcvvsareaevenness'].add_trace(go.Scatter(
            # x=[row.areaeveness],
            # y=[np.log(df.wdiststd/df.wdistmean)],
            # mode='markers', marker_size=10, name=row.city,))

    # figs['meanpathcvvsareaevenness'].add_trace(go.Scatter(
        # x=[row.areaeveness],
        # y=[np.log(df.wdiststd/df.wdistmean)],
        # mode='markers', marker_size=10, name=row.city,))

    entropylogpearson = pearsonr(df.areaeveness, np.log(df.wdiststd/df.wdistmean))

    figs['meanpathcvvsareaevenness'].add_trace(go.Scatter(
        x=[df.areaeveness],
        y=[np.log(df.wdiststd/df.wdistmean)],
        mode='markers', marker_size=10,
        line=go.scatter.Line(color="gray"),
    ))

    figs['meanpathcvvsareaevenness'].update_layout(
            title="Logarithm of the coefficient of variation of the average path length vs. evenness of block areas (Pearson: {:.2f})".\
                    format(entropylogpearson[0]),
            xaxis_title="Evenness of block areas",
            yaxis_title="Logarithm of the coefficient of variation of the average path length",
            )

    print('foo')

    df2 = df.copy()
    import plotly.express as px
    df2['x'] = df.areaeveness
    df2['y'] = np.log(df.wdiststd/df.wdistmean)
    fig = px.scatter(df2, x='x', y='y', trendline="ols")
    # fig.show()
    plotly.offline.plot(fig,
                        filename='/tmp/foo.html',
                        auto_open=False)

    figs['entropyhist'].add_trace(
            go.Histogram(x=df.areadiventropy,
                histnorm='probability',
                ))
    figs['entropyhist'].update_layout(
            title="Entropy of block area",
            xaxis_title="Entropy",
            yaxis_title="Relative frequency",
            )

    figs['areaevennesshist'].add_trace(
            go.Histogram(x=df.areaeveness,
                histnorm='probability',
                ))
    figs['areaevennesshist'].update_layout(
            title="Evenness of block areas",
            xaxis_title="Evenness of block areas",
            yaxis_title="Relative frequency",
            )
    ########################################################## Save to file
    for k in figs.keys():
        plotly.offline.plot(figs[k],
                            filename=pjoin(outdir, k + '.html'),
                            auto_open=False)

##########################################################
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('graphsdir', help='Graphs directory')
    parser.add_argument('outdir', default='/tmp/', help='Output directory')
    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s] %(message)s',
    datefmt='%Y%m%d %H:%M', level=logging.INFO)

    skeldir = pjoin(args.outdir, 'skel/')
    compdir = pjoin(args.outdir, 'comp/')
    weightdir = pjoin(args.outdir, 'weighted/')

    # generate_test_graphs(args.graphsdir)

    # add_weights_to_edges(args.graphsdir, weightdir)
    lonlatranges = get_maps_ranges(weightdir, args.outdir)
    # plot_graph_raster(weightdir, skeldir)
    # components = get_components_from_raster(skeldir, args.outdir)
    # generate_components_vis(components, compdir)
    # allareas = calculate_block_areas(components, lonlatranges, args.outdir)
    # compute_statistics(weightdir, allareas, args.outdir)
    plot_distributions(args.outdir, lonlatranges)

    # areasentropy = compute_areas_entropy(args.outdir)
    # print(list(areasentropy))
    # z = pd.DataFrame.from_dict()

    # diams = {}
    # summary = pd.read_csv('/tmp/summary.csv')
    # i = 0
    # for c in summary.city:
        # # if 'Anaheim' in c: break
        # print(i, c)
        # i += 1
        # graphpath = pjoin('/tmp/20191108-cadataverse/weighted/', c+'.graphml')
        # diams[c] = igraph.Graph.Read(graphpath).diameter(weights='weight')

    # z = pd.DataFrame.from_dict(diams, orient='index')
    # z.index.names = ['city']
    # summary.set_index('city', drop=True, inplace=True)
    # zz = pd.concat([summary, z], axis = 1, sort=True)
    # zz.to_csv('/tmp/foo.csv', index_label='city')


if __name__ == "__main__":
    main()

