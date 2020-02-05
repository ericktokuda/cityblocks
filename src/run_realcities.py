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
def colorize_random(labels):
    ulabels = np.unique(labels)
    _colors = np.random.randint(0, 255, (len(ulabels), 3))

    labeled_img = np.zeros((labels.shape[0], labels.shape[1], 3), dtype=np.uint8)
    for i, lab in enumerate(ulabels[1:]): # 0: skeleton, 2: external part
        z = np.where(labels == lab)
        labeled_img[z] = _colors[i]

    return labeled_img

##########################################################
def colorize_by_size(labels):
    ulabels, counts = np.unique(labels, return_counts=True)
    m = (np.max(counts))

    labeled_img = np.zeros((labels.shape[0], labels.shape[1], 3), dtype=np.uint8)

    # Label regions
    for i, lab in enumerate(ulabels): # 0: skeleton, 1: external part
        if i < 2: continue
        z = np.where(labels == lab)
        labeled_img[z] = [0, int(((counts[i])/m) * 240) + 15, 0]

    # Label skeleton
    z = np.where(labels == 0)
    labeled_img[z] = [255]*3
    return labeled_img

##########################################################
def compute_raster_real_conversion(rasterranges, lonlatranges):
    conversionfactors = {}
    for k, c in lonlatranges.items():
        if k == '0662602_Rolling_Hills': continue 
        realbboxcoords = np.array([[c[0], c[1]], [c[0], c[3]],
                                   [c[2], c[3]], [c[2], c[1]], ])
        real = calculate_real_area(realbboxcoords)
        r = rasterranges[k]
        raster = (r[2] - r[0]) * (r[3] - r[1])
        conversionfactors[k] = real / raster
    return conversionfactors

##########################################################
def calculate_raster_areas(labels, outdir): # It contains both the skeleton and outer
    info('Computing block areas from components ...')

    allareas = {}
    skelranges = {}
    for k, label in labels.items():
        info(' *' + k)
        if k == '0662602_Rolling_Hills': continue #TODO: fix this
        _, areas = np.unique(label, return_counts=True)
        allareas[k] = np.array(areas) # Eliminate skeleton and external part later
        skelids = np.where(label == 0) # skeleton (cv2)
        skelranges[k] = np.array([ np.min(skelids[0]),  np.min(skelids[1]),
            np.max(skelids[0]), np.max(skelids[1]) ])
    return allareas, skelranges

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

#############################################################
def get_ref_params(graphsdir):
    for filepath in os.listdir(graphsdir):
        if not filepath.endswith('.graphml'): continue

        g = igraph.Graph.Read(pjoin(graphsdir, filepath))

        if 'x' in g.vertex_attributes():
            xs = np.array(g.vs['x']).astype(float)
        else:
            xs = np.array(g.vs['posx']).astype(float)

        if 'y' in g.vertex_attributes():
            ys = np.array(g.vs['y']).astype(float)
        else:
            ys = np.array(g.vs['posy']).astype(float)

        mins = np.array([np.min(xs), np.min(ys)])
        maxs = np.array([np.max(xs), np.max(ys)])
        return g.vcount(), g.ecount(), mins, maxs, filepath

#############################################################
def get_4connected_neighbours_2d(i, j, nx, ny, thoroidal=False):
    """Get 4-connected neighbours. It does not check if there are repeated entries (2x2 or 1x1)
    Args:
    i(int): row of the matrix
    j(int): column of the matrix
    n(int): side of the square matrix
    Returns:
    ndarray 4x2: 4 neighbours indices
    """
    inds = []
    if j > 0: # left
        inds.append([i, j-1])
    elif thoroidal:
        inds.append([i, ny-1])

    if j < ny-1: # right
        inds.append([i, j+1])
    elif thoroidal:
        inds.append([i, 0])

    if i > 0: # top
        inds.append([i-1, j])
    elif thoroidal:
        inds.append([nx-1, j])

    if i < nx-1: # bottom
        inds.append([i+1, j])
    elif thoroidal:
        inds.append([0, j])

    return np.array(inds)

##########################################################
def generate_lattice(ncols, nrows, thoroidal=False, s=10):
    """Generate 2d lattice of side n
    Args:
    n(int): side of the lattice
    thoroidal(bool): thoroidal lattice
    s(float): edge size
    Returns:
    ndarray nx2, ndarray nxn: positions and adjacency matrix (triangular)
    """
    n2 = nrows*ncols
    pos = np.ndarray((n2, 2), dtype=float)
    adj = np.zeros((n2, n2), dtype=int)

    k = 0
    for i in range(nrows):
        for j in range(ncols): # Set positions
            pos[k] = [j*s, i*s]
            k += 1

    for i in range(nrows): # Set connectivity
        for j in range(ncols):
            neighs2d = get_4connected_neighbours_2d(i, j, nrows,
                                                    ncols, thoroidal)
            neighids = np.ravel_multi_index((neighs2d[:, 0],
                                             neighs2d[:, 1]),
                                            (nrows, ncols))
            curidx = np.ravel_multi_index((i, j), (nrows, ncols))

            for neigh in neighids:
                adj[curidx, neigh] = 1
                adj[neigh, curidx] = 1
    return pos, adj

##########################################################
def generate_test_graphs(graphsdir, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    nvertices, nedges, mins, maxs, refcity = get_ref_params(graphsdir)
    range_ = (maxs - mins)

    yxratio = range_[1]  / range_[0]
    info('Generating toy graphs with {} vertices (ref:{})'.format(
        nvertices, refcity))

    gridnx = np.around(np.sqrt(nvertices/yxratio)).astype(int)
    gridny = np.around(gridnx*yxratio).astype(int)
    pos, adj = generate_lattice(gridny, gridnx)
    g = igraph.Graph.Adjacency(adj.tolist(), mode=igraph.ADJ_UNDIRECTED)
    coords = pos
    normalized = (coords - np.min(coords, 0)) / \
        (np.max(coords, 0) - np.min(coords, 0))
    scaled = normalized * (maxs - mins) + mins
    g.vs['x'] = scaled[:, 0]
    g.vs['y'] = scaled[:, 1]

    aux = [ (s[0], s[1]) for s in scaled ]

    g = add_weights_to_edges(g)
    outfilename = pjoin(outdir, 'lattice.pkl')
    pkl.dump(g, open(outfilename, 'wb'))

    delratio = .25
    nedges = g.ecount()
    aux = np.random.rand(int(nedges*delratio))*nedges
    g.delete_edges(aux.astype(int))
    g = g.components(mode='weak').giant()
    g = add_weights_to_edges(g)
    outfilename = pjoin(outdir, 'lattice_{}.pkl'.format(int(delratio*100)))
    pkl.dump(g, open(outfilename, 'wb'))

##########################################################
def add_weights_to_edges(g):
    for e in g.es:
        try: # Assuming we are using harvard dataset
            e['weight'] = float(e['length']) / 1000 # we want in km
        except:
            coordu = np.array([ g.vs[e.source]['y'], g.vs[e.source]['x'] ])
            coordv = np.array([ g.vs[e.target]['y'], g.vs[e.target]['x'] ])
            e['weight'] = haversine(coordu, coordv, unit='km')
    return g

#############################################################
def parse_graphml(graphsdir, weightdir):
    if not os.path.exists(weightdir):
        # info('{} already exists. Skipping ...'.format(weightdir))
        os.makedirs(weightdir)

    for filepath in os.listdir(graphsdir):
        outpath = pjoin(weightdir, filepath.replace('.graphml', '.pkl'))
        if not filepath.endswith('.graphml') or os.path.exists(outpath): continue

        info(' *' + filepath)
        g = igraph.Graph.Read(pjoin(graphsdir, filepath))
        g = g.components(mode='weak').giant()

        if 'x' in g.vs.attributes():
            g.vs['x'] = np.array(g.vs['x']).astype(float)
            g.vs['y'] = np.array(g.vs['y']).astype(float)
        else:
            g.vs['x'] = np.array(g.vs['posx']).astype(float)
            g.vs['y'] = np.array(g.vs['posy']).astype(float)

        if 'id' in g.vertex_attributes():
            del g.vs['id'] # Avoid future warnings

        g = add_weights_to_edges(g)

        pkl.dump(g, open(outpath, 'wb'))

##########################################################
def get_maps_ranges(graphsdir, outdir):
    outpath = pjoin(outdir, 'lonlatranges.pkl')

    if os.path.exists(outpath):
        info('{} already exists. Skipping ...'.format(outpath))
        return

    info('Getting map ranges ...')
    ranges = {}
    for filepath in os.listdir(graphsdir):
        k = os.path.splitext(filepath)[0]
        if not filepath.endswith('.pkl'): continue
        # g = igraph.Graph.Read(pjoin(graphsdir, filepath))
        g = pkl.load(open(pjoin(graphsdir, filepath), 'rb'))
        lon = g.vs['x']
        lat = g.vs['y']
        ranges[k] = np.array([np.min(lon), np.min(lat),
                              np.max(lon), np.max(lat)])

    pkl.dump(ranges, open(outpath, 'wb'))

#########################################################
def plot_graph_raster(graphsdir, skeldir, figscale=20000):
    if not os.path.exists(skeldir):
        # info('{} already exists. Skipping ...'.format(skeldir))
        os.makedirs(skeldir)

    info('Reading graphs from {} ...'.format(graphsdir))
    allareas = {}
    for filepath in os.listdir(graphsdir):
        info(' *' + filepath)
        outpath = pjoin(skeldir, os.path.splitext(filepath)[0] + '.png')
        if not filepath.endswith('.pkl') or os.path.exists(outpath): continue
        g = pkl.load(open(pjoin(graphsdir, filepath), 'rb'))
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

        lonfigsize = (lonrange[1] - lonrange[0])*figscale
        latfigsize = (latrange[1] - latrange[0])*figscale

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
    info('Generating components visualization ...')
    if not os.path.exists(compdir):
        os.makedirs(compdir)

    for k, labels in components.items():
        info(' *' + k)
        outpath = pjoin(compdir, k + '_labels.png')
        if k == '0662602_Rolling_Hills': continue #TODO: fix this
        if os.path.exists(outpath): continue
        # labeled_img = colorize_random(labels)
        labeled_img = colorize_by_size(labels)
        cv2.imwrite(outpath, labeled_img)

##########################################################
def calculate_block_areas(labels, outdir):
    outpath = pjoin(outdir, 'areas.pkl')
    lonlatrangespath = pjoin(outdir, 'lonlatranges.pkl')
    lonlatranges =  pkl.load(open(lonlatrangespath, 'rb'))

    if os.path.exists(outpath):
        info('{} already exists. Skipping ...'.format(outpath))
        return pkl.load(open(outpath, 'rb'))

    areas, skelranges = calculate_raster_areas(labels, outdir)
    conversionfactors = compute_raster_real_conversion(skelranges, lonlatranges)

    for k in areas.keys():
        areas[k] = areas[k] * conversionfactors[k]
    pkl.dump(areas, open(outpath, 'wb'))
    return areas

##########################################################
def compute_statistics(graphsdir, blockareas, blockminarea, outdir):
    """Compute statistics from areas

    Args:
    areas(dict): city as keys and list of areas as values
    """
    outpath = pjoin(outdir, 'results.csv')
    if os.path.exists(outpath):
        info('{} already exists. Skipping ...'.format(outpath))
        return

    info('Computing graph statistics...')
    errors = []

    header = 'city,nvertices,nedges,degreeentropy,diameter,' \
        'nblocksall,nblocksvalid,areasum,' \
        'areamean,areastd,areacv,areamin,areamax,areasentropy0001,' \
        'areasentropy001,areasentropy01,areasentropy1,areasentropy10,' \
        'areadiventropy,areaeveness,segmean,segstd,udistmean,udiststd,' \
        'wdistmean,wdiststd,betwvmean,betwvstd'

    df = pd.DataFrame(columns=header.split(','),
                          index=blockareas.keys())

    for k in sorted(blockareas.keys()):
        info(' *' + k)
        # if True:
        try:
            filepath = pjoin(graphsdir, k + '.pkl')
            # g = igraph.Graph.Read(filepath)
            g = pkl.load(open(filepath, 'rb'))

            nvertices = g.vcount()
            nedges = g.ecount()
            ndists = int((nvertices * (nvertices-1)) / 2)

            udsum = 0.0
            ud2sum = 0.0
            wdsum = 0.0
            wd2sum = 0.0

            ########################################################## avg paths
            for i in range(nvertices): # Calling with all vertice at once crashes
                # Unweighted paths
                aux = np.array(g.shortest_paths(source=i, mode='ALL'))[0][i+1:]
                udsum += np.sum(aux)
                ud2sum += np.sum(np.square(aux))

                # Weighted paths
                aux = np.array(g.shortest_paths(source=i, mode='ALL',
                    weights=g.es['weight']))[0][i+1:]
                wdsum += np.sum(aux)
                wd2sum += np.sum(np.square(aux))

            udistmean = udsum / ndists
            udiststd = ( np.sum(ud2sum) - ((np.sum(udsum)**2)/ndists)) / ndists
            wdistmean = wdsum / ndists
            wdiststd = ( np.sum(wd2sum) - ((np.sum(wdsum)**2)/ndists)) / ndists

            ########################################################## segments
            diameter = g.diameter(weights='weight')
            segmean = np.mean(g.es['weight'])
            segstd = np.std(g.es['weight'])
            
            ########################################################## betweenness
            bv = g.betweenness()
            betwvmean = np.mean(bv)
            betwvstd = np.std(bv)

            ########################################################## node degrees
            degrees = g.degree()
            N = np.sum(degrees)
            freq = np.array([len(list(group)) for key, group in groupby(degrees)]) / N
            degreeentropy = -np.sum(freq * np.log(freq))

            ########################################################## block areas
            areas = blockareas[k][2:] # 0: skeleton, 1: background
            nblocksall = len(areas)
            validind = areas >= blockminarea
            nblocksvalid = np.sum(validind)
            areas = areas[validind]
            arel = areas / np.sum(areas)
            diventropy = -np.sum(arel * np.log(arel))
            evenness = np.exp(diventropy) / len(arel)

            ##########################################################
            binsize = [0.001, 0.01, 0.1, 1, 10]
            areasentropy = {}
            for b in binsize:
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

            ########################################################## exporting
            a = areas
            df.loc[k] = [k, nvertices, nedges, degreeentropy, diameter,
                         nblocksall, nblocksvalid,
                         np.sum(a), np.mean(a), np.std(a), np.std(a)/np.mean(a),
                         np.min(a), np.max(a),
                         areasentropy[0.001], areasentropy[0.01], areasentropy[0.1],
                         areasentropy[1], areasentropy[10], diventropy, evenness,
                         segmean, segstd, udistmean, udiststd, wdistmean, wdiststd,
                         betwvmean, betwvstd,
                         ]
        except:
            errors.append(k)

    df.to_csv(outpath, index=False)

    info('ATTRIB\t\tMEAN(STD)')
    for col in df.columns:
        if col == 'city': continue
        info('{}:\t{:.3f} ({:.3f})'.format(col, np.mean(df[col]), np.std(df[col])))
    info('Errors:')
    info(errors)

##########################################################
def plot_distributions(outdir):
    """Plot distributions for each array of areas

    Args:
    allareas(dict): filename as key and areas as values
    epsilon(float): value to consider as invalid area
    """

    info('Plotting distributions of data from {} ...'.format(pjoin(outdir, 'results.csv')))
    figs = {}
    for k in ['avgpathlengthvsblocksdiventr', 'avgpathlengthvsblocksevenness',
              'avgpathlengthvsdegreestd', 'avgpathlengthvsdegreesnonnullstd',
              'avgpathlengthvsdegreeentr','degreeentrvsblocksdiventr',
              ]:
        figs[k] = go.Figure()

    df = pd.read_csv(pjoin(outdir, 'results.csv'))

    ########################################################## Entropy plots
    for i, row in df.iterrows():
         figs['avgpathlengthvsblocksdiventr'].add_trace(go.Scatter(
            x=[row.avgpathlength],
            y=[row.blocksdiventr],
            mode='markers', marker_size=10, name=str(i),))
    # entropylogpearson = pearsonr(df.areadiventropy, np.log(df.wdistmean/df.segmean))
    figs['avgpathlengthvsblocksdiventr'].update_layout(
            title="Average path length vs. block areas divisional entropy",
            yaxis_title="Average path length",
            xaxis_title="Divisional entropy of block areas",
            # yaxis_type="log"
            )

    ########################################################## Entropy plots
    for i, row in df.iterrows():
         figs['avgpathlengthvsblocksevenness'].add_trace(go.Scatter(
            x=[row.avgpathlength],
            y=[row.blocksevenness],
            mode='markers', marker_size=10, name=str(i),))
    # entropylogpearson = pearsonr(df.areadiventropy, np.log(df.wdistmean/df.segmean))
    figs['avgpathlengthvsblocksevenness'].update_layout(
            title="Average path length vs. block areas evenness",
            yaxis_title="Average path length",
            xaxis_title="Evenness of block areas",
            # yaxis_type="log"
            )
    ########################################################## Entropy plots
    for i, row in df.iterrows():
         figs['avgpathlengthvsdegreestd'].add_trace(go.Scatter(
            x=[row.avgpathlength],
            y=[row.degreestd],
            mode='markers', marker_size=10, name=str(i),))
    # entropylogpearson = pearsonr(df.areadiventropy, np.log(df.wdistmean/df.segmean))
    figs['avgpathlengthvsdegreestd'].update_layout(
            title="Average path length vs. degrees standard deviation (all)",
            yaxis_title="Average path length",
            xaxis_title="Degrees standard deviation",
            # yaxis_type="log"
            )

    ########################################################## Entropy plots
    for i, row in df.iterrows():
         figs['avgpathlengthvsdegreesnonnullstd'].add_trace(go.Scatter(
            x=[row.avgpathlength],
            y=[row.degreesnonnullstd],
            mode='markers', marker_size=10, name=str(i),))
    # entropylogpearson = pearsonr(df.areadiventropy, np.log(df.wdistmean/df.segmean))
    figs['avgpathlengthvsdegreesnonnullstd'].update_layout(
            title="Average path length vs. degrees standard deviation (positive degrees)",
            yaxis_title="Average path length",
            xaxis_title="Degrees standard deviation",
            # yaxis_type="log"
            )

    ########################################################## Entropy plots
    for i, row in df.iterrows():
         figs['avgpathlengthvsdegreeentr'].add_trace(go.Scatter(
            x=[row.avgpathlength],
            y=[row.degreesentr],
            mode='markers', marker_size=10, name=str(i),))
    # entropylogpearson = pearsonr(df.areadiventropy, np.log(df.wdistmean/df.segmean))
    figs['avgpathlengthvsdegreeentr'].update_layout(
            title="Average path length vs. degrees entropy (positive degrees)",
            yaxis_title="Average path length",
            xaxis_title="Entropy of degrees",
            # yaxis_type="log"
            )

    ########################################################## Entropy plots
    for i, row in df.iterrows():
         figs['degreeentrvsblocksdiventr'].add_trace(go.Scatter(
            x=[row.degreesentr],
            y=[row.blocksdiventr],
            mode='markers', marker_size=10, name=str(i),))
    # entropylogpearson = pearsonr(df.areadiventropy, np.log(df.wdistmean/df.segmean))
    figs['degreeentrvsblocksdiventr'].update_layout(
            title="Degrees entropy vs. block areas divisional entropy",
            yaxis_title="Entropy of degrees",
            xaxis_title="Block areas divisional entropy",
            # yaxis_type="log"
            )
    ########################################################## Save to file
    for k in figs.keys():
        plotly.offline.plot(figs[k],
                            filename=pjoin(outdir, k + '.html'),
                            auto_open=False)

##########################################################
def plot_areas_distrib(areaspath, outdir):
    import matplotlib.pyplot as plt
    areas = pkl.load(open(areaspath, 'rb'))
    for k, v in areas.items():
        a = v[2:] # Ignore skeleton and outer part
        bins = np.arange(np.min(a), np.max(a)+0.1, 0.01)
        plt.clf()
        plt.hist(np.around(a, decimals=3), bins=bins)
        plt.title('Histogram for {}'.format(k))
        plt.xlabel('Block area')
        plt.savefig(pjoin(outdir, k + '.png'))

##########################################################
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('graphsdir', help='Graphs directory')
    parser.add_argument('--outdir', required=False, default='/tmp/', help='Output directory')
    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s] %(message)s',
    datefmt='%Y%m%d %H:%M', level=logging.INFO)

    skeldir = pjoin(args.outdir, 'skel/')
    compdir = pjoin(args.outdir, 'comp/')
    weightdir = pjoin(args.outdir, 'weighted/')

    np.random.seed(0)
    figscale = 20000
    blockminarea = 0.0004

    # plot_areas_distrib('/tmp/20200203-CAdataverse_sample/areas.pkl', '/tmp/')
    parse_graphml(args.graphsdir, weightdir)
    generate_test_graphs(args.graphsdir, weightdir)
    get_maps_ranges(weightdir, args.outdir)
    plot_graph_raster(weightdir, skeldir, figscale)
    components = get_components_from_raster(skeldir, args.outdir)
    generate_components_vis(components, compdir)
    areas = calculate_block_areas(components, args.outdir)
    compute_statistics(weightdir, areas, blockminarea, args.outdir)
    return
    # plot_distributions(args.outdir)

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

