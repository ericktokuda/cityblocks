#!/usr/bin/env python3
"""Plot blocks from toy experiments
"""

import argparse
import logging
from os.path import join as pjoin
from logging import debug, info
import colorsys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random

def generate_colors(n):
    """Generate colors

    Args:
    n(int): number of colors to be generated

    Returns:
    list: set of colors
    """

    HSV_tuples = [(x*1.0/n, 0.5, 0.5) for x in range(n)]
    RGB_tuples = list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))
    random.shuffle(RGB_tuples)
    return np.array(RGB_tuples)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', default='/tmp/', help='Output directory')
    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s] %(message)s',
    datefmt='%Y%m%d %H:%M', level=logging.DEBUG)

    dsvpath = '/tmp/blocks.dsv'
    fh = open(dsvpath)

    fblockcols = 30
    plotrate = 10
    colors = generate_colors(600)

    t = 0
    for line in fh:
        blocks = line.strip().split(';')[:-1]
        # print('t:{}, nblocks:{}'.format(t, len(blocks)))

        for bl in blocks:
            blid, fblockstr = bl.split(':')
            blid = int(blid)
            fblockids = [ int(xx) for xx in fblockstr.split(',')[:-1]]
            fblockcoords = [ [fblid//fblockcols, fblid%fblockcols] for fblid in fblockids] 
            fblockcoords = np.array(fblockcoords)

            if t % plotrate == 0:
                plt.scatter(fblockcoords[:, 1], fblockcoords[:, 0],
                            c=colors[blid], marker='s', s=60)


        if t % plotrate == 0:
            plt.yticks([])
            plt.savefig('/tmp/{:03d}.png'.format(t), dpi=200)
            plt.clf()
        t += 1

    fh.close()

if __name__ == "__main__":
    main()

