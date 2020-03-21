#!/usr/bin/env python3
"""Understand the evenness concept
"""

import argparse
import logging
import time
import os
from os.path import join as pjoin
from logging import debug, info
import matplotlib.pyplot as plt

import numpy as np

def normalize_0_1(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def main():
    t0 = time.time()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', default='/tmp/', help='Output directory')
    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s] %(message)s',
    datefmt='%Y%m%d %H:%M', level=logging.INFO)

    if not os.path.isdir(args.outdir): os.mkdir(args.outdir)

    samplesz = 80
    data = {}
    data['exponential'] = np.random.exponential(size=samplesz)
    data['uniform'] = np.ones(samplesz)*.3
    data['uniform'][1] = 0.5

    figsz = 5
    fig, axs = plt.subplots(1, 2, figsize=(2*figsz, 1*figsz))
    epsilon = .01

    for i, k in enumerate(data.keys()):
        d = data[k]
        p = normalize_0_1(d) + epsilon
        # print(p)
        diventropy = -np.sum(p * np.log(p)) + epsilon
        evenness = np.exp(diventropy) / len(p)
        axs[i].hist(d, bins=10, density=True)
        axs[i].set_title('{} evenness: {:.3f}'.format(k, evenness))
    plt.savefig(pjoin(args.outdir, 'evenness.pdf'))
    # plt.show()

    info('Elapsed time:{}'.format(time.time()-t0))
if __name__ == "__main__":
    main()

