#!/usr/bin/env python3

from numpy import loadtxt
from glob import glob1
from pylab import plot, grid, xlabel, ylabel, title, savefig, clf
from os.path import basename, exists, join, splitext
from multiprocessing import Pool

import os
import argparse


def make_argument_parser():
    parser = argparse.ArgumentParser(description="Plot graphics using .dat files in directory")
    parser.add_argument('input_directory', help="input directory")
    parser.add_argument('-j', '--processes', help="number of processes in the pool", type=int)
    return parser


def parse_num(string):
    return int(string[5])


def plot_figure(index, table, basenames, ans_directory_name):
    clf()
    xlabel("x")
    ylabel(basenames[index][0], dict(rotation=0,
                                     horizontalalignment='right',
                                     verticalalignment='center',
                                     x=-0.01))
    title("%s(x, %d)" % (basenames[index][0], parse_num(basenames[index]) + 1))
    plot(table[:,0], table[:,1])
    grid(True)
    filename = join(ans_directory_name, '%s.png' % basenames[index])
    savefig(filename)
    print('Saved %s' % filename)


def main():
    args = make_argument_parser().parse_args()
    results_dir = args.input_directory

    for subdir in os.listdir(results_dir):
        files = glob1(join(results_dir, subdir), '*.dat')
        print('Files in %s: %s' % (subdir, ', '.join(files) or 'no files'))
        ans_directory_name = join(results_dir, subdir, 'png')

        basenames = [splitext(f)[0] for f in files]

        if not exists(ans_directory_name):
            os.makedirs(ans_directory_name)

        with Pool(processes=args.processes) as pool:
            for index, filename in enumerate(files):
                f = join(results_dir, subdir, filename)
                pool.apply_async(plot_figure,
                                 (index, loadtxt(f), basenames, ans_directory_name))
            pool.close()
            pool.join()


if __name__ == '__main__':
    main()
