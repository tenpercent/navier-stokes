#!/usr/bin/env python3

from numpy import loadtxt
from glob import glob
from pylab import plot, grid, xlabel, ylabel, title, savefig, clf
from os.path import basename, exists, join, splitext

import os
import argparse


def make_argument_parser():
    parser = argparse.ArgumentParser(description="Plot graphics using .dat files in directory")
    parser.add_argument('input_directory', help="input directory")
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
    title("%s(x, t) at t = %d" % (basenames[index][0], parse_num(basenames[index]) + 1))
    plot(table[:,0], table[:,1])
    grid(True)
    filename = join(ans_directory_name, '%s.png' % basenames[index])
    savefig(filename)
    print('Saved %s' % filename)


def main():
    args = make_argument_parser().parse_args()
    results_dir = args.input_directory

    results_subdirs = next(os.walk(results_dir))[1]

    for subdir in results_subdirs:
        files = glob(join(results_dir, subdir, '*.dat'))
        print('Files in %s: %s' % (subdir, ', '.join(files) or 'no files'))
        ans_directory_name = join(results_dir, subdir, 'png')

        basenames = [splitext(basename(f))[0] for f in files]

        if not exists(ans_directory_name):
            os.makedirs(ans_directory_name)

        for index, f in enumerate(files):
            plot_figure(index, loadtxt(f), basenames, ans_directory_name)


if __name__ == '__main__':
    main()
