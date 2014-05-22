from numpy import loadtxt
from glob import glob
from pylab import plot, grid, xlabel, ylabel, title, savefig, clf

import os
import argparse
import sys


def make_argument_parser():

    parser = argparse.ArgumentParser(description="Plot graphics using .dat files in directory")
    parser.add_argument('-i', help="input directory")

    if not len(sys.argv) > 1:
        msg = "empty argument list"
        parser.print_help()
        raise argparse.ArgumentTypeError(msg)

    return parser


def parse_num(string):
    return int (string[5])


def plot_figure (index, table, basenames, ans_directory_name):
    clf()
    xlabel("x")
    ylabel(basenames[index][0], dict(rotation=0,
                                     horizontalalignment='right',
                                     verticalalignment='center',
                                     x=-0.01))
    title("%s(x, t) at t = %d" % (basenames[index][0], parse_num(basenames[index]) + 1))
    plot(table[:,0], table[:,1])
    grid(True)
    savefig('%s/%s.png' % (ans_directory_name, basenames[index]))
    return


def main():
    args = make_argument_parser().parse_args()
    results_dir = args.i

    results_subdirs = os.walk(results_dir).next()[1]

    for subdir in results_subdirs:

        files = glob(results_dir + subdir + '/' + "*.dat")
        ans_directory_name = results_dir + subdir + '/png'

        tables = [loadtxt(f) for f in files]
        basenames = [os.path.splitext(os.path.basename(f))[0] for f in files]

        if not os.path.exists(ans_directory_name):
            os.makedirs(ans_directory_name)

        for (index, table) in enumerate(tables):
            plot_figure(index, table, basenames, ans_directory_name)
    return


main()
