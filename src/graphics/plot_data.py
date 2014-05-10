from numpy import loadtxt
from glob import glob
from pylab import plot, grid, xlabel, ylabel, title, savefig, clf
from os.path import splitext, basename, exists
from os import makedirs

ans_directory_name = "png"
results_path = "../../build/results/"


def plot_figure (index, table, basenames):
    clf()
    xlabel("x")
    ylabel("values")
    title(basenames[index])
    plot(table[:,0], table[:,1])
    grid(True)
    savefig('%s/%s.png' % (ans_directory_name, basenames[index]))
    return


def main():
    files = glob(results_path + "*.dat")

    tables = [loadtxt(f) for f in files]
    basenames = [splitext(basename(f))[0] for f in files]

    if not exists(ans_directory_name):
        makedirs(ans_directory_name)

    for (index, table) in enumerate(tables):
        plot_figure(index, table, basenames)
    return


main()
