from numpy import loadtxt
from glob import glob
from pylab import plot, savefig, clf
from os.path import splitext, basename

def main():
  results_path = "../../build/results/"
  files = glob(results_path + "*.dat")
  tables = [loadtxt(f, skiprows=1) for f in files]
  basenames = [splitext(basename(f))[0] for f in files]

  for (index, table) in enumerate(tables):
    plot(table[:,0], table[:,1])
    savefig('%s.png' % basenames[index])
    clf()
  return

main()
