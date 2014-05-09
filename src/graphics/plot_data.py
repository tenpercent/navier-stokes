from numpy import loadtxt
from glob import glob
from pylab import plot, savefig, clf
from os.path import splitext, basename, exists
from os import makedirs

def main():
  results_path = "../../build/results/"
  ans_directory_name = "png"
  files = glob(results_path + "*.dat")
  tables = [loadtxt(f, skiprows=1) for f in files]
  basenames = [splitext(basename(f))[0] for f in files]

  if not exists(ans_directory_name):
    makedirs(ans_directory_name)

  for (index, table) in enumerate(tables):
    plot(table[:,0], table[:,1])
    savefig('%s/%s.png' % (ans_directory_name, basenames[index]))
    clf()
  return

main()
