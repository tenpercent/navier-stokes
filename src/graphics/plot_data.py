from numpy import loadtxt
from glob import glob
from pylab import plot, savefig

def main():
  results_path = "../../build/results/"
  files = glob(results_path + "*.dat")
  tables = [loadtxt(f, skiprows=1) for f in files]
  for (index, table) in enumerate(tables):
    plot(table[:,0], table[:,1])
    savefig('%d.png' % index)
  return

main()
