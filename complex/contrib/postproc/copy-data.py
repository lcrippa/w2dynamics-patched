#!/usr/bin/env python
import numpy as np
import h5py as hdf5
import argparse
import sys

__version__ = "0.1"


def parser():
    p = argparse.ArgumentParser(
            description="Copy elements from one HDF5 file to another")
    p.add_argument('infile', type=str)
    p.add_argument('outfile', type=str)
    p.add_argument('path', type=str, nargs='+')
    return p


def main(args):
    args = parser().parse_args(args)
    with hdf5.File(args.infile, 'r') as infile, \
            hdf5.File(args.outfile, 'r+') as outfile:
        for path in args.path:
            sys.stderr.write(f"{path}\n")
            outfile[path]['value'][...] = infile[path]['value'][...]
            outfile[path]['error'][...] = infile[path]['error'][...]


if __name__ == '__main__':
    main(sys.argv[1:])
