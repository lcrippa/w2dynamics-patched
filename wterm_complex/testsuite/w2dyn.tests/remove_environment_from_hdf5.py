#!/usr/bin/env python3
import tempfile
from sys import argv
import h5py

hf = h5py.File(argv[1], "r+")
del hf['/.environment']
hf.close()

print("Group /.environment removed. Use h5repack to actually remove associated data.")
