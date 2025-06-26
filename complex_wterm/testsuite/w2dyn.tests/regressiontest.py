import os
import subprocess
import argparse
import h5py as hdf5
import numpy as np


def match(x, y):
    """Return indices `i` such that `y[i]` is closest to `x`."""
    x = np.asarray(x)
    y = np.asarray(y)
    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("arguments must be vectors")

    diff = np.abs(x[:, None] - y[None, :])
    return diff.argmin(1)


parser = argparse.ArgumentParser()
parser.add_argument("project_source_dir", help="directory of the DMFT.py and cthyb scripts")
parser.add_argument("test_source_dir", help="directory of the input files and reference data for the regression tests")
parser.add_argument("executable", help="executable: cthyb or DMFT.py")
args = parser.parse_args()

print('args.project_source_dir ', os.path.realpath(args.project_source_dir))
print('args.test_source_dir ', os.path.realpath(args.test_source_dir))
print('args.executable ', args.executable)
print('os.environ["PYTHON_EXECUTABLE"] ', os.environ["PYTHON_EXECUTABLE"] )

os.chdir(args.test_source_dir)

executable = args.project_source_dir + "/" + args.executable
para_in = "Parameters.in"
#print('executable ', executable )
#print('para_in', para_in)

reffilename= "reference_data.hdf5"
testfilename = args.test_source_dir + ".hdf5"
print('testfilename', testfilename)

try:
    os.remove(testfilename)
except:
    print('no old testfile existing!')

arch = subprocess.run([os.environ["PYTHON_EXECUTABLE"],
                       executable,
                       para_in], check=True)

reffile = hdf5.File(reffilename, "r")
testfile = hdf5.File(testfilename, "r")

def compare(obj):
    return np.allclose(reffile[obj][()], testfile[obj][()], rtol=1e-9, atol=1e-9)

with open("axes.dat", "r") as axes:
    for ax in axes:
        if not compare(ax.strip()):
            raise ValueError("Axis " + ax.strip() + " does not match!")

with open("quantities.dat", "r") as observables:
    mismatch = []
    print(f"{'':=<45s} {'':=<10s} {'':=<10s}")
    print(f"{'OBSERVABLE':45s} {'ABS. DEV.':>10s} {'REL. DEV.':>10s}")
    for obs in observables:
        obs = obs.strip()
        cmp_eigenvalues = 'densitymatrix' in obs.split('/')

        ref_data = reffile[obs][...]
        test_data = testfile[obs][...]

        # In this case, only the eigenvalues are meaningful, because the
        # eigenvectors depend on the exact details of the diagonalization
        # and thus on the LAPACK implementation
        if cmp_eigenvalues:
            ref_data = np.linalg.eigvals(ref_data)
            test_data = np.linalg.eigvals(test_data)
            test_data = test_data[match(ref_data, test_data)]

        diff = np.abs(ref_data - test_data)
        maxdiff = diff.max()
        magn = np.maximum(np.abs(ref_data).max(), np.finfo(float).tiny)
        print(f"{obs + (' * ' if cmp_eigenvalues else ' '):.<45s}",
              f"{maxdiff:10.3e} {maxdiff/magn:10.3e}")

        thr = max(1e-9, 1e-9 * magn)
        if (maxdiff <= thr):
            continue

        any_nonzero = ref_data.astype(bool) + test_data.astype(bool)

        # Something is wrong, print everything
        if any_nonzero.ndim:
            print("-----")
            print(f"{'index':20s} {'reference value':>25s} {'this test':>25s}")
            for i in zip(*any_nonzero.nonzero()):
                idx = ",".join(map(lambda i: "%5s" % i, i))
                ref_val = ref_data[i]
                test_val = test_data[i]
                print(f"{idx:20s} {ref_val:25.16g} {test_val:25.16g}",
                      f"{ref_val - test_val:25.16g}")

            if not any_nonzero.all():
                print("(all other elements are exactly zero in both cases)")
            print("-----")

        mismatch.append(obs)

    print(f"{'':=<45s} {'':=<10s} {'':=<10s}")
    if mismatch:
        raise RuntimeError(f"Observables {', '.join(mismatch)} does not match!")
