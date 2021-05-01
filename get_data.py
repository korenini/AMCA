import sys
import numpy as np
import itertools
import settings as pset


data_files = pset.data_files

data = []
dist = []

for f in data_files:
    try:
        raw = np.genfromtxt(f, delimiter=pset.delimiter)
    except (IOError, OSError):
        print("Data file(s) not found.")
        sys.exit(1)
    if raw.shape.__len__() == 2:
        pass
    elif raw.shape.__len__() == 1:
        raw.shape = (raw.size, 1)
    else:
        raise ValueError('Data not properly formatted.')

    data.append(raw)

    dist_euclid = np.zeros([len(raw), len(raw)], dtype=np.float64)

    for i in itertools.combinations(range(len(raw)), 2):
        dist_pairs = raw[i[0],] - raw[i[1],]
        dist_euclid[i[0],i[1]] = dist_euclid[i[1],i[0]] = (dist_pairs*dist_pairs).sum()

    dist.append(dist_euclid)


npdist = np.array(dist)


# Make distance array read-only
npdist.setflags(write=False)

# Number of criteria
ncrits = len(data)
mcrits = tuple(range(len(pset.data_files)))

# Number of cases
nunits = len(npdist[0])