#!/usr/bin/env python3

import numpy as np
import sys
from genice import yaplotlib as yp

def LoadGRO(file):
    headline = file.readline()
    natom = int(file.readline())
    atoms = []
    extras   = []
    for i in range(natom):
        line = file.readline()
        cols = line[20:].split()[:3]
        xyz = [float(x) for x in cols]
        extra = line[:20]
        extras.append(extra)
        atoms.append(xyz)
    #the last line is cell shape
    line  = file.readline()
    cols  = line.split()
    dimen = np.array([float(x) for x in cols])
    if dimen.shape[0] == 9:
        cell = [[dimen[0], dimen[3], dimen[4]],
                [dimen[5], dimen[1], dimen[6]],
                [dimen[7], dimen[8], dimen[2]]]
        cell  = np.array(cell)
    else:
        cell = np.diag(dimen)
    return headline, cell, atoms, extras


def ToGRO(headline, cell, atoms, extras):
    s = headline
    s += "{0}\n".format(len(atoms))
    for position,extra in zip(atoms, extras):
        s += "{0}{1:8.3f}{2:8.3f}{3:8.3f}\n".format(extra,position[0],position[1],position[2])
    s += "    {0} {1} {2}\n".format(cell[0,0],cell[1,1],cell[2,2])
    return s


headline, cell, atoms, extras = LoadGRO(sys.stdin)
#from absolute to relative
atoms = np.dot(atoms, np.linalg.inv(cell))
x,y,z = cell
xL = np.dot(x,x)**0.5
yL = np.dot(y,y)**0.5
zL = np.dot(z,z)**0.5
#assume the cell is almost orthogonal.
cell = np.diag((xL,yL,zL))
atoms = np.dot(atoms, cell)
print(ToGRO(headline, cell, atoms, extras), end="")
