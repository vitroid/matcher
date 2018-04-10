#!/usr/bin/env python

__version__ = "0.1"

import sys
import numpy as np
from csmatcher import smatcher


#note: this is modified. do not reuse.
def LoadGRO(file):
    line = file.readline()
    natom = int(file.readline())
    Os = []
    for i in range(natom):
        line = file.readline()
        cols = line[20:].split()[:3]
        xyz = [float(x) for x in cols]
        if " O" in line:
            Os.append(xyz)
    Os = np.array(Os)*10 #in angstrom
    #the last line is cell shape
    line  = file.readline()
    cols  = line.split()
    dimen = np.array([float(x) for x in cols])
    cell  = dimen*10 #in angstrom
    return cell, Os


def test():
    """
    emulate smatcher.c
    """
    cell, Oatoms = LoadGRO(open(sys.argv[1]))
    radius = float(sys.argv[2])
    rmsdmax = float(sys.argv[3])
    every = 100000
    print(cell.shape);
    print(Oatoms.shape)
    smatch = smatcher(Oatoms, cell, radius, rmsdmax, every)
    print(smatch)


if __name__ == "__main__":
    test()
