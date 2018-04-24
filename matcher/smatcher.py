#!/usr/bin/env python


import sys
import numpy as np
from matcher.csmatcher import smatcher


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
    Os = np.array(Os) #in angstrom
    #the last line is cell shape
    line  = file.readline()
    cols  = line.split()
    dimen = np.array([float(x) for x in cols])
    cell  = dimen #in angstrom
    return cell, Os


def hook1(lattice):
    """
    Slide and match. A format plugin for GenIce.
    
    usage: genice 7 -f smatcher[radius:rmsdmax]

    matcher plugin accepts several parameters:
    1st [0.8]    : Radius of matching. (0.8 nm)
    2nd [0.06]   : Maximum rmsd between ice and template structures. (0.06 nm)

    """
    lattice.logger.info("Hook1: Translational matching.")
    cellmat = lattice.repcell.mat
    # simulation box must be orthogonal.
    assert cellmat[1,0] == 0 and cellmat[2,0] == 0 and cellmat[2,1] == 0
    assert cellmat[0,2] == 0 and cellmat[0,1] == 0 and cellmat[1,2] == 0
    cell = np.array([cellmat[0,0], cellmat[1,1], cellmat[2,2]])
    positions = lattice.reppositions - np.floor(lattice.reppositions)
    apos = np.dot(positions, cellmat)
    N = len(positions)**2
    every = 1
    if N > 100000:
        #//evaluate only 10 M pairs.
        every = N//(2*100000)+1
        lattice.logger.info("Too many combinations! Will skip every {0} to reduce time.".format(every))

    smatches = smatcher(apos, cell, lattice.smatcher_radius, lattice.smatcher_rmsdmax, every)
    for smatch in smatches:
        i, j, rad, d, rmsd = smatch
        s = "{0} {1} {2:.4f} ".format(i, j, rad)
        for i in d:
            s += "{0:.5f} ".format(i)
        s += "{0:.5f} ".format(rmsd)
        print(s)
    lattice.logger.info("Hook1: end.")

    

def hook0(lattice, arg):
    lattice.smatcher_radius  = 0.8  # nm
    lattice.smatcher_rmsdmax = 0.06 # nm
    if arg != "":
        cols = arg.split(":")
        lattice.smatcher_radius  = float(cols[0])
        if len(cols)>1:
            lattice.smatcher_rmsdmax = float(cols[1])




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


hooks = {0:hook0, 1:hook1}
