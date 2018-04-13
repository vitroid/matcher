#!/usr/bin/env python

__version__ = "0.1.1"

import sys
import numpy as np
from cmatcher import matcher


#note: this is modified. do not reuse.
def LoadGRO(file, rel=False):
    line = file.readline()
    natom = int(file.readline())
    Os = []
    for i in range(natom):
        line = file.readline()
        cols = line[20:].split()[:3]
        xyz = [float(x) for x in cols]
        if " O" in line:
            Os.append(xyz)
    Os = np.array(Os) #in nm
    #the last line is cell shape
    line  = file.readline()
    cols  = line.split()
    cell  = np.array([float(x) for x in cols]) # in nm
    if rel:
        Os /= cell
        Os -= np.floor(Os)
    return cell, Os




def hook1(lattice):
    """
    Match a given small unit cell with an ice structure. A format plugin for GenIce.
    
    usage: genice 7 -f matcher[1c.gro:0.03:0.06:1]

    matcher plugin accepts several parameters:
    1st [1c.gro] : the template structure.
    2nd [0.03]   : Tolerance for geometrical hashing; Less than 0.03 but non-zero is relevant.
    3rd [0.06]    : Maximum rmsd between ice and template structures. [0.06 nm]
    4th [1]      : Size of the template will be scaled to adjust the density when it is 1.

    3D Geometric Hashing is utilized.
    https://en.wikipedia.org/wiki/Geometric_hashing
    Nussinov, R.; Wolfson, H. J., ProNAS 88, 10495–10499 (1991). doi:10.1073/pnas.88.23.10495. 
    """
    lattice.logger.info("Hook1: Template matching.")
    cellmat = lattice.repcell.mat
    # simulation box must be orthogonal.
    assert cellmat[1,0] == 0 and cellmat[2,0] == 0 and cellmat[2,1] == 0
    assert cellmat[0,2] == 0 and cellmat[0,1] == 0 and cellmat[1,2] == 0
    cell = np.array([cellmat[0,0], cellmat[1,1], cellmat[2,2]])
    positions = lattice.reppositions - np.floor(lattice.reppositions)
    apos = np.dot(positions, cellmat)
    matches = matcher(apos, cell, unitatoms, unitcell, err, rmsdmax, adjdens)
    for match in matches:
        rmsd, a, au, mat, atoms = match
        s = "{0:.5f} {1} {2} ".format(rmsd, a, au)
        for i in mat:
            s += "{0:.5f} ".format(i)
        s += "{0} ".format(len(atoms))
        for i in atoms:
            s += "{0} ".format(i)
        print(s)
    lattice.logger.info("Hook1: end.")

    

def argparser(arg):
    global err, rmsdmax, adjdens, unitatoms, unitcell
    cols = arg.split(":")
    unitcell, unitatoms = LoadGRO(open(cols[0]), rel=True)
    err = 0.03
    rmsdmax = 0.06 # nm
    adjdens = False
    if len(cols)>1:
        err = float(cols[1])
        if len(cols)>2:
            rmsdmax = float(cols[2])
            if len(cols)>3:
                adjdens = cols[3] in ("1", "TRUE", "True", "T", "t")
    
    
def test():
    """
    emulate matcher.c
    """
    cell, Oatoms = LoadGRO(open(sys.argv[1]))
    unitcell, runitatoms = LoadGRO(open(sys.argv[2]), rel=True)
    err = float(sys.argv[3])
    rmsdmax = float(sys.argv[4])
    adjdens = False
    matches = matcher(Oatoms, cell, runitatoms, unitcell, err, rmsdmax, adjdens)
    print(matches)


if __name__ == "__main__":
    test()


hooks = {1:hook1}
