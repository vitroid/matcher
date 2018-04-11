#!/usr/bin/env python
# coding: utf-8
# slide and overlay the lattice manually.


# written completely. for slide-matcher2

# data format: i j radius dx dy dz rmsd

import sys
import numpy as np
import yaplotlib as yp
import pairlist
import networkx as nx

def LoadGRO(file):
    line = file.readline()
    natom = int(file.readline())
    atoms = []
    for i in range(natom):
        line = file.readline()
        atomname = line[10:15].replace(" ", "")
        cols = line[20:].split()
        if atomname[0] == "O":
            atoms.append([float(x) for x in cols[:3]])
    box = np.array([float(x) for x in file.readline().split()[:3]])
    return box,np.array(atoms)
        

def main():
    gro = open(sys.argv[1])
    box, atoms = LoadGRO(gro)
    ratoms = atoms/box
    cellmat = np.diag(box)
    page1 = ""
    page2 = ""
    page3 = ""
    page4 = ""
    page1 += yp.Size(0.05)
    grid = None
    while True:
        line = sys.stdin.readline() #expect *.smatch result of slide-matcher2
        if len(line) == 0:
            break
        cols = line.split()
        if len(cols)<7:
            break # broken line; maybe end of the file.
        i,j = [int(x) for x in cols[:2]]
        rad, dx,dy,dz, rmsd = [float(x) for x in cols[2:]]
        if grid is None:
            grid = pairlist.determine_grid(cellmat, rad)
            lastrad = rad
            p,q,r = pairlist.pairs_fine(ratoms, rad, cellmat, grid, distance=True, raw=True)
            g = nx.Graph()
            for k,pq in enumerate(zip(p,q)):
                p,q = pq
                g.add_edge(p,q,length=r[k])
        else:
            assert lastrad == rad
        d = np.array([dx,dy,dz])
        if rmsd < 0.09: # 0.09 for hyper ice T
            # page 1: displacement vectors
            page1 += yp.Line(atoms[i], atoms[i]+d)
            page1 += yp.Circle(atoms[j])
            # pages 2: displacement vectors (centered)
            page2 += yp.Size(0.05)
            page2 += yp.Line(np.zeros(3), d)
            page2 += yp.Circle(d)
            page2 += yp.Text(d, "{0} {1}".format(i,j))
            # page 3..: matching
            s = ""
            s += yp.Size(0.05)
            s += yp.Layer(1)
            s += yp.Color(3)
            for ni in g[i]:
                s += yp.Circle(atoms[ni])
            s += yp.Layer(2)
            s += yp.Color(4)
            for nj in g[j]:
                s += yp.Circle(atoms[nj]-atoms[j]+atoms[i])
            page3 += s
            s = ""
            s += yp.Size(0.05)
            s += yp.Layer(1)
            s += yp.Color(3)
            for ni in g[i]:
                s += yp.Circle(atoms[ni]-atoms[i])
            s += yp.Layer(2)
            s += yp.Color(4)
            for nj in g[j]:
                s += yp.Circle(atoms[nj]-atoms[j])
            page4 += s
            page4 += yp.NewPage()

    print(page1)
    print(page2)
    print(page3)
    print(page4)
    

if __name__ == "__main__":
    main()
