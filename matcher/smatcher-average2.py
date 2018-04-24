#!/usr/bin/env python
# coding: utf-8
# smatcherの結果を平均して単位胞内の密度分布関数を作る。
# 次にCrystalAnalysis/symm.pyで対称性分析を行う。

# written completely. for slide-matcher2

# data format: i j radius dx dy dz rmsd

import sys
import numpy as np
import yaplotlib as yp
import pairlist as pl
from collections import defaultdict
import logging

def LoadGRO(file):
    line = file.readline()
    natom = int(file.readline())
    Oatoms = []
    Hatoms = []
    for i in range(natom):
        line = file.readline()
        atomname = line[10:15].replace(" ", "")
        cols = line[20:].split()
        if atomname[0] == "O":
            Oatoms.append([float(x) for x in cols[:3]])
        if atomname[0] == "H":
            Hatoms.append([float(x) for x in cols[:3]])
    box = np.array([float(x) for x in file.readline().split()[:3]])
    return box,np.array(Oatoms),np.array(Hatoms)
        
def loadBOX9(file):
    line = file.readline()
    A = np.array([float(x) for x in file.readline().split()])
    B = np.array([float(x) for x in file.readline().split()])
    C = np.array([float(x) for x in file.readline().split()])
    return A,B,C


def distribute(gridsize, centers, rOatoms, cellmat, uniti, g, gh=None, rHatoms=None):
    grid = np.zeros(gridsize)
    for i in centers: # center
        # 周囲の分子を出力する。
        for j in g[i]:
            d = rOatoms[j] - rOatoms[i]
            d -= np.floor(d+0.5)      #wrap
            pos = np.dot(d, cellmat)  #rel to abs in cell
            rpos = np.dot(pos, uniti) #unit-cell relative
            if -0.5 < rpos[0] < 0.5:
                if -0.5 < rpos[1] < 0.5:
                    if -0.5 < rpos[2] < 0.5:
                        rpos -= np.floor(rpos)
                        ipos = [int(x) for x in np.floor(rpos*np.array(gridsize))]
                        grid[ipos[0], ipos[1], ipos[2]] += 1
    if gh is None:
        return grid
    gridH = np.zeros(gridsize)
    for i in centers: # center
        # O周囲のH分子を出力する。
        for j in gh[i]:
            d = rHatoms[j] - rOatoms[i]
            d -= np.floor(d+0.5)      #wrap
            pos = np.dot(d, cellmat)  #rel to abs in cell
            rpos = np.dot(pos, uniti) #unit-cell relative
            if -0.5 < rpos[0] < 0.5:
                if -0.5 < rpos[1] < 0.5:
                    if -0.5 < rpos[2] < 0.5:
                        rpos -= np.floor(rpos)
                        ipos = [int(x) for x in np.floor(rpos*np.array(gridsize))]
                        gridH[ipos[0], ipos[1], ipos[2]] += 1
    return grid,gridH
    


def slid_add(a,b,sli=None):
    """
    Slide b to maximize the overlap, and add
    """
    logger = logging.getLogger()
    if a is None:
        return b, (0.,0.,0.)
    assert a.shape == b.shape
    if sli is None:
        vmax   = 0
        for ix in range(a.shape[0]):
            c = np.roll(b,ix,axis=0)
            for iy in range(a.shape[1]):
                d = np.roll(c,iy,axis=1)
                for iz in range(a.shape[2]):
                    e = np.roll(d,iz,axis=2)
                    v = np.sum(a*e)
                    # logger.info("v {0}".format(v))
                    if vmax < v:
                        vmax = v
                        emax = e
                        sli  = (ix,iy,iz)
    else:
        emax = np.roll(np.roll(np.roll(b,sli[0],axis=0),sli[1],axis=1),sli[2],axis=2)
    logger.info("Slid {0} {1} {2}".format(*sli))
    return a+emax, sli


def neighborlist(z, hetero=False):
    nei=dict()
    if hetero:
        for i,j in z:
            if i not in nei:
                nei[i] = []
            nei[i].append(j)
    else:
        for i,j in z:
            if i not in nei:
                nei[i] = []
            nei[i].append(j)
            if j not in nei:
                nei[j] = []
            nei[j].append(i)
    return nei

def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    logger = logging.getLogger()
    gro = open(sys.argv[1])
    unitinfo = open(sys.argv[2])
    logger.info("Loading {0}".format(sys.argv[1]))
    box, Oatoms, Hatoms = LoadGRO(gro)
    rOatoms = Oatoms/box
    rHatoms = Hatoms/box
    cellmat = np.diag(box)
    celli = np.linalg.inv(cellmat)
    a,b,c = loadBOX9(unitinfo)
    unitcell = np.vstack((a,b,c))
    # print(unitcell)
    uniti = np.linalg.inv(unitcell)
    D = np.linalg.norm((a+b+c)/2)+0.3 # diag of the half cell.
    grid = pl.determine_grid(cellmat, D)
    # print(grid)
    logger.info("Making the neighbor list 1 (O).")
    g = neighborlist( pl.pairs_fine(rOatoms, D, cellmat, grid, distance=False))
    logger.info("Making the neighbor list 2 (H).")
    gh = neighborlist(pl.pairs_fine_hetero(rOatoms, rHatoms, D, cellmat, grid, distance=False), hetero=True)
    logger.info("Reading smatch file.")
    matched = defaultdict(list)
    while True:
        line = sys.stdin.readline() #expect *.smatch result of slide-matcher2
        if len(line) == 0:
            break
        cols = line.split()
        if len(cols)<7:
            break # broken line; maybe end of the file.
        i,j = [int(x) for x in cols[:2]]
        rad, dx,dy,dz, rmsd = [float(x) for x in cols[2:]]
        # d = np.array([dx,dy,dz])
        # dL = np.linalg.norm(d)
        if rmsd < 0.085: # 0.08 for T144
            # 距離によらず、よく一致したペアのラベルを保存しておく。
            matched[i].append(j)
            matched[j].append(i)
    # 一番いろんなところと似ていた原子kをさがす。
    #上位10個について、それぞれ周辺原子分布を作り、平行移動して重ねる。
    # うまくいくかどうかはわからない。
    logger.info("Overlaying atomic configurations.")
    keys = matched.keys()
    logger.info(keys)
    order = sorted(keys, key=lambda x:len(matched[x]), reverse=True)
    logger.info([len(matched[x]) for x in order[:10]])
    gridsize = (32,32,12)
    sumgrid = None
    sumgridH = None
    accum = 0
    for k in order[:10]:
        matched[k].append(k) # add itself
        kL = len(matched[k])
        logger.info("Champion: {0} with {1} neighbors.".format(k, kL))
        grid, gridH = distribute(gridsize, matched[k], rOatoms, cellmat, uniti, g, gh=gh, rHatoms=rHatoms)
        sumgrid, sli  = slid_add(sumgrid, grid)
        sumgridH, sli = slid_add(sumgridH, gridH, sli)
        accum += kL
    logger.info("{0} configurations collected.".format(accum))
    # ここである程度対称性を整えてから出力するとあとが楽。
    print("@GRID")
    print(*gridsize)
    for i in range(gridsize[0]):
        for j in range(gridsize[1]):
            for k in range(gridsize[2]):
                print(sumgrid[i,j,k] / accum)
    print("@GRID")
    print(*gridsize)
    for i in range(gridsize[0]):
        for j in range(gridsize[1]):
            for k in range(gridsize[2]):
                print(sumgridH[i,j,k] / accum)
    logger.info("Done.")
    

def test():
    a = np.zeros((5,5))
    b = np.zeros
if __name__ == "__main__":
    main()
