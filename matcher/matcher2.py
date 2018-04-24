#!/usr/bin/env python

# 近くから順にマッチングさせる、順当なやりかた。
# 巨大な情報を準備する必要がなく、かなり速い。
# Cにすればもっと速くなるだろう。


import sys
import logging
import numpy as np
import yaplotlib as yp
from collections import defaultdict
import itertools as it
from matcher import cmatcher2

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
    assert cell.shape[0] == 3
    if rel:
        Os /= cell
        Os -= np.floor(Os+0.5)
    return cell, Os

def rot_trans(uv, gv):
    ucom = np.mean(uv, axis=0)
    gcom = np.mean(gv, axis=0)
    H = np.zeros((3,3))
    for i in range(uv.shape[0]):
        H += (uv[i] -ucom).reshape((3,1))*(gv[i] - gcom).reshape((1,3))
    U,S,Vt = np.linalg.svd(H.T)
    # R = Vt.T * U.T
    R = np.dot(U,Vt).T
    t = -np.dot(ucom,R) + gcom
    I = 0.0
    for i in range(uv.shape[0]):
        d = np.dot(uv[i],R)+t - gv[i]
        I += np.dot(d,d)
    rmsd = (I/uv.shape[0])**0.5
    #print("##", rmsd,R,t)
    return R,t,rmsd


def address(points, rprox):
    ipoints = np.floor(points/rprox)
    addressbook = defaultdict(list)
    for i, ip in enumerate(ipoints):
        addressbook[tuple(ip)].append(i)
    return addressbook


def nearby(addr, addressbook):
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                prox = (addr[0]+i, addr[1]+j, addr[2]+k)
                for resident in addressbook[prox]:
                    yield resident

def find_nearest(pos, rprox, addressbook, points):
    addr = np.floor(pos / rprox)
    dmin = 1e99
    rmin = -1
    for resident in nearby(addr, addressbook):
        d = pos - points[resident]
        dd = np.dot(d,d)
        if dd < dmin:
            dmin = dd
            rmin = resident
    return rmin, dmin



def matcher2_python(gatoms, gcell, uatoms, ucell, adjdens=True, nostore=True):
    """
    gatoms, uatoms: relative
    gcell, ucell: orthorhombic
    """
    logger = logging.getLogger()
    ng = len(gatoms)
    nu = len(uatoms)
    #uRは単位胞の外接球の半径
    uR = np.sum(ucell)/2.0
    #rproxは近接距離。グリッドのサイズ
    rprox = 0.4
    #まず重ねる点を決める。
    ucenter = 0
    matches = []
    for gcenter in range(ng):
        if int(gcenter*100//ng) != int((gcenter-1)*100/ng):
            logger.info("Progress {0}%".format(int(gcenter*100//ng)))
        # それらの点を原点に平行移動する。
        sgatoms = gatoms - gatoms[gcenter]
        sgatoms -= np.floor(sgatoms+0.5)
        suatoms = uatoms - uatoms[ucenter]
        suatoms -= np.floor(suatoms+0.5)
        # 絶対座標に変換する
        sgatoms *= gcell # orthorhombic, absolute position
        suatoms *= ucell # orthorhombic, absolute position
        # 住所録を作る
        addressbook = address(sgatoms, rprox)
        # 原点に近い順にソートする。
        sgR = np.linalg.norm(sgatoms,axis=1)
        gorder = sorted(list(range(ng)), key=lambda x:sgR[x])
        suR = np.linalg.norm(suatoms,axis=1)
        uorder = sorted(list(range(nu)), key=lambda x:suR[x])
        # uの周囲の最近接3点を、groの最近接3点に重ねる。(6通り)
        for ni,nj,nk in it.permutations(gorder[1:4]):
            glist = [gcenter, ni,nj,nk]
            #その状態で並進と回転を最適化する。
            uv = suatoms[uorder[0:4]]
            gv = sgatoms[glist]
            #uvをRだけ回転しtだけ並進するとgvに重なる。
            R,t,rmsd = rot_trans(uv,gv)
            for N in range(4,nu):
                #次のuはソートされたリストから選ぶ。
                uv = suatoms[uorder[0:N+1]]
                #gは近接点のなかからさがす。重複してもよい。
                uvRt = np.dot(uv[N],R)+t
                nearest, dmin = find_nearest(uvRt, rprox, addressbook, sgatoms)
                #print(nearest, dmin)
                glist.append(nearest)
                gv = sgatoms[glist]
                #uvをRだけ回転しtだけ並進するとgvに重なる。
                R,t,rmsd = rot_trans(uv,gv)
                # print("#",gcenter,rmsd)
                if rmsd > 0.1:
                    break
            if rmsd < 0.1:
                corr = [0]*nu
                for i in range(nu):
                    corr[uorder[i]] = glist[i]
                if nostore:
                    print("{0:.5f}".format(rmsd), gcenter, ucenter, end=" ")
                    for i in range(3):
                        for j in range(3):
                            print("{0:.5f}".format(R[i,j]), end=" ")
                    print(nu, *corr)
                else:
                    matches.append((rmsd, gcenter, ucenter, R, corr))
    return matches


def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(message)s")
    gcell, gatoms = LoadGRO(open(sys.argv[1]), rel=True)
    ucell, uatoms = LoadGRO(open(sys.argv[2]), rel=True)
    if len(sys.argv) > 3 and sys.argv[3] == "py":
        match_iter = matcher2_python(gatoms, gcell, uatoms, ucell, adjdens=True, nostore=True)
    else:
        match_iter = cmatcher2.matcher2(gatoms, gcell, uatoms, ucell, True, True)
        
    for match in match_iter:
        rmsd, gcenter, ucenter, R, corr = match
        print("{0:.5f}".format(rmsd), gcenter, ucenter, end=" ")
        for i in range(3):
            for j in range(3):
                print("{0:.5f}".format(R[i,j]), end=" ")
        print(len(corr), corr)
        
def hook1(lattice):
    """
    Slide and match. A format plugin for GenIce.
    
    usage: genice 7 -f smatcher[radius:rmsdmax]

    matcher plugin accepts several parameters:
    1st [0.8]    : Radius of matching. (0.8 nm)
    2nd [0.06]   : Maximum rmsd between ice and template structures. (0.06 nm)

    """
    lattice.logger.info("Hook1: Template matching.")
    cellmat = lattice.repcell.mat
    # simulation box must be orthogonal.
    assert cellmat[1,0] == 0 and cellmat[2,0] == 0 and cellmat[2,1] == 0
    assert cellmat[0,2] == 0 and cellmat[0,1] == 0 and cellmat[1,2] == 0
    gcell = np.array([cellmat[0,0], cellmat[1,1], cellmat[2,2]])
    gatoms = lattice.reppositions - np.floor(lattice.reppositions)
    # apos = np.dot(positions, cellmat)
    adjdens = False
    nostore = True
    matches = cmatcher2.matcher2(gatoms, gcell, lattice.uatoms, lattice.ucell, adjdens, nostore)
    lattice.logger.info("Hook1: end.")

    

def hook0(lattice, arg):
    args = arg.split(":")
    lattice.ucell, lattice.uatoms = LoadGRO(open(args[0]), rel=True)

hooks={0:hook0, 1:hook1}
        
def test():
    H = np.array([2.0, 3.0, 5.0, 1.0, 4.0, 2.0, 3.0, 5.0, 1.0]).reshape((3,3))
    U,S,Vt = np.linalg.svd(H)
    # R = Vt.T * U.T
    R = np.dot(U,Vt)
    print("U")
    print(U)
    print("Vt")
    print(Vt)
    print("R")
    print(R)

    
if __name__ == "__main__":
    main()
