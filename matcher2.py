#!/usr/bin/env python

import sys
import numpy as np
import yaplotlib as yp
from collections import defaultdict
import itertools as it

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

def main():
    gcell, gatoms = LoadGRO(open(sys.argv[1]), rel=True)
    ucell, uatoms = LoadGRO(open(sys.argv[2]), rel=True)
    ng = len(gatoms)
    nu = len(uatoms)
    #uRは単位胞の外接球の半径
    uR = sum(ucell)/2.0
    #rproxは近接距離。グリッドのサイズ
    rprox = 0.4
    #まず重ねる点を決める。
    ucenter = 145
    for gcenter in range(ng):
        # それらの点を原点に平行移動する。
        sgatoms = gatoms - gatoms[gcenter]
        sgatoms -= np.floor(sgatoms+0.5)
        suatoms = uatoms - uatoms[ucenter]
        suatoms -= np.floor(suatoms+0.5)
        # 絶対座標に変換する
        sgatoms = sgatoms * gcell # orthorhombic, absolute position
        suatoms = suatoms * ucell # orthorhombic, absolute position
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
                print("{0:.5f}".format(rmsd), gcenter, ucenter, end=" ")
                for i in range(3):
                    for j in range(3):
                        print("{0:.5f}".format(R[i,j]), end=" ")
                corr = [0]*nu
                for i in range(nu):
                    corr[uorder[i]] = glist[i]
                print(nu, *corr)

if __name__ == "__main__":
    main()
    