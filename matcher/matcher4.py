"""
テンプレートを平行移動+回転して重ねる。
最近接3点で先に照合し、高速化をはかる

/Users/matto/Unison/github/analysis/methane/match より移動
"""

import logging
import numpy as np
from collections import defaultdict
import itertools as it


def rot_trans(uv, gv):
    """
    特異値分解により、与えられた点雲同士を最適に重ねるアフィン変換を求める。
    """
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
        d = uv[i] @ R + t - gv[i]
        I += d @ d
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



def matcher4(gatoms, gcell, uatoms, rc=0.35*1.2, rmsd_max=0.1):
    """
    gatoms: マッチされる点雲、セル内相対座標
    uatoms: テンプレートの点雲絶対座標, in nm
    uatoms[0]はテンプレートの中心にあり、原点でなければならない。(点雲同士の位置あわせのため)
    gcell: セル行列, in nm
    rc: 最近接とみなす最大距離 in nm
    """

    logger = logging.getLogger()
    ng = len(gatoms)
    nu = len(uatoms)

    # suRはuatomsのラベル、原点に近い順。
    suR = np.linalg.norm(uatoms,axis=1)
    uorder = np.argsort(suR)
    # 原子0が原点になければならない。
    assert uorder[0] == 0
    i,j,k = uorder[1:4]
    udeltas = [np.linalg.norm(uatoms[i]-uatoms[j]),
               np.linalg.norm(uatoms[j]-uatoms[k]),
               np.linalg.norm(uatoms[k]-uatoms[i])]
    logger.info(udeltas)
    #rproxは近接距離。グリッドのサイズ
    rprox = 0.4 # nm
    #まず重ねる点を決める。
    for gcenter in range(ng):
        # 進捗を表示
        if int(gcenter*20//ng) != int((gcenter-1)*20/ng):
            logger.info("Progress {0}%".format(int(gcenter*100//ng)))

        # 原点をgcenterとして、gの点を平行移動する。
        sgatoms = gatoms - gatoms[gcenter]
        sgatoms -= np.floor(sgatoms+0.5)
        # 絶対座標に。
        sgatoms = sgatoms @  gcell

        # 住所録(グリッド位置表)を作る
        addressbook = address(sgatoms, rprox)
        # gの点を原点に近い順にソートする。
        sgR = np.linalg.norm(sgatoms,axis=1)
        gnei = set(np.argwhere(sgR<rc).reshape(-1))
        # logger.info(gnei)
        gnei -= {gcenter}
        # logger.info(len(gnei))
        # 2022
        # テンプレート側では、中心に最も近い3点を代表点とし、それらの間の距離を計算しておく。
        # サンプル側では、中心に近いすべての点のなかから、3つを任意に選び、それらの相互距離が
        # テンプレートの3点間距離に近い場合に、照合を先に進める。

        # uの原点最近接3点を、groの最近接3点に重ねる。(6通り)
        for ni,nj,nk in it.permutations(gnei, 3):
            gdeltas = [np.linalg.norm(sgatoms[ni]-sgatoms[nj]),
                       np.linalg.norm(sgatoms[nj]-sgatoms[nk]),
                       np.linalg.norm(sgatoms[nk]-sgatoms[ni])]
            if not np.allclose(udeltas, gdeltas, rtol=5e-2):
                # logger.info(f"{ni} {nj} {nk}\nU {udeltas}\nG {gdeltas}\n")
                continue
            glist = [gcenter, ni,nj,nk]
            #その状態で並進と回転を最適化する。
            uv = uatoms[uorder[0:4]]
            gv = sgatoms[glist]
            #uvをRだけ回転しtだけ並進するとgvに重なる。
            R,t,rmsd = rot_trans(uv,gv)
            # 残りの点を1つずつ対応させる。
            for N in range(4,nu):
                #次のuはソートされたリストから選ぶ。
                uv = uatoms[uorder[0:N+1]]
                #gは近接点のなかからさがす。重複してもよい。
                uvRt = uv[N] @ R + t
                nearest, dmin = find_nearest(uvRt, rprox, addressbook, sgatoms)
                #print(nearest, dmin)
                glist.append(nearest)
                gv = sgatoms[glist]
                # 1点増やすごとに、アフィン変換を再調整する。
                #uvをRだけ回転しtだけ並進するとgvに重なる。
                R,t,rmsd = rot_trans(uv,gv)
            if rmsd < rmsd_max:
                # 点の対応表。
                # corr[i] は、uのi番目の点に対応するgの点のラベル
                corr = [0]*nu
                for i in range(nu):
                    corr[uorder[i]] = glist[i]
                ucenter = 0
                yield (rmsd, gcenter, ucenter, R, corr)
