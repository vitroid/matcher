#!/usr/bin/env python3

import numpy as np
import sys
import yaplotlib as yp

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
    cell  = np.diag(dimen)*10 #in angstrom
    return cell, Os



#Put different colors for different vectors
def direction2color(v, digitize=4):
    x = np.array([1.0,0.0,0.0])
    y = np.array([0.0,1.0,0.0])
    z = np.array([0.0,0.0,1.0])
    e = v / np.linalg.norm(v)
    R = int(abs(np.dot(x,e))*digitize)
    G = int(abs(np.dot(y,e))*digitize)
    B = int(abs(np.dot(z,e))*digitize)
    return R,G,B


def drawbox(origin, box, halfshift=True):
    s = ""
    if halfshift:
        center = (box[0] + box[1] + box[2])/2
        ori = origin - center
    else:
        ori = origin
    for v in (np.zeros(3), box[1], box[2], box[1]+box[2]):
        s += yp.Line(v+ori, v+ori+box[0])
    for v in (np.zeros(3), box[0], box[2], box[0]+box[2]):
        s += yp.Line(v+ori, v+ori+box[1])
    for v in (np.zeros(3), box[0], box[1], box[0]+box[1]):
        s += yp.Line(v+ori, v+ori+box[2])
    return s

def drawbox2(origin, box):
    s = ""
    # center = (box[0] + box[1] + box[2])/2
    ori = origin #- center
    s += yp.Color(2)
    s += yp.Polygon([ori,ori+box[0],ori+box[1]])
    s += yp.Color(3)
    s += yp.Polygon([ori,ori+box[1],ori+box[2]])
    s += yp.Color(4)
    s += yp.Polygon([ori,ori+box[2],ori+box[0]])
    return s

def drawatoms(atoms, members=None):
    s = ""
    if members is None:
        for a in atoms:
            s += yp.Circle(a)
    else:
        for i in members:
            s += yp.Circle(atoms[i])
    return s

def usage():
    print("usage: {0} [-a][-e every][-v maxval] traj.gro pattern.gro < matcher.output")
    sys.exit(1)
    
def main():
    every = 1
    maxval = 1.0
    adjdens = False
    while sys.argv[1][0] == "-":
        if sys.argv[1] == "-e":
            every = int(sys.argv[2])
            sys.argv.pop(1)
            sys.argv.pop(1)
        elif sys.argv[1] == "-v":
            maxval = float(sys.argv[2])
            sys.argv.pop(1)
            sys.argv.pop(1)
        elif sys.argv[1] == "-a":
            adjdens = True
            sys.argv.pop(1)
        else:
            usage()

    Cell, Oatoms = LoadGRO(open(sys.argv[1]))  # in AA, in AA
    Unitcell, Unitatoms = LoadGRO(open(sys.argv[2])) # in AA, in rel
    dens0 = Oatoms.shape[0] / np.linalg.det(Cell)
    dens1 = Unitatoms.shape[0] / np.linalg.det(Unitcell)
    ratio = (dens1/dens0)**(1./3.)
    if adjdens:
        Unitcell *= ratio

    mode = ""
    if len(sys.argv) > 3:
        mode = sys.argv[3]

    s = ""
    s += yp.Color(2)
    s += yp.Layer(1)
    origin = np.zeros(3)
    s += drawbox(origin,Cell,halfshift=False)
    #in absolute coord
    # unitatoms = np.dot(unitatoms, Unitcell)
    #s += unitatoms)

    #s = ""
    palette = dict()
    matched = set()
    nline = 0
    while True:
        line = sys.stdin.readline()
        if len(line) == 0:
            break
        nline += 1
        #parse the line
        cols = line.split()
        if len(cols) < 13:
            break
        # 2018-4-9 New output format of matcher.c
        msd     = float(cols[0])
        Origin  = Oatoms[int(cols[1])].copy()  #atom at the matching center
        # Origin -= np.floor(Origin/Cell+0.5)*Cell
        celli = np.linalg.inv(Cell)
        relorigin = np.dot(Origin, celli)
        relorigin -= np.floor(relorigin)
        Origin = np.dot(relorigin, Cell)
        center  = int(cols[2])
        rotmat  = np.array([float(x) for x in cols[3:12]]).reshape((3,3))
        N = int(cols[12])
        if len(cols) < 13+N:
            break
        irot = np.linalg.inv(rotmat)
        members = [int(x) for x in cols[13:N+13]]
        #draw matched box

        # roll the unit cell to center (abs)
        rel = Unitatoms - Unitatoms[center]
        rel -= np.dot(np.floor(np.dot(rel, np.linalg.inv(Unitcell))+0.5),Unitcell)
        # rel to abs
        Slid = rel
        # rotate box
        RotUnitcell       = np.dot(Unitcell, rotmat)
        # 
        Boxslide = np.dot(-Unitatoms[center], rotmat) + Origin
        # rotate atoms in the unit cell
        Slidunit    = np.dot(Slid, rotmat) + Origin
        #s += yp.Color(3)
        if mode == "R":
            color = direction2color(rotmat[0]+rotmat[1]+rotmat[2])
        elif mode == "T":
            color = direction2color(rotmat[2])
        else:
            color = (0,3,0) #green
        if color not in palette:
            palette[color] = len(palette)+5
            s += yp.SetPalette(palette[color],color[0]*255//3,color[1]*255//3,color[2]*255//3)
        if msd < maxval:
            matched |= set(members)
        if every != 0 and nline % every == 0 and msd < maxval:
            s += yp.Color(palette[color])
            # matched box
            s += yp.Layer(4)
            s += drawbox(Origin,RotUnitcell,halfshift=True)
            # unit cell
            s += yp.Layer(2)
            s += drawbox(Boxslide,RotUnitcell,halfshift=False)

            #s += yp.Layer(2)
            #s += drawbox2(Boxslide,RotUnitcell)

            #s += yp.Layer(3)
            #s += yp.Size(0.4)

            #s += yp.Color(4)
            #s += drawatoms(Slidunit)
            #for i in range(len(Slidunit)):
            #    g = Oatoms[members[i]] - Origin
            #    g -= np.floor(g/Cell+0.5)*Cell
            #    g += Origin
            #    d = Slidunit[i] - g
            #    d -= np.floor(d/Cell+0.5)*Cell
            #    s += yp.Line(g, g+d)
            #    # 変位ベクトルはどうやってもうまくいかないので、やめる。

    s += yp.Size(0.4)
    s += yp.Color(4)
    s += yp.Layer(3)
    s += drawatoms(Oatoms, members=matched)
    print(s) # end of frame

if __name__ = "__main__":
    main()
