#!/usr/bin/env python3

import numpy as np
import sys
from genice import yaplotlib as yp

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


def LoadAR3R(file):
    while True:
        line = file.readline()
        if len(line) == 0:
            break
        if "@BOX3" == line[:5]:
            line  = file.readline()
            cols  = line.split()
            dimen = np.array([float(x) for x in cols])
            cell  = np.diag(dimen)
        elif "@AR3R" == line[:5]:
            natom = int(file.readline())
            Os = []
            for i in range(natom):
                line = file.readline()
                cols = line.split()
                xyz = np.array([float(x) for x in cols])
                xyz -= np.floor( xyz + 0.5 )
                Os.append(xyz)
            Os = np.array(Os)
    #the last line is cell shape
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
    #print(box)
    if halfshift:
        center = (box[0] + box[1] + box[2])/2
        ori = origin - center
    else:
        ori = origin
    for v in (np.zeros(3), box[1], box[2], box[1]+box[2]):
        print(yp.Line(v+ori, v+ori+box[0]), end="")
        #s += yp.Line(v+ori, v+ori+box[0])
    for v in (np.zeros(3), box[0], box[2], box[0]+box[2]):
        print(yp.Line(v+ori, v+ori+box[1]), end="")
        #s += yp.Line(v+ori, v+ori+box[1])
    for v in (np.zeros(3), box[0], box[1], box[0]+box[1]):
        print(yp.Line(v+ori, v+ori+box[2]), end="")
        #s += yp.Line(v+ori, v+origin+box[2])

def drawbox2(origin, box):
    center = (box[0] + box[1] + box[2])/2
    ori = origin - center
    print(yp.Color(2), end="")
    print(yp.Polygon([ori,ori+box[0],ori+box[1]]), end="")
    print(yp.Color(3), end="")
    print(yp.Polygon([ori,ori+box[1],ori+box[2]]), end="")
    print(yp.Color(4), end="")
    print(yp.Polygon([ori,ori+box[2],ori+box[0]]), end="")

def drawatoms(atoms, members=None):
    if members is None:
        for a in atoms:
            print(yp.Circle(a),end="")
    else:
        for i in members:
            print(yp.Circle(atoms[i]),end="")

every = 1
if sys.argv[1] == "-e":
    every = int(sys.argv[2])
    sys.argv.pop(1)
    sys.argv.pop(1)
    
cell, atoms = LoadGRO(open(sys.argv[1]))
unitcell, unitatoms = LoadAR3R(open(sys.argv[2]))
mode = "R"
if len(sys.argv) > 3:
    mode = sys.argv[3]

print(yp.Color(2),end="")
print(yp.Layer(1),end="")
origin = np.zeros(3)
drawbox(origin,cell,halfshift=False)
#in absolute coord
unitatoms = np.dot(unitatoms, unitcell)
#print(unitatoms)
dmin = 1e99
orig = None
for a in unitatoms:
    L = np.dot(a,a)
    if L < dmin:
        dmin = L
        orig = a
    
#s = ""
palette = dict()
matched = set()
for nline,line in enumerate(sys.stdin):
    #parse the line
    cols = line.split()
    N = int(cols[0])
    members = [int(x) for x in cols[1:N+1]]
    msd     = float(cols[N+1])
    origin  = atoms[int(cols[N+2])].copy()  #atom at the matching center
    rotmat  = np.array([float(x) for x in cols[N+3:N+12]]).reshape((3,3))
    #draw matched box
    boxorigin = np.dot(orig, rotmat)
    origin   -= boxorigin #corner of the cell
    box       = np.dot(unitcell, rotmat)
    uatoms    = np.dot(unitatoms, rotmat) + origin
    #print(yp.Color(3),end="")
    if mode == "R":
        color = direction2color(rotmat[0]+rotmat[1]+rotmat[2])
    elif mode == "T":
        color = direction2color(rotmat[2])
    else:
        color = (0,3,0) #green
    if color not in palette:
        palette[color] = len(palette)+5
        print(yp.SetPalette(palette[color],color[0]*255//3,color[1]*255//3,color[2]*255//3),end="")
    matched |= set(members)
    if every != 0 and nline % every == 0:
        print(yp.Color(palette[color]), end="")
        print(yp.Layer(1),end="")
        drawbox(origin,box)

        print(yp.Layer(2),end="")
        drawbox2(origin,box)

        print(yp.Layer(3),end="")
        print(yp.Size(0.15),end="")
    
        print(yp.Color(4), end="")
        drawatoms(uatoms)
        for i in range(len(uatoms)):
            print(yp.Line(uatoms[i], atoms[members[i]]),end="")

print(yp.Size(0.1),end="")
print(yp.Color(4),end="")
print(yp.Layer(3),end="")
drawatoms(atoms, members=matched)
