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
    center = (box[0] + box[1] + box[2])/2
    ori = origin - center
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

s = ""
s += yp.Color(2)
s += yp.Layer(1)
origin = np.zeros(3)
s += drawbox(origin,cell,halfshift=False)
#in absolute coord
unitatoms = np.dot(unitatoms, unitcell)
#s += unitatoms)
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
    matched |= set(members)
    if every != 0 and nline % every == 0:
        s += yp.Color(palette[color])
        s += yp.Layer(1)
        s += drawbox(origin,box)

        s += yp.Layer(2)
        s += drawbox2(origin,box)

        s += yp.Layer(3)
        s += yp.Size(0.15)
    
        s += yp.Color(4)
        s += drawatoms(uatoms)
        for i in range(len(uatoms)):
            s += yp.Line(uatoms[i], atoms[members[i]])

s += yp.Size(0.1)
s += yp.Color(4)
s += yp.Layer(3)
s += drawatoms(atoms, members=matched)
print(s) # end of frame
