#!/usr/bin/env python3

#make images
#a point consists of three tuples
#a tuple contains an axisname, its sign, and offset


import sympy
from math import floor

def screwz(points):
    for p in points:
        x,y,z = p
        yield (x,y,z)
        yield (y,1/2-x,1/4+z)
        yield(1/2-x,1/2-y,1/2+z)
        yield(1/2-y,x,3/4+z)


def rotz(points):
    for p in points:
        x,y,z = p
        yield (x,y,z)
        yield (-x,-y,z)
        
def roty(points):
    for p in points:
        x,y,z = p
        yield (x,y,z)
        yield (1/2-x,y,-z)
        
def rotx(points):
    for p in points:
        x,y,z = p
        yield (x,y,z)
        yield (x,1/4-y,-z)

def center(points):
    for p in points:
        x,y,z = p
        yield (x,y,z)
        yield (-x,1/2-y,-z)

def center2(points):
    for p in points:
        x,y,z = p
        yield (x,y,z)
        yield (1/4-x,-y,-z)

def normalize(point):
    def n0(x):
        eq = x.as_coefficients_dict()
        #print(type(sympy.N(eq[1])))
        f = floor(sympy.N(eq[1]))
        if f:
            return x-f
        return x
    x,y,z = point
    return (n0(x),n0(y),n0(z))

x,y,z = sympy.symbols("x y z")
origin = (x,y,z)
points = set([origin])
points = set([p for p in screwz(points)])
#points = set([p for p in rotz(points)])
#points = set([p for p in roty(points)])
points = set([p for p in rotx(points)])
#points = set([p for p in center(points)])
points = set([normalize(p) for p in center2(points)])
points = set([p for p in screwz(points)])
points = set([p for p in rotx(points)])
points = set([normalize(p) for p in center2(points)])
for p in points:
    x,y,z = p
    print("'{0},{1},{2}'".format(x,y,z))
#print(len(points))
        
