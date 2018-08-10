#!/usr/bin/env python

# Copyright 2008, Angel Yanguas-Gil
# Rev. 2014, Ningchuan Xiao
# The revised version does not require Xygraph

# Rev. 2018, sort for python3

# only the classes defined here are accessible
__all__ = ['Dcel', 'Vertex', 'Hedge', 'Face']

import math as m

import sys
sys.path.append('..')
from geom.point import *

class DcelError(Exception): pass

class Vertex(Point):
    """Minimal implementation of a vertex of a 2D dcel"""

    def __init__(self, px, py):
        Point.__init__(self, px, py)
        self.hedgelist = []

    def sortincident(self):
        print(self.hedgelist)
        # self.hedgelist.sort(hsort, reverse=True)
        self.hedgelist.sort(key=lambda a:a.angle, reverse=True)
        # h1.angle < h2.angle:

class Vertexx:
    """Minimal implementation of a vertex of a 2D dcel"""

    def __init__(self, px, py):
        self.x = px
        self.y = py
        self.hedgelist = []

    def sortincident(self):
        # self.hedgelist.sort(hsort, reverse=True)
        self.hedgelist.sort(keh=lambda a:a.angle, reverse=True)
    def __eq__(self, other):
        if isinstance(other, Vertex):
            return self.x==other.x and self.y==other.y
        return NotImplemented
    def __repr__(self):
        return "({0},{1})".format(self.x, self.y)

class Hedge:
    """Minimal implementation of a half-edge of a 2D dcel"""

    def __init__(self,v1,v2):
        #The origin is defined as the vertex it points to (v2)
        self.origin = v2
        self.twin = None
        self.face = None
        self.nexthedge = None
        self.angle = hangle(v2.x-v1.x, v2.y-v1.y)
        self.prevhedge = None
        self.length = m.sqrt((v2.x-v1.x)**2 + (v2.y-v1.y)**2)
    def __eq__(self, other):
        return self.origin == other.origin and \
            self.nexthedge.origin == other.nexthedge.origin
    def __repr__(self):
        if self.nexthedge is not None:
            return "({0},{1})->({2},{3})".format(self.origin.x, self.origin.y,
                                             self.nexthedge.origin.x,
                                             self.nexthedge.origin.y)
        else:
            return "({0},{1})->()".format(self.origin.x, self.origin.y)

class Face:
    """Implements a face of a 2D dcel"""

    def __init__(self):
        self.wedge = None      # outterComponent, wall edge?
        self.data = None
        self.InnerComponents = None
        self.external = None

    def area(self):
        h = self.wedge
        a = 0
        while(not h.nexthedge is self.wedge):
            p1 = h.origin
            p2 = h.nexthedge.origin
            a += p1.x*p2.y - p2.x*p1.y
            h = h.nexthedge

        p1 = h.origin
        p2 = self.wedge.origin
        a = (a + p1.x*p2.y - p2.x*p1.y)/2
        return a

    def perimeter(self):
        h = self.wedge
        p = 0
        while (not h.nexthedge is self.wedge):
            p += h.length
            h = h.nexthedge
        return p

    def vertexlist(self):
        h = self.wedge
        pl = [h.origin]
        while(not h.nexthedge is self.wedge):
            h = h.nexthedge
            pl.append(h.origin)
        return pl

    def isinside(self, p):
        """Determines whether a point is inside a face"""

        h = self.wedge
        inside = False
        if lefton(h, p):
            while(not h.nexthedge is self.wedge):
                h = h.nexthedge
                if not lefton(h, p):
                    return False
            return True
        else:
            return False

class Dcel():
    """
    Implements a doubly-connected edge list
    """

    def __init__(self, vl=[], el=[]):
        # Xygraph.__init__(self, vl, el)
        self.vl = vl
        self.el = el
        self.vertices = []
        self.hedges = []
        self.faces = []
        if vl != []:
            self.build_dcel()

    def build_dcel(self):
        """
        Creates the dcel from the list of vertices and edges
        """

#Step 1: vertex list creation
        for v in self.vl:
            self.vertices.append(Vertex(v[0], v[1]))

#Step 2: hedge list creation. Assignment of twins and
#vertices
        for e in self.el:
            if e[0] >= 0 and e[1] >= 0:
                h1 = Hedge(self.vertices[e[0]], self.vertices[e[1]])
                h2 = Hedge(self.vertices[e[1]], self.vertices[e[0]])
                h1.twin = h2
                h2.twin = h1
                self.vertices[e[1]].hedgelist.append(h1)
                self.vertices[e[0]].hedgelist.append(h2)
                self.hedges.append(h2)
                self.hedges.append(h1)

        #Step 3: Identification of next and prev hedges
        for v in self.vertices:
            v.sortincident()
            l = len(v.hedgelist)
            if l < 2:
                raise DcelError(
                    "Badly formed dcel: less than two hedges in vertex")
            else:
                """
                # The original code is incorrect. NCX 07/06/2014
                for i in range(l-1):
                    v.hedgelist[i].nexthedge = v.hedgelist[i+1].twin
                    v.hedgelist[i+1].prevhedge = v.hedgelist[i]
                v.hedgelist[l-1].nexthedge = v.hedgelist[0].twin
                v.hedgelist[0].prevhedge = v.hedgelist[l-1]
                """
                for i in range(l-1):
                    v.hedgelist[i].prevhedge = v.hedgelist[i+1].twin
                    v.hedgelist[i].twin.nexthedge = v.hedgelist[i+1]
                    v.hedgelist[i+1].prevhedge = v.hedgelist[i].twin
                    v.hedgelist[i+1].twin.nexthedge = v.hedgelist[i]
                v.hedgelist[l-1].prevhedge = v.hedgelist[0].twin
                v.hedgelist[l-1].twin.nexthedge = v.hedgelist[0]
                v.hedgelist[0].prevhedge = v.hedgelist[l-1].twin
                v.hedgelist[0].twin.nexthedge = v.hedgelist[l-1]

        #Step 4: Face assignment
        provlist = self.hedges[:]
        nf = 0
        nh = len(self.hedges)

        while nh > 0:
            h = provlist.pop()
            nh -= 1
            #We check if the hedge already points to a face
            if h.face == None:
                f = Face()
                nf += 1
                #We link the hedge to the new face
                f.wedge = h
                f.wedge.face = f
                #And we traverse the boundary of the new face
                while (not h.nexthedge is f.wedge):
                    h = h.nexthedge
                    h.face = f
                self.faces.append(f)
        #And finally we have to determine the external face
        for f in self.faces:
            f.external = f.area() < 0

    def findpoints(self, pl, onetoone=False):
        """Given a list of points pl, returns a list of
        with the corresponding face each point belongs to and
        None if it is outside the map.
        """
        ans = []
        if onetoone:
            fl = self.faces[:]
            for p in pl:
                found = False
                for f in fl:
                    if f.external:
                        continue
                    if f.isinside(p):
                        fl.remove(f)
                        found = True
                        ans.append(f)
                        break
                if not found:
                    ans.append(None)

        else:
            for p in pl:
                found = False
                for f in self.faces:
                    if f.external:
                        continue
                    if f.isinside(p):
                        found = True
                        ans.append(f)
                        break
                if not found:
                    ans.append(None)

        return ans

    def loadx(self, filename):
        """reads a dcel from file using xygraph file format"""
        a = Xygraph.load(self, filename)
        self.build_dcel()
        return a

    def load(self, pgon, edges):
        """reads a dcel from from pgon"""
        data = [[len(pgon)]]
        for i in pgon:
            data.append(i)
        for i in edges:
            data.append(i)
        nv = data[0][0]
        self.vl = data[1:nv+1]
        self.el = data[nv+1:]
        self.build_dcel()
        #return a

    def areas(self):
        return [f.area() for f in self.faces if not f.external]

    def perimeters(self):
        return [f.perimeter() for f in self.faces if not f.external]

    def nfaces(self):
        return len(self.faces)

    def nvertices(self):
        return len(self.vertices)

    def nedges(self):
        return len(self.hedges)/2

    # Functions added by NCX 2014:

    def findvertex(self, p):
        for v in self.vertices:
            if v.x == p.x and v.y == p.y:
                return v
        return None

    def findhedges(self, v1, v2): # NCX 7/4/2014
        """Finds half edges that have v1 and v2 as endpoints.
        Returns edges where v1 and v2 are the origins, respectively."""
        for h in self.hedges:
            if h.origin == v1 and h.twin.origin == v2:
                return h, h.twin
            if h.origin == v2 and h.twin.origin == v1:
                return h.twin, h
        return None, None

    def printallhedges(self):
        for e in self.hedges:
            print(e)

#Misc. functions

def hsort(h1, h2):
    """Sorts two half edges counterclockwise"""
    if h1.angle < h2.angle:
        return -1
    elif h1.angle > h2.angle:
        return 1
    else:
        return 0

def checkhedges(hl):
    """Consistency check of a hedge list: nexthedge, prevhedge"""
    for h in hl:
        if h.nexthedge not in hl or h.prevhedge not in hl:
            raise DcelError("Problems with an orphan hedge...")

def area2(hedge, point):
    """Determines the area of the triangle formed by a hedge and
    an external point"""
    pa = hedge.twin.origin
    pb=hedge.origin
    pc=point
    return (pb.x - pa.x)*(pc[1] - pa.y) - (pc[0] - pa.x)*(pb.y - pa.y)

def lefton(hedge, point):
    """Determines if a point is to the left of a hedge"""
    return area2(hedge, point) >= 0

def hangle(dx,dy):
    """Determines the angle with respect to the x axis of a segment
    of coordinates dx and dy
    """
    l = m.sqrt(dx*dx + dy*dy)
    if dy > 0:
        return m.acos(dx/l)
    else:
        return 2*m.pi - m.acos(dx/l)

if __name__=='__main__':
    d = Dcel()
    pgon = [ [0,10], [5,0], [10,10], [15,0], [20,10], [25,0], [30,20],
             [40,20], [45,0], [50,50], [40,40], [30,50], [25,20],
             [20,50], [15,10], [10,50], [8, 8], [4,50] ]
    edges = [ [0,1], [1,2], [2,3], [3,4], [4,5],[5,6],[6,7],[7,8],
              [8,9],[9,10],[10,11],[11,12],[12,13],
              [13,14],[14,15],[15,16],[16,17],[17,0] ]
    d.load(pgon, edges)
    for a,p in zip(d.areas(), d.perimeters()):
        print(a,p)
