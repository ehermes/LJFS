#! /usr/bin/env python

import sys
import re
import numpy as np
import readarc

fileprefix = sys.argv[1]
cutoff = float(sys.argv[2])

filename = fileprefix + '.arc'
ljoname = fileprefix + '.ljo'

arc = readarc.readarc(filename)

n = arc[0]
atomname = arc[1]
x = arc[2]
y = arc[3]
z = arc[4]
atomtype = arc[5]
bonds = arc[6]

ljofile = open(ljoname, 'r')
nsolatoms = int(re.split('\s+',ljofile.readline().strip())[0])
ljofile.close()

#cutoff = 4.0

newatomname = []
newx = []
newy = []
newz = []
newatomtype = []
newbonds = []

for i in xrange(n):
    newatomnum = 0
    newatomname.append([])
    newx.append([])
    newy.append([])
    newz.append([])
    newatomtype.append([])
    newbonds.append([])
    newatomname[i].append(atomname[i][0])
    newx[i].append(x[i][0])
    newy[i].append(y[i][0])
    newz[i].append(z[i][0])
    newatomtype[i].append(atomtype[i][0])
    newbonds[i].append(bonds[i][0])
    for j in xrange(nsolatoms,len(x[i]),4):
        rmag = np.sqrt((x[i][0] - x[i][j])**2 + (y[i][0] - y[i][j])**2 + 
                (z[i][0] - z[i][j])**2)
        if rmag < cutoff:
            newatomnum += 1
            for k in xrange(j, j+4):
                newatomname[i].append(atomname[i][k])
                newx[i].append(x[i][k])
                newy[i].append(y[i][k])
                newz[i].append(z[i][k])
                newatomtype[i].append(atomtype[i][k])
            newbonds[i].append([4*newatomnum-1,4*newatomnum,4*newatomnum+1])
            newbonds[i] += 3*[[4*newatomnum-2]]

out = open(fileprefix + '_cluster.arc','w')

for i in xrange(n):
    numatoms = len(newx[i])
    header = '{0:4d} {1:10} cluster\n'.format(numatoms, fileprefix)
    out.write(header)
    for j in xrange(numatoms):
        atomline = '{0:4d} {1:3} {2:.6f} {3:.6f} {4:.6f} {5:4d}'.format(j+1, 
                newatomname[i][j], newx[i][j], newy[i][j], newz[i][j], newatomtype[i][j])
        for k in xrange(len(newbonds[i][j])):
            atomline += ' ' + str(newbonds[i][j][k])
        atomline += '\n'
        out.write(atomline)

out.close()
