#! /usr/bin/env python

import os
import sys
import math
import re
import readarc

fileprefix = sys.argv[1]

arcname = fileprefix + '.arc'
ljoname = fileprefix + '.ljo'

arc = readarc.readarc(arcname)

ljofile = open(ljoname, 'r')
nsolatoms = int(re.split('\s+',ljofile.readline().strip())[0])
charge = re.split('\s+',ljofile.readline().strip())[0]
ljofile.close()

n = arc[0]
x = arc[2]
y = arc[3]
z = arc[4]
atomtype = arc[5]

symbols = {900:'C', 901:'C', 902:'C', 910:'H', 911:'H', 912:'H', 920:'Cl', 
        921:'Cl', 922:'Cl', 930:'O', 931:'O', 932:'O', 940:'H', 941:'H', 942:'H', 
        55:'O', 56:'H', 57:''}
chgmult = str(int(float(charge))) + ' 1 ' + str(int(float(charge))) + ' 1 0 1'

dirname = 'gauss_' + fileprefix

if not os.path.exists(dirname):
    os.makedirs(dirname)

#footer = '\n@/home/ehermes/local/basis/def2-svpd.gbs\n\n'
footer = '\n@def2-svpd.gbs\n\n'

for i in xrange(n):
    solutename = fileprefix + '_' + str(i)
    clustername = fileprefix + '_cluster_' + str(i)
    watername = fileprefix + '_water_' + str(i)
    solutefile = open(dirname + '/' + solutename + '.com', 'w')
    clusterfile = open(dirname + '/' + clustername + '.com','w')
    waterfile = open(dirname + '/' + watername + '.com','w')
    headersolute = ('# wb97xd/gen freq=hpmodes IOp(2/15=1)\n\nconfiguration ' +str(i) + 
            ' with only solute\n\n' + charge.split('.')[0] + ' 1\n')
    headercluster = ('# wb97xd/gen counterpoise=2 force IOp(2/15=1) pop=chelpg' + 
            '\n\nconfiguration ' + str(i) + ' with solute and water\n\n' + chgmult + '\n')
    headerwater = ('# wb97xd/gen IOp(2/15=1)\n\nconfiguration ' + str(i) + 
            ' water only\n\n0 1\n')
    solutefile.write(headersolute)
    clusterfile.write(headercluster)
    waterfile.write(headerwater)
    for j in xrange(nsolatoms):
        atomsymbol = symbols[atomtype[i][j]]
        soluteline = '{0:15} {1:.6f} {2:.6f} {3:.6f}\n'.format(atomsymbol, x[i][j],
                y[i][j], z[i][j])
        clusterline = '{0:15} {1:.6f} {2:.6f} {3:.6f}\n'.format(atomsymbol +
                '(Fragment=1)', x[i][j], y[i][j], z[i][j])
        solutefile.write(soluteline)
        clusterfile.write(clusterline)
    for j in xrange(nsolatoms,len(x[i])):
        atomsymbol = symbols[atomtype[i][j]]
        if atomsymbol != '':
            clusterline = '{0:15} {1:.6f} {2:.6f} {3:.6f}\n'.format(atomsymbol + 
                    '(Fragment=2)', x[i][j], y[i][j], z[i][j])
            waterline = '{0:15} {1:.6f} {2:.6f} {3:.6f}\n'.format(atomsymbol, 
                    x[i][j], y[i][j], z[i][j])
            clusterfile.write(clusterline)
            waterfile.write(waterline)
    solutefile.write(footer)
    clusterfile.write(footer)
    waterfile.write(footer)
    solutefile.close()
    clusterfile.close()
    waterfile.close()

