#! /usr/bin/env python

import os
import sys
import math
import readarc

fileprefix = sys.argv[1]

filename = fileprefix + '.arc','r'

arc = readarc.readarc(filename)

n = arc[0]
x = arc[2]
y = arc[3]
z = arc[4]
atomtype = arc[5]

lastsolute = 1
symbols = {'900':'C', '901':'C', '902':'C', '910':'H', '911':'H', '912':'H', '920':'Cl', '921':'Cl', '922':'Cl', '55':'O', '56':'H', '57':''}
chgmult = '-1 1 -1 1 0 1'

dirname = 'gauss_' + fileprefix

if not os.path.exists(dirname):
    os.makedirs(dirname)

footer = '\n@def2-svpd.gbs\n\n'

for i in xrange(n):
    clustername = fileprefix + '_cluster_' + str(i)
    watername = fileprefix + '_water_' + str(i)
    clusterfile = open(dirname + '/' + clustername + '.com','w')
    waterfile = open(dirname + '/' + watername + '.com','w')
    headercluster = '%chk=' + clustername + '.chk\n# wb97xd/gen counterpoise=2\n\nconfiguration ' + str(i) + ' with solute and water\n\n' + chgmult + '\n'
    headerwater =  '%chk=' + watername + '.chk\n# wb97xd/gen\n\nconfiguration ' + str(i) + ' water only\n\n0 1\n'
    clusterfile.write(headercluster)
    waterfile.write(headerwater)
    for j in xrange(lastsolute):
        atomsymbol = symbols[atomtype[i][j]]
        clusterline = '{0:15} {1:.6f} {2:.6f} {3:.6f}\n'.format(atomsymbol +
                '(Fragment=1)', x[i][j], y[i][j], z[i][j])
        clusterfile.write(clusterline)
    for j in xrange(lastsolute,len(x[i])):
        atomsymbol = symbols[atomtype[i][j]]
        if atomsymbol != '':
            clusterline = '{0:15} {1:.6f} {2:.6f} {3:.6f}\n'.format(atomsymbol + 
                    '(Fragment=2)', x[i][j], y[i][j], z[i][j])
            waterline = '{0:15} {1:.6f} {2:.6f} {3:.6f}\n'.format(atomsymbol, 
                    x[i][j], y[i][j], z[i][j])
            clusterfile.write(clusterline)
            waterfile.write(waterline)
    clusterfile.write(footer)
    waterfile.write(footer)
    clusterfile.close()
    waterfile.close()

