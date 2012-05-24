#! /usr/bin/env python

import sys
import copy
import numpy as np
import readarc
import readkey
import analyze

fileprefix = sys.argv[1]

arcname = fileprefix + '.arc'
keyname = fileprefix + '.key'

arcdata = readarc.readarc(arcname)

n = arcdata[0]
x = arcdata[2]
y = arcdata[3]
z = arcdata[4]
atomtype = arcdata[5]
bonds = arcdata[6]
newbonds = copy.deepcopy(bonds)

for i in xrange(n):
    for j in xrange(len(x[i])):
        for k in bonds[i][j]:
            newbonds[i][k-1] += bonds[i][j]

keydata = readkey.readkey(keyname)

prm_sigma = keydata[0]
prm_ep = keydata[1]
prm_ch = keydata[2]
conv = keydata[3]

def eandf(atomnum,n,x,y,z,atomtypes,prm_sigma,prm_ep,prm_ch,conv):
    elist = []
    delist = []
    flist = []
    dflist = []
    for i in xrange(n):
    energy = 0.0
    gradenergy = 0.0
    force = np.array([0.0,0.0,0.0])
    gradforce = np.array([0.0,0.0,0.0])
        for j in xrange(len(x[i])):
            if j != atomnum and j+1 not in newbonds[i][j]:
                ep1 = prm_ep[atomtype[i][atomnum]]
                ep2 = prm_ep[atomtype[i][j]]
                sig1 = prm_sigma[atomtype[i][atomnum]]
                sig2 = prm_sigma[atomtype[i][j]]
                r = np.array([x[i][atomnum]-x[i][j],y[i][atomnum]-y[i][j],
                    z[i][atomnum]-z[i][j]])
                q1 = prm_ch[atomtype[i][atomnum]]
                q2 = prm_ch[atomtype[i][j]]
                if ep1 != 0 and ep2 != 0:
                    energy += analyze.elj12(ep1,ep2,sig1,sig2,r)
                    energy += analyze.elj12(ep1,ep2,sig1,sig2,r)
                    gradenergy += analyze.ddepelj12(ep1,ep2,sig1,sig2,r)
                    gradenergy += analyze.ddsigelj12(ep1,ep2,sig1,sig2,r)
                    gradenergy += analyze.ddepelj6(ep1,ep2,sig1,sig2,r)
                    gradenergy += analyze.ddsigelj6(ep1,ep2,sig1,sig2,r)
                    force += analyze.flj12(ep1,ep2,sig1,sig2,r)
                    force += analyze.flj6(ep1,ep2,sig1,sig2,r)
                    gradforce += analyze.ddepflj12(ep1,ep2,sig1,sig2,r)
                    gradforce += analyze.ddsigflj12(ep1,ep2,sig1,sig2,r)
                    gradforce += analyze.ddepflj6(ep1,ep2,sig1,sig2,r)
                    gradforce += analyze.ddsigflj6(ep1,ep2,sig1,sig2,r)
                if q1 != 0 and q2 != 0:
                    energy += analyze.ecoul(conv,q1,q2,r)
                    gradenergy += analyze.ddqecoul(conv,q1,q2,r)
                    force += analyze.fcoul(conv,q1,q2,r)
                    gradforce += ddfcoul(conv,q1,q2,r)
        elist.append(energy)
        delist.append(gradenergy)
        flist.append(force)
        dflist.append(gradforce)

    return [elist,delist,flist,dflist]
