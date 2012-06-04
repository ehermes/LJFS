#! /usr/bin/env python

import sys
import copy
import numpy as np
import readarc
import readkey
import analyze

fileprefix = sys.argv[1]
optatom = int(sys.argv[2])

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

def chisq(e,s):
    epsilon = prm_ep
    sigma = prm_sigma
    charge = prm_ch
    epsilon[atomtype[optatom - 1]] = e
    sigma[atomtype[optatom - 1]] = s
    [energy,force] = e_f(epsilon,sigma,charge)

def gradchisq(e,s):
    epsilon = prm_ep
    sigma = prm_sigma
    charge = prm_ch
    epsilon[atomtype[optatom - 1]] = e
    sigma[atomtype[optatom - 1]] = s
    [denergy,dforce] = de_df(epsilon,sigma,charge)

def e_f(epsilon=prm_ep,sigma=prm_sigma,charge=prm_ch):
    elist = []
    flist = []
    for i in xrange(n):
    energy = 0.0
    force = np.array([0.0,0.0,0.0])
        for j in xrange(len(x[i])):
            if j != atomnum and j+1 not in newbonds[i][j]:
                ep1 = epsilon[atomtype[i][atomnum]]
                ep2 = epsilon[atomtype[i][j]]
                sig1 = sigma[atomtype[i][atomnum]]
                sig2 = sigma[atomtype[i][j]]
                r = np.array([x[i][atomnum]-x[i][j],y[i][atomnum]-y[i][j],
                    z[i][atomnum]-z[i][j]])
                q1 = charge[atomtype[i][atomnum]]
                q2 = charge[atomtype[i][j]]
                if ep1 != 0 and ep2 != 0:
                    energy += analyze.elj12(ep1,ep2,sig1,sig2,r)
                    energy += analyze.elj12(ep1,ep2,sig1,sig2,r)
                    force += analyze.flj12(ep1,ep2,sig1,sig2,r)
                    force += analyze.flj6(ep1,ep2,sig1,sig2,r)
                if q1 != 0 and q2 != 0:
                    energy += analyze.ecoul(conv,q1,q2,r)
                    force += analyze.fcoul(conv,q1,q2,r)
        elist.append(energy)
        flist.append(force)
    return [elist,flist]

def de_df(epsilon=prm_ep,sigma=prm_sigma,charge=prm_ch):
    delist = []
    dflist = []
    for i in xrange(n):
    dedep = 0.0
    dedsig = 0.0
    dedq = 0.0
    dfdep = np.array([0.0,0.0,0.0])
    dfdsig = np.array([0.0,0.0,0.0])
    dfdq = np.array([0.0,0.0,0.0])
        for j in xrange(len(x[i])):
            if j != atomnum and j+1 not in newbonds[i][j]:
                ep1 = epsilon[atomtype[i][atomnum]]
                ep2 = epsilon[atomtype[i][j]]
                sig1 = sigma[atomtype[i][atomnum]]
                sig2 = sigma[atomtype[i][j]]
                r = np.array([x[i][atomnum]-x[i][j],y[i][atomnum]-y[i][j],
                    z[i][atomnum]-z[i][j]])
                q1 = charge[atomtype[i][atomnum]]
                q2 = charge[atomtype[i][j]]
                if ep1 != 0 and ep2 != 0:
                    dedep += analyze.ddepelj12(ep1,ep2,sig1,sig2,r) + \
                            analyze.ddepelj6(ep1,ep2,sig1,sig2,r)
                    dedsig += analyze.ddsigelj12(ep1,ep2,sig1,sig2,r) + \
                            analyze.ddsigelj6(ep1,ep2,sig1,sig2,r)
                    dfdep += analyze.ddepflj12(ep1,ep2,sig1,sig2,r) + \
                            analyze.ddepflj6(ep1,ep2,sig1,sig2,r)
                    dfdsig += analyze.ddsigflj12(ep1,ep2,sig1,sig2,r) + \
                            analyze.ddsigflj6(ep1,ep2,sig1,sig2,r)
                if q1 != 0 and q2 != 0:
                    dedq += analyze.ddqecoul(conv,q1,q2,r)
                    dfdq += ddfcoul(conv,q1,q2,r)
        delist.append(np.array([dedep,dedsig,dedq]))
        dflist.append(np.array([dfdep,dfdsig,dfdq]))
    return [delist,dflist]
