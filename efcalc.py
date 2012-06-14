#! /usr/bin/env python

import sys
import copy
import numpy as np
import readarc
import readkey
import analyze

def e_f(nsolatoms,optlist,epsilon,sigma,charge,n,x,y,z,atomtype,newbonds,conv):
    elist = []
    flist = []
    for i in xrange(n):
        energy = 0.0
        force = np.zeros((len(optlist),3))
        for j in xrange(nsolatoms):
            for l in xrange(len(optlist)):
                if j == optlist[l]:
                    optatomnum = l
                    break
            else:
                optatomnum = -1
            for k in xrange(nsolatoms,len(x[i])):
                ep1 = epsilon[atomtype[i][j]]
                ep2 = epsilon[atomtype[i][k]]
                sig1 = sigma[atomtype[i][j]]
                sig2 = sigma[atomtype[i][k]]
                r = np.array([x[i][j]-x[i][k],y[i][j]-y[i][k],z[i][j]-z[i][k]])
                q1 = charge[atomtype[i][j]]
                q2 = charge[atomtype[i][k]]
                if ep1 != 0 and ep2 != 0:
                    energy += analyze.elj12(ep1,ep2,sig1,sig2,r)
                    energy += analyze.elj6(ep1,ep2,sig1,sig2,r)
                    if optatomnum != -1:
                        force[optatomnum] += analyze.flj12(ep1,ep2,sig1,sig2,r)
                        force[optatomnum] += analyze.flj6(ep1,ep2,sig1,sig2,r)
                if q1 != 0 and q2 != 0:
                    energy += analyze.ecoul(conv,q1,q2,r)
                    if optatomnum:
                        force[optatomnum] += analyze.fcoul(conv,q1,q2,r)
        elist.append(energy)
        flist.append(force)
    return [elist,flist]

def de_df(nsolatoms,optatoms,optlist,epsilon,sigma,charge,n,x,y,z,atomtype,newbonds,conv):
    delist = []
    dflist = []
    for i in xrange(n):
        de = np.zeros(3*len(optatoms))
        df = np.zeros((len(optlist),3*len(optatoms),3))
        for j in xrange(len(optlist)):
            l = optlist[j]
            for k in xrange(nsolatoms,len(x[i])):
                ep1 = epsilon[atomtype[i][l]]
                ep2 = epsilon[atomtype[i][k]]
                sig1 = sigma[atomtype[i][l]]
                sig2 = sigma[atomtype[i][k]]
                r = np.array([x[i][l]-x[i][k],y[i][l]-y[i][k],z[i][l]-z[i][k]])
                q1 = charge[atomtype[i][l]]
                q2 = charge[atomtype[i][k]]
                if ep1 != 0 and ep2 != 0:
                    de[3*j] += analyze.ddepelj12(ep1,ep2,sig1,sig2,r) + \
                            analyze.ddepelj6(ep1,ep2,sig1,sig2,r)
                    de[3*j+1] += analyze.ddsigelj12(ep1,ep2,sig1,sig2,r) + \
                            analyze.ddsigelj6(ep1,ep2,sig1,sig2,r)
                    df[j][3*j] += analyze.ddepflj12(ep1,ep2,sig1,sig2,r) + \
                            analyze.ddepflj6(ep1,ep2,sig1,sig2,r)
                    df[j][3*j+1] += analyze.ddsigflj12(ep1,ep2,sig1,sig2,r) + \
                            analyze.ddsigflj6(ep1,ep2,sig1,sig2,r)
                if q1 != 0 and q2 != 0:
                    de[3*j+2] += analyze.ddqecoul(conv,q1,q2,r)
                    df[j][3*j+2] += analyze.ddqfcoul(conv,q1,q2,r)
        delist.append(de)
        dflist.append(df)
    return [delist,dflist]
