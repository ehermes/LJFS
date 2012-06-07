#! /usr/bin/env python

import sys
import copy
import numpy as np
import readarc
import readkey
import analyze

def e_f(optatoms, epsilon, sigma, charge, n, x, y, z, atomtype, newbonds, conv):
    elist = []
    flist = []
    force_tmp = []
    for i in xrange(n):
        energy = 0.0
        force = np.zeros((len(optatoms),3))
        for j in xrange(len(x[i])):
            for k in xrange(len(optatoms)):
                if j not in optatoms and j+1 not in newbonds[i][optatoms[k]]:
#            if j not in optatoms:
#                for k in xrange(len(optatoms)):
                    ep1 = epsilon[atomtype[i][optatoms[k]]]
                    ep2 = epsilon[atomtype[i][j]]
                    sig1 = sigma[atomtype[i][optatoms[k]]]
                    sig2 = sigma[atomtype[i][j]]
                    r = np.array([x[i][optatoms[k]]-x[i][j],y[i][optatoms[k]]-y[i][j],
                        z[i][optatoms[k]]-z[i][j]])
                    q1 = charge[atomtype[i][optatoms[k]]]
                    q2 = charge[atomtype[i][j]]
                    if ep1 != 0 and ep2 != 0:
                        energy += analyze.elj12(ep1,ep2,sig1,sig2,r)
                        energy += analyze.elj6(ep1,ep2,sig1,sig2,r)
                        force[k] += analyze.flj12(ep1,ep2,sig1,sig2,r)
                        force[k] += analyze.flj6(ep1,ep2,sig1,sig2,r)
                    if q1 != 0 and q2 != 0:
                        energy += analyze.ecoul(conv,q1,q2,r)
                        force[k] += analyze.fcoul(conv,q1,q2,r)
        elist.append(energy)
        flist.append(force)
    return [elist,flist]

def de_df(optatoms, epsilon, sigma, charge, n, x, y, z, atomtype, newbonds, conv):
    delist = []
    dflist = []
    for i in xrange(n):
        de = np.zeros(3)
        df = np.zeros((len(optatoms),3,3))
        for j in xrange(len(x[i])):
            for k in xrange(len(optatoms)):
                if j not in optatoms and j+1 not in newbonds[i][optatoms[k]]:
#            if j not in optatoms:
#                for k in xrange(len(optatoms)):
                    ep1 = epsilon[atomtype[i][optatoms[k]]]
                    ep2 = epsilon[atomtype[i][j]]
                    sig1 = sigma[atomtype[i][optatoms[k]]]
                    sig2 = sigma[atomtype[i][j]]
                    r = np.array([x[i][optatoms[k]]-x[i][j],y[i][optatoms[k]]-y[i][j],
                        z[i][optatoms[k]]-z[i][j]])
                    q1 = charge[atomtype[i][optatoms[k]]]
                    q2 = charge[atomtype[i][j]]
                    if ep1 != 0 and ep2 != 0:
                        de[0] += analyze.ddepelj12(ep1,ep2,sig1,sig2,r) + \
                                analyze.ddepelj6(ep1,ep2,sig1,sig2,r)
                        de[1] += analyze.ddsigelj12(ep1,ep2,sig1,sig2,r) + \
                                analyze.ddsigelj6(ep1,ep2,sig1,sig2,r)
                        df[k][0] += analyze.ddepflj12(ep1,ep2,sig1,sig2,r) + \
                                analyze.ddepflj6(ep1,ep2,sig1,sig2,r)
                        df[k][1] += analyze.ddsigflj12(ep1,ep2,sig1,sig2,r) + \
                                analyze.ddsigflj6(ep1,ep2,sig1,sig2,r)
                    if q1 != 0 and q2 != 0:
                        de[2] += analyze.ddqecoul(conv,q1,q2,r)
                        df[k][2] += analyze.ddqfcoul(conv,q1,q2,r)
        delist.append(de)
        dflist.append(df)
    return [delist,dflist]
