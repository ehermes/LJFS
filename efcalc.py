#! /usr/bin/env python

import sys
import copy
import operator
import numpy as np
import readarc
import readkey
from analyze2 import analyze2 as analyze

def e_f(nsolatoms,optlist,epsilon,sigma,charge,n,x,y,z,atomtype,conv):
    elist = np.zeros(n)
    flist = np.zeros((n,len(optlist),3))
    for i in xrange(n):
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
                    elist[i] += analyze.e_lj(ep1,ep2,sig1,sig2,r)
                    if optatomnum != -1:
                        flist[i][optatomnum] += analyze.f_lj(ep1,ep2,sig1,sig2,r)
                if q1 != 0 and q2 != 0:
                    elist[i] += analyze.ecoul(conv,q1,q2,r)
                    if optatomnum != -1:
                        flist[i][optatomnum] += analyze.fcoul(conv,q1,q2,r)
    return elist, flist

def de_df(nsolatoms,optatoms,optlist,epsilon,sigma,charge,n,x,y,z,atomtype,conv):
    delist = np.zeros((n,3*len(optatoms)))
    dflist = np.zeros((n,3*len(optatoms),len(optlist),3))
    for i in xrange(n):
        for j in xrange(len(optlist)):
            h = optlist[j]
            g = -1
            for u in xrange(len(optatoms)):
                if atomtype[i][h] == optatoms[u]:
                    g = u
                    break
            else:
                print "Something terrible has happened."
                sys.exit()
            for k in xrange(nsolatoms,len(x[i])):
                ep1 = epsilon[atomtype[i][h]]
                ep2 = epsilon[atomtype[i][k]]
                sig1 = sigma[atomtype[i][h]]
                sig2 = sigma[atomtype[i][k]]
                r = np.array([x[i][h]-x[i][k],y[i][h]-y[i][k],z[i][h]-z[i][k]])
                q1 = charge[atomtype[i][h]]
                q2 = charge[atomtype[i][k]]
                if ep1 != 0 and ep2 != 0:
                    (delist[i][3*g], delist[i][3*g+1], dflist[i][3*g][j], 
                            dflist[i][3*g+1][j]) = map(operator.add, (delist[i][3*g],
                                delist[i][3*g+1], dflist[i][3*g][j], dflist[i][3*g+1][j]),
                                analyze.dedf_lj(ep1,ep2,sig1,sig2,r))
                if q1 != 0 and q2 != 0:
                    (delist[i][3*g+2], dflist[i][3*g+2][j]) = map(operator.add,
                            (delist[i][3*g+2], dflist[i][3*g+2][j]), 
                            analyze.dedf_coul(conv,q1,q2,r))
    return delist, dflist

def de_df_chcomp(nsolatoms,chcomptype,epsilon,sigma,charge,n,x,y,z,atomtype,conv):
    delist = []
    dflist = []
    for i in xrange(n):
        de = 0
        df = []
        for j in xrange(nsolatoms):
            if atomtype[i][j] == chcomptype:
                dfi = np.zeros(3)
                for k in xrange(nsolatoms,len(x[i])):
                    r = np.array([x[i][j]-x[i][k],y[i][j]-y[i][k],z[i][j]-z[i][k]])
                    q1 = charge[atomtype[i][j]]
                    q2 = charge[atomtype[i][k]]
                    if q1 != 0 and q2 != 0:
                        (de, dfi) = map(operator.add, (de, dfi), 
                                analyze.dedf_coul(conv,q1,q2,r))
                df.append(dfi)
            else:
                df.append(np.zeros(3))
        delist.append(de)
        dflist.append(df)
    return np.array(delist), np.array(dflist)

def dipquad(charge,n,nsolatoms,x,y,z,atomtype):
    debye = 4.803205
    diplist = np.zeros((n,3))
    quadlist = np.zeros((n,6))
    for i in xrange(n):
        for j in xrange(nsolatoms):
            diplist[i] += charge[atomtype[i][j]]*np.array([x[i][j],y[i][j],z[i][j]])
            rsqi = x[i][j]**2 + y[i][j]**2 + z[i][j]**2
            quadlist[i][0] += charge[atomtype[i][j]]*(3*x[i][j]**2 - rsqi)
            quadlist[i][1] += charge[atomtype[i][j]]*(3*y[i][j]**2 - rsqi)
            quadlist[i][2] += charge[atomtype[i][j]]*(3*z[i][j]**2 - rsqi)
            quadlist[i][3] += charge[atomtype[i][j]]*(3*x[i][j]*y[i][j])
            quadlist[i][4] += charge[atomtype[i][j]]*(3*x[i][j]*z[i][j])
            quadlist[i][5] += charge[atomtype[i][j]]*(3*y[i][j]*z[i][j])
    diplist *= debye
    quadlist *= debye
    return diplist, quadlist

def ddipdquad(optatoms,optlist,n,x,y,z,atomtype):
    debye = 4.803205
    ddiplist = np.zeros((n,3*len(optatoms),3))
    dquadlist = np.zeros((n,3*len(optatoms),6))
    for i in xrange(n):
        for j in xrange(len(optlist)):
            h = optlist[j]
            g = -1
            for u in xrange(len(optatoms)):
                if atomtype[i][h] == optatoms[u]:
                    g = u
                    break
            else:
                print "Something terrible has happened."
                sys.exit()
            rsqi = x[i][h]**2 + y[i][h]**2 + z[i][h]**2
            ddiplist[i][3*g+2] += np.array([x[i][h],y[i][h],z[i][h]])
            dquadlist[i][3*g+2][0] += 3*x[i][h]**2 - rsqi
            dquadlist[i][3*g+2][1] += 3*y[i][h]**2 - rsqi
            dquadlist[i][3*g+2][2] += 3*z[i][h]**2 - rsqi
            dquadlist[i][3*g+2][3] += 3*x[i][h]*y[i][h]
            dquadlist[i][3*g+2][4] += 3*x[i][h]*z[i][h]
            dquadlist[i][3*g+2][5] += 3*y[i][h]*z[i][h]
    ddiplist *= debye
    dquadlist *= debye
    return ddiplist, dquadlist

def ddip_dquad_chcomp(nsolatoms,chcomptype,charge,n,x,y,z,atomtype):
    debye = 4.803205
    ddiplist = np.zeros((n,3))
    dquadlist = np.zeros((n,6))
    for i in xrange(n):
        for j in xrange(nsolatoms):
            if atomtype[i][j] == chcomptype:
                rsqi = x[i][j]**2 + y[i][j]**2 + z[i][j]**2
                ddiplist[i] += np.array([x[i][j],y[i][j],z[i][j]])
                dquadlist[i][0] += 3*x[i][j]**2 - rsqi
                dquadlist[i][1] += 3*y[i][j]**2 - rsqi
                dquadlist[i][2] += 3*z[i][j]**2 - rsqi
                dquadlist[i][3] += 3*x[i][j]*y[i][j]
                dquadlist[i][4] += 3*x[i][j]*z[i][j]
                dquadlist[i][5] += 3*y[i][j]*z[i][j]
    ddiplist *= debye
    dquadlist *= debye
    return ddiplist, dquadlist

