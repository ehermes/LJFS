#! /usr/bin/env python

import sys
import copy
import numpy as np
import readarc
import readkey
import readgauss
import analyze
import efcalc
from scipy import optimize

fileprefix = sys.argv[1]
nsolatoms = int(sys.argv[2])
optatomname = sys.argv[3]

atom_tinker = {900:'C', 901:'C', 902:'C', 910:'H', 911:'H', 912:'H', 920:'Cl', 
        921:'Cl', 922:'Cl', 55:'O', 56:'H', 57:''}
atom_gauss = {1:'H', 6:'C', 8:'O', 17:'Cl'}

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

optatoms = []
for i in xrange(nsolatoms):
    if atom_tinker[atomtype[0][i]] == optatomname:
        optatoms.append(i)

if optatoms == []:
    print "No atoms to optimize parameters for!"
    sys.exit()

[e_solute, f_solute] = readgauss.readsolute(fileprefix)
[e_cluster, f_cluster] = readgauss.readcluster(fileprefix,n)
e_water = readgauss.readwater(fileprefix,n)

e_qm = []
f_qm = []
for i in xrange(n):
    e_tmp = e_cluster[i] - e_water[i] - e_solute
    e_tmp *= 627.509
    e_qm.append(e_tmp)
    f_tmp = []
    for j in xrange(len(f_cluster[i])):
        if atom_gauss[f_cluster[i][j][0]] == optatomname:
            fi_tmp = np.array([f_cluster[i][j][1], f_cluster[i][j][2], 
                f_cluster[i][j][3]]) 
            fi_tmp -= np.array([f_solute[j][1], f_solute[j][2], f_solute[j][3]])
            fi_tmp *= 1185.821
            f_tmp.append(fi_tmp)
    f_qm.append(f_tmp)

e_std = np.std(e_qm)
f_std = np.std(f_qm,axis=0)

optnum = 921

def chisq(prms):
    epsilon = prm_ep
    sigma = prm_sigma
    charge = prm_ch
    epsilon[optnum] = prms[0]
#    sigma[optnum] = 4.05562627
    sigma[optnum] = prms[1]
    [energy, force] = efcalc.e_f(optatoms,epsilon,sigma,charge,n,x,y,z,atomtype,
            newbonds,conv)
    chi_e = 0.0
    chi_f = np.zeros((len(optatoms),3))
    for i in xrange(n):
        chi_e += (energy[i] - e_qm[i])**2
        for j in xrange(len(optatoms)):
            chi_f += (force[i][j] - f_qm[i][j])**2
    chi_e /= e_std**2
    chi_f /= f_std**2
    chi_sq = 0.5 * chi_e + 0.5 * (np.sum(chi_f)/(3*len(optatoms)))
    chi = np.sqrt(chi_sq)
    return chi

def gradchisq(prms):
    epsilon = prm_ep
    sigma = prm_sigma
    charge = prm_ch
    epsilon[optnum] = prms[0]
    sigma[optnum] = prms[1]
    [energy, force] = efcalc.e_f(optatoms,epsilon,sigma,charge,n,x,y,z,atomtype,
            newbonds,conv)
    [denergy, dforce] = efcalc.de_df(optatoms,epsilon,sigma,charge,n,x,y,z,
            atomtype,newbonds,conv)
    dchi_e = np.zeros(2)
    dchi_f = np.zeros((len(optatoms),2,3))
    for i in xrange(n):
        dchi_e += (energy[i] - e_qm[i]) * denergy[i][:2]
        for j in xrange(len(optatoms)):
            dchi_f[j] += (force[i][j] - f_qm[i][j]) * dforce[i][j][:2]
    dchi_e *= 2/(e_std**2)
    dchi_f *= 2/(f_std**2)
    dchi_sq = 0.5 * dchi_e + 0.5 * (np.sum(np.sum(dchi_f,axis=0),axis=1)/(3*len(optatoms)))
    return dchi_sq

#chi = chisq(0.11800,4.41720)
#dchi = gradchisq(0.11800,4.41720)

#print chi
#print dchi

#initprm = np.array([0.10905819,  4.05590243])
initprm = np.array([0.11099003,  4.06886509])
#initprm = np.array([0.103,4.1])
#initprm = 0.108
prmbounds=[(0.09,0.12),(3.5,5.0)]
lower = np.array([0.09,3.5])
upper = np.array([0.12,5.0])
bruterange = ((0.09,0.12),(3.5,5.0))

#optvalues = optimize.fmin(chisq,initprm)
#optvalues = optimize.anneal(chisq,initprm,lower=lower,upper=upper,schedule='boltzmann')
#optvalues = optimize.fmin_tnc(chisq,initprm,gradchisq,bounds=prmbounds)
optvalues = optimize.fmin_bfgs(chisq,initprm)

print optvalues
