#! /usr/bin/env python

import sys
import copy
import re
import numpy as np
import readarc
import readkey
import readgauss
import analyze
import efcalc
from scipy import optimize

kb = 1.9872041e-3
T = 298
weight_e = 1.0
weight_f = 1.0

fileprefix = sys.argv[1]

atom_tinker = {901:'C', 911:'H', 921:'Cl', 55:'O', 56:'H', 57:''}
atom_gauss = {1:'H', 6:'C', 8:'O', 17:'Cl'}

arcname = fileprefix + '.arc'
keyname = fileprefix + '.key'
ljoname = fileprefix + '.ljo'
outname = fileprefix + '.out'

arcdata = readarc.readarc(arcname)
n = arcdata[0]
x = arcdata[2]
y = arcdata[3]
z = arcdata[4]
atomtype = arcdata[5]
bonds = arcdata[6]

keydata = readkey.readkey(keyname)
prm_sigma = keydata[0]
prm_ep = keydata[1]
prm_ch = keydata[2]
conv = keydata[3]

ljofile = open(ljoname, 'r')

nsolatoms = int(re.split('\s+',ljofile.readline().strip())[0])
chargedata = re.split('\s+',ljofile.readline().strip())
istranstate = bool(int(re.split('\+',ljofile.readline().strip())[0]))
syscharge = float(chargedata[0])
chcomptype = int(chargedata[1])

optatoms = []
initprms = []
prmbounds = []

compcharge = False

while True:
    line = ljofile.readline().strip()
    if line == '': break
    data = re.split('\s+', line)
    optatoms.append(int(data[0]))
    initprms.append(prm_ep[int(data[0])])
    if int(data[1]):
        prmbounds.append((prm_ep[int(data[0])]*0.000001,prm_ep[int(data[0])]*1000.0))
    else:
        prmbounds.append((prm_ep[int(data[0])],prm_ep[int(data[0])]))
    initprms.append(prm_sigma[int(data[0])])
    if int(data[2]):
        prmbounds.append((prm_sigma[int(data[0])]*0.001,prm_sigma[int(data[0])]*1000.0))
    else:
        prmbounds.append((prm_sigma[int(data[0])],prm_sigma[int(data[0])]))
    initprms.append(prm_ch[int(data[0])])
    if int(data[3]):
        compcharge = True
        prmbounds.append((prm_ch[int(data[0])]-1.0,prm_ch[int(data[0])]+1.0))
    else:
        prmbounds.append((prm_ch[int(data[0])],prm_ch[int(data[0])]))

ljofile.close()

numoptatom = np.zeros(len(optatoms))
numchcompatom = 0

for i in xrange(nsolatoms):
    for j in xrange(len(optatoms)):
        if atomtype[0][i] == optatoms[j]:
            numoptatom[j] += 1
            break
        elif atomtype[0][i] == chcomptype:
            numchcompatom += 1
            break

optlist = []
for i in xrange(nsolatoms):
    if atomtype[0][i] in optatoms:
        optlist.append(i)

e_solute, f_solute = readgauss.readsolute(fileprefix,n)
e_cluster, f_cluster = readgauss.readcluster(fileprefix,n)
e_water = readgauss.readwater(fileprefix,n)

e_qm = []
f_qm = []
for i in xrange(n):
    e_tmp = e_cluster[i] - e_water[i] - e_solute[i]
    e_tmp *= 627.509
    e_qm.append(e_tmp)
    f_tmp = []
    for j in xrange(nsolatoms):
        if atom_gauss[f_cluster[i][j][0]] in [atom_tinker[k] for k in optatoms]:
            fi_tmp = np.array([f_cluster[i][j][1], f_cluster[i][j][2], 
                f_cluster[i][j][3]]) 
            fi_tmp -= np.array([f_solute[i][j][1], f_solute[i][j][2], f_solute[i][j][3]])
            fi_tmp *= 1185.821
            f_tmp.append(fi_tmp)
    f_qm.append(f_tmp)

e_qm = np.array(e_qm)
f_qm = np.array(f_qm)
f_qm_com = np.sum(f_qm,axis=1)

nwatermol = []
for i in xrange(n):
    nwatermol.append((len(x[i]) - nsolatoms)/4)

boltz = np.zeros(n)
for i in xrange(n):
    boltz[i] = np.exp(-e_qm[i]/(nwatermol[i]*kb*T))

Q = np.sum(boltz)

e_avg = np.average(np.absolute(e_qm))
f_avg = np.average(np.absolute(f_qm),axis=0)

def chisq(prmlist):
    epsilon = prm_ep
    sigma = prm_sigma
    charge = prm_ch
    for i in xrange(len(optatoms)):
        epsilon[optatoms[i]] = prmlist[3*i]
        sigma[optatoms[i]] = prmlist[3*i+1]
        charge[optatoms[i]] = prmlist[3*i+2]
    if compcharge:
        charge[chcomptype] = syscharge
        for i in xrange(len(optatoms)):
            charge[chcomptype] -= charge[optatoms[i]] * numoptatom[i]
        charge[chcomptype] /= numchcompatom
    energy, force = efcalc.e_f(nsolatoms,optlist,epsilon,sigma,charge,n,x,y,z,
            atomtype,conv)
    denergy, dforce = efcalc.de_df(nsolatoms,optatoms,optlist,epsilon,sigma,charge,n,
            x,y,z,atomtype,conv)
    deltae = energy - e_qm
    force_com = np.sum(force,axis=1)
    deltaf_com = force_com - f_qm_com
    chi_e = 0.0
    chi_f = np.zeros((len(optlist),3))
    chi_e = np.sum(deltae**2 * boltz)/Q
    for i in xrange(n):
        chi_f += (deltaf_com[i])**2 * boltz[i]
    chi_f /= Q
    chi_sq = weight_e * chi_e + weight_f * np.average(chi_f)
    if compcharge:
        de_chcomp, df_chcomp = efcalc.de_df_chcomp(nsolatoms,chcomptype,epsilon,sigma,
                charge,n,x,y,z,atomtype,conv)
        for i in xrange(n):
            for j in xrange(len(optatoms)):
                denergy[i][3*j+2] -= de_chcomp[i] * numoptatom[j] / numchcompatom
                dforce[i][3*j+2] -= df_chcomp[i] * numoptatom[j] / numchcompatom
    dforce_com = np.sum(dforce,axis=2)
    dchi_e = np.zeros(3*len(optatoms))
    dchi_f = np.zeros((3*len(optatoms),3))
    for i in xrange(n):
        dchi_e += deltae[i] * denergy[i] * boltz[i]
        for j in xrange(len(dchi_f)):
            dchi_f[j] += (force_com[i] - f_qm_com[i]) * dforce_com[i][j] * boltz[i]
    dchi_e *= 2./Q
    dchi_f *= 2./Q
    dchi_sq = weight_e * dchi_e + weight_f * np.average(dchi_f,axis=1)
    return chi_sq, dchi_sq

optvalues = optimize.fmin_l_bfgs_b(chisq,initprms,bounds=prmbounds)

outfile = open(outname, 'w')
for i in xrange(len(optatoms)):
    outputline = atom_tinker[optatoms[i]] + ': {0:.6f} {1:.6f} {2:.6f}\n'.format(
            optvalues[0][3*i],optvalues[0][3*i+1],optvalues[0][3*i+2])
    outfile.write(outputline)

print optvalues
