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
#newbonds = copy.deepcopy(bonds)

#for i in xrange(n):
#    for j in xrange(len(x[i])):
#        for k in bonds[i][j]:
#            newbonds[i][k-1] += bonds[i][j]

keydata = readkey.readkey(keyname)
prm_sigma = keydata[0]
prm_ep = keydata[1]
prm_ch = keydata[2]
conv = keydata[3]

ljofile = open(ljoname, 'r')

nsolatoms = int(re.split('\s+',ljofile.readline().strip())[0])
chargedata = re.split('\s+',ljofile.readline().strip())
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
        prmbounds.append((prm_ep[int(data[0])]*0.1,prm_ep[int(data[0])]*10.0))
    else:
        prmbounds.append((prm_ep[int(data[0])],prm_ep[int(data[0])]))
    initprms.append(prm_sigma[int(data[0])])
    if int(data[2]):
        prmbounds.append((prm_sigma[int(data[0])]*0.1,prm_sigma[int(data[0])]*2.0))
    else:
        prmbounds.append((prm_sigma[int(data[0])],prm_sigma[int(data[0])]))
    initprms.append(prm_ch[int(data[0])])
    if int(data[3]):
        compcharge = True
        prmbounds.append((prm_ch[int(data[0])]-0.5,prm_ch[int(data[0])]+0.5))
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
    for j in xrange(nsolatoms):
        if atom_gauss[f_cluster[i][j][0]] in [atom_tinker[k] for k in optatoms]:
            fi_tmp = np.array([f_cluster[i][j][1], f_cluster[i][j][2], 
                f_cluster[i][j][3]]) 
            fi_tmp -= np.array([f_solute[j][1], f_solute[j][2], f_solute[j][3]])
            fi_tmp *= 1185.821
            f_tmp.append(fi_tmp)
    f_qm.append(f_tmp)

e_std = np.std(e_qm)
f_std = np.std(f_qm,axis=0)

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
    chi_e = 0.0
    chi_f = np.zeros((len(optlist),3))
    for i in xrange(n):
        chi_e += (energy[i] - e_qm[i])**2
        for j in xrange(len(optlist)):
            chi_f += (force[i][j] - f_qm[i][j])**2
    chi_e /= e_std**2
    chi_f /= f_std**2
    chi_sq = 1.00 * chi_e + 0.00 * (np.sum(chi_f)/3)
    if compcharge:
        dedq_chcomp = efcalc.dedq_chcomp(nsolatoms,chcomptype,epsilon,sigma,charge,n,
                x,y,z,atomtype,conv)
        for i in xrange(n):
            for j in xrange(len(optatoms)):
                denergy[i][3*j+2] -= dedq_chcomp[i] * numoptatom[j] / numchcompatom
    dchi_e = np.zeros(3*len(optatoms))
    dchi_f = np.zeros((len(optlist),3*len(optatoms),3))
    for i in xrange(n):
        dchi_e += (energy[i] - e_qm[i]) * denergy[i]
        for j in xrange(len(dchi_f)):
            for k in xrange(len(dchi_f[j])):
                dchi_f[j][k] += (force[i][j] - f_qm[i][j]) * dforce[i][j][k]
    dchi_e *= 2/(e_std**2)
    for i in xrange(len(dchi_f)):
        for j in xrange(len(dchi_f[i])):
            dchi_f[i][j] *= 2/(f_std[i]**2)
    dchi_sq = 1.00 * dchi_e + 0.00 * (np.sum(np.sum(dchi_f,axis=0),axis=1)/
            (3*len(optlist)))
    return chi_sq, dchi_sq

optvalues = optimize.fmin_l_bfgs_b(chisq,initprms,bounds=prmbounds)

outfile = open(outname, 'w')
for i in xrange(len(optatoms)):
    outputline = atom_tinker[optatoms[i]] + ': {0:.6f} {1:.6f} {2:.6f}\n'.format(
            optvalues[0][3*i],optvalues[0][3*i+1],optvalues[0][3*i+2])
    outfile.write(outputline)

print optvalues
