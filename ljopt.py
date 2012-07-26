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
weight_e = 15.19562
weight_f = 1.293817
weight_t = 1.0

fileprefix = sys.argv[1]

atom_tinker = {901:'C', 911:'H', 921:'Cl', 55:'O', 56:'H', 57:''}
atom_gauss = {1:'H', 6:'C', 8:'O', 17:'Cl'}
mass = {'H':1.008,'C':12.000,'O':15.999,'Cl':35.453}

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
istranstate = (re.split('\s+',ljofile.readline().strip())[0].lower() == 'true')
syscharge = float(chargedata[0])
chcomptype = int(chargedata[1])

solutemass = 0.0

for i in xrange(nsolatoms):
    solutemass += mass[atom_tinker[atomtype[0][i]]]

com = np.zeros((n,3))

for i in xrange(n):
    for j in xrange(nsolatoms):
        com[i] += np.array([x[i][j],y[i][j],z[i][j]])*mass[atom_tinker[atomtype[i][j]]]

com /= solutemass

r_com = np.zeros((n,nsolatoms,3))

for i in xrange(n):
    for j in xrange(nsolatoms):
        r_com[i][j] = np.array([x[i][j],y[i][j],z[i][j]]) - com[i]

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
        if atomtype[0][i] == chcomptype:
            numchcompatom += 1
            break
        elif atomtype[0][i] == optatoms[j]:
            numoptatom[j] += 1
            break

optlist = []
for i in xrange(nsolatoms):
    if atomtype[0][i] in optatoms:
        optlist.append(i)

e_solute, f_solute, imode = readgauss.readsolute(fileprefix,n,istranstate,nsolatoms)
e_cluster, f_cluster = readgauss.readcluster(fileprefix,n,nsolatoms)
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

t_qm = np.zeros((n,3))

for i in xrange(n):
    for j in xrange(nsolatoms):
        t_qm[i] += np.cross(r_com[i][j],f_qm[i][j])

nwatermol = []
for i in xrange(n):
    nwatermol.append((len(x[i]) - nsolatoms)/4)

boltz = np.zeros(n)
for i in xrange(n):
    boltz[i] = np.exp(-e_qm[i]/(nwatermol[i]*kb*T))

Q = np.sum(boltz)

e_avg = np.average(np.absolute(e_qm))
f_avg = np.average(np.absolute(f_qm_com))
t_avg = np.average(np.absolute(t_qm))

if t_avg == 0.0:
    t_avg = 1.0

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
    torque = np.zeros((n,3))
    for i in xrange(n):
        for j in xrange(nsolatoms):
            torque[i] += np.cross(r_com[i][j],force[i][j])
    deltat = torque - t_qm
    chi_e = np.sum(deltae**2 * boltz)/(Q*e_avg**2)
    chi_f = np.zeros(3)
    chi_t = np.zeros(3)
    for i in xrange(n):
        chi_f += deltaf_com[i]**2 * boltz[i]
        chi_t += deltat[i]**2 * boltz[i]
    chi_f /= Q * f_avg**2
    chi_t /= Q * t_avg**2
    print chi_e, chi_f, chi_t
    chi_sq = (weight_e * chi_e + weight_f * np.average(chi_f) + 
            weight_t * np.average(chi_t)) 
    if compcharge:
        de_chcomp, df_chcomp = efcalc.de_df_chcomp(nsolatoms,chcomptype,epsilon,sigma,
                charge,n,x,y,z,atomtype,conv)
        for i in xrange(n):
            for j in xrange(len(optatoms)):
                denergy[i][3*j+2] -= de_chcomp[i] * numoptatom[j] / numchcompatom
                dforce[i][3*j+2] -= df_chcomp[i] * numoptatom[j] / numchcompatom
    dforce_com = np.sum(dforce,axis=2)
    dtorque = np.zeros((n,3*len(optatoms),3))
    for i in xrange(n):
        for j in xrange(3*len(optatoms)):
            for k in xrange(len(optlist)):
                dtorque[i][j] += np.cross(r_com[i][k],dforce[i][j][k])
    dchi_e = np.zeros(3*len(optatoms))
    dchi_f = np.zeros((3*len(optatoms),3))
    dchi_t = np.zeros((3*len(optatoms),3))
    for i in xrange(n):
        dchi_e += deltae[i] * denergy[i] * boltz[i]
        for j in xrange(3*len(optatoms)):
            dchi_f[j] += deltaf_com[i] * dforce_com[i][j] * boltz[i]
            dchi_t[j] += deltat[i] * dtorque[i][j] * boltz[i]
    dchi_e *= 2./(Q*e_avg**2)
    dchi_f *= 2./(Q*f_avg**2)
    dchi_t *= 2./(Q*t_avg**2)
    dchi_sq = (weight_e * dchi_e + weight_f * np.average(dchi_f,axis=1) +
            weight_t * np.average(dchi_t,axis=1))
    return chi_sq, dchi_sq

optvalues = optimize.fmin_l_bfgs_b(chisq,initprms,bounds=prmbounds)

outfile = open(outname, 'w')
for i in xrange(len(optatoms)):
    outputline = atom_tinker[optatoms[i]] + ': {0:.6f} {1:.6f} {2:.6f}\n'.format(
            optvalues[0][3*i],optvalues[0][3*i+1],optvalues[0][3*i+2])
    outfile.write(outputline)

print optvalues
