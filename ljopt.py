#! /usr/bin/env python

import sys
import copy
import numpy as np
import readarc
import readkey
import readgauss
import analyze
import efcalc

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
    if atomtype[0][i] == optatomname:
        optatoms.append(i)

if optatoms == []:
    print "No atoms to optimize parameters for!"
    sys.exit()

[e_solute, f_solute] = readgauss.readsolute(fileprefix)
[e_cluster, f_cluster] = readgauss.readcluster(fileprefix,len(n))
e_water = readgauss.readwater(fileprefix,len(n))

e_qm = []
f_qm = []
for i in xrange(len(n)):
    e_tmp = e_cluster[i] - e_water[i] - e_solute
    e_tmp *= 627.509
    e_qm.append(e_tmp)
    f_tmp = []
    for j in xrange(len(f_cluster[i])):
        if atom_gauss[f_cluster[i][j][0]] == optatomname:
            fi_tmp = np.array([f_cluster[i][j][1], f_cluster[i][j][2], 
                f_cluster[i][j][3]]) - np.array([f_solute[i][j][1], 
                    f_solute[i][j][2], f_solute[i][j][3]])
            fi_tmp *= 1185.821
            f_tmp.append(fi_tmp)
    f_qm.append(f_tmp)

optnum = 921

def chisq(e=prm_ep[optnum],s=prm_sigma[optnum]):
    epsilon = prm_ep
    sigma = prm_sigma
    charge = prm_ch
    epsilon[optnum] = e
    sigma[optnum] = s
    [energy, force] = e_f(epsilon,sigma,charge)

def gradchisq(e=prm_ep[optnum],s=prm_sigma[optnum]):
    epsilon = prm_ep
    sigma = prm_sigma
    charge = prm_ch
    epsilon[optnum] = e
    sigma[optnum] = s
    [denergy, dforce] = de_df(epsilon,sigma,charge)
