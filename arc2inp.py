#/usr/bin/env python

import os
import sys
import math
import re
import readarc
import numpy as np

fileprefix = sys.argv[1]

keywords = '$rem\nJOBTYPE SP\nEXCHANGE omegaB97X-D\nBASIS GENERAL\nPURECART 2222\nEFP TRUE\nCHELPG TRUE\n$end\n\n'

arcname = fileprefix + '.arc'
ljoname = fileprefix + '.ljo'

arc = readarc.readarc(arcname)

ljofile = open(ljoname, 'r')
nsolatoms = int(re.split('\s+',ljofile.readline().strip())[0])
charge = re.split('\s+',ljofile.readline().strip())[0]
ljofile.close()

n = arc[0]
x = arc[2]
y = arc[3]
z = arc[4]
atomtype = arc[5]

symbols = {900:'C', 901:'C', 902:'C', 910:'H', 911:'H', 912:'H', 920:'Cl',
        921:'Cl', 922:'Cl', 930:'O', 931:'O', 932:'O', 940:'H', 941:'H', 942:'H',
        55:'O', 56:'H', 57:''}
chgmult = str(int(float(charge))) + ' 1'

dirname = 'qchem_' + fileprefix

if not os.path.exists(dirname):
    os.makedirs(dirname)


water_ref = np.array([[0.,0.,0.0664326840],[0.,0.7532,-0.5271673160],
    [0.,-0.7532,-0.5271673160]])
water_mass = np.array([[15.999,15.999,15.999],[1.008,1.008,1.008],[1.008,1.008,1.008]])
water_totalmass = 18.015
tmp_x = np.cross(water_ref[2],water_ref[1])
x_ref = tmp_x / np.sqrt(np.dot(tmp_x,tmp_x))
z_ref = water_ref[0] / np.sqrt(np.dot(water_ref[0],water_ref[0]))
for i in xrange(n):
    nwater = (len(x[i]) - nsolatoms)/4
    efpname = fileprefix + '_' + str(i)
    efpfile = open(dirname + '/' + efpname + '.inp','w')
    efpfile.write('$molecule\n' + chgmult + '\n')
    for j in xrange(nsolatoms):
        efpfile.write(' {0:s}{1:d} {2:.10f} {3:.10f} {4:.10f}\n'.format(
            symbols[atomtype[i][j]], j, x[i][j], y[i][j], z[i][j]))
    efpfile.write('$end\n\n' + keywords + '$efp_fragments\n')
    for j in xrange(nwater):
        iwater = np.zeros((3,3))
        iwater[0] = np.array([x[i][4*j+nsolatoms],y[i][4*j+nsolatoms],
            z[i][4*j+nsolatoms]])
        iwater[1] = np.array([x[i][4*j+nsolatoms+1],y[i][4*j+nsolatoms+1],
            z[i][4*j+nsolatoms+1]])
        iwater[2] = np.array([x[i][4*j+nsolatoms+2],y[i][4*j+nsolatoms+2],
            z[i][4*j+nsolatoms+2]])
        iwater_com = np.sum((iwater * water_mass)/water_totalmass,axis=0)
        iwater_center = iwater - iwater_com
        tmp_iwx = np.cross(iwater_center[2],iwater_center[1])
        iwater_x = tmp_iwx / np.sqrt(np.dot(tmp_iwx,tmp_iwx))
        iwater_z = iwater_center[0]/np.sqrt(np.dot(iwater_center[0],iwater_center[0]))
        tmp_n = np.cross(z_ref,iwater_z)
        norm = tmp_n / np.sqrt(np.dot(tmp_n,tmp_n))
        tmp_anum = np.cross(x_ref,norm)
        alpha_num = np.sqrt(np.dot(tmp_anum,tmp_anum)) * np.sign(np.vdot(tmp_anum,z_ref))
        alpha_den = np.dot(x_ref,norm)
        alpha = np.arctan2(alpha_num,alpha_den)
        tmp_bnum = np.cross(z_ref,iwater_z)
        beta_num = np.sqrt(np.dot(tmp_bnum,tmp_bnum))
        beta_den = np.dot(z_ref,iwater_z)
        beta = np.arctan2(beta_num,beta_den)
        tmp_gnum = np.cross(norm,iwater_x)
        gamma_num = np.sqrt(np.dot(tmp_gnum,tmp_gnum)) * np.sign(np.vdot(tmp_gnum,iwater_z))
        gamma_den = np.dot(norm,iwater_x)
        gamma = np.arctan2(gamma_num,gamma_den)
        iwaterline = 'WATER_L {0:.8f} {1:.8f} {2:.8f} {3:.8f} {4:.8f} {5:.8f}\n'.format(
                iwater_com[0], iwater_com[1], iwater_com[2], alpha, beta, gamma)
        efpfile.write(iwaterline)
    efpfile.write('$end\n\n')
    efpfile.close()
