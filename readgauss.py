#! /usr/bin/env python

import re
import sys
import numpy as np

def readsolute(fileprefix,nfiles,istranstate,nsolatoms):
    forces = []
    energy = []
    imode = np.zeros((nfiles,nsolatoms,3))
    for i in xrange(nfiles):
        count = 0
        nforces = []
        forcematrix = 0
        freqblock = 0
        nmodes = 0
        solutefilename = 'gauss_' + fileprefix + '/' + fileprefix + '_' + str(i) + '.log'
        solutefile = open(solutefilename, 'r')
        for line in solutefile:
            line = line.strip()
            if not forcematrix:
                if re.search('Hartrees/Bohr', line):
                    forcematrix = 1
                elif not freqblock:
                    if re.match('SCF Done', line):
                        edata = re.split('\s+', line)
                    elif istranstate and re.match('Harmonic frequencies', line):
                        freqblock = 1
                elif not nmodes:
                    line = re.split('\s+',line)
                    if line[0] == 'Frequencies':
                        assert float(line[2]) < 0
                    elif line[0] == 'Coord':
                        nmodes += 1
                elif (0 < nmodes < 3*nsolatoms+1):
                    line = re.split('\s+',line)
                    imode[i][int(line[1])-1][int(line[0])-1] = float(line[3])
                    nmodes += 1
            else:
                if re.match('---.', line):
                    count += 1
                elif count == 1:
                    tmpforce = re.split('\s+', line)
                    nforces.append([int(tmpforce[1]), float(tmpforce[2]),
                        float(tmpforce[3]), float(tmpforce[4])])
                elif count == 2: break
        solutefile.close()
        energy.append(float(edata[4]))
        forces.append(nforces)
    return np.array(energy), forces, np.array(imode)

def readimplicit(fileprefix,nfiles,nsolatoms):
    dipole = np.zeros((nfiles,3))
    quadrupole = np.zeros((nfiles,6))
    traceless = False
    for i in xrange(nfiles):
        count = 0
        nforces = []
        forcematrix = 0
        freqblock = 0
        nmodes = 0
        multipole = False
        implicitfilename = ('gauss_' + fileprefix + '/' + fileprefix + '_implicit_' + 
                str(i) + '.log')
        implicitfile = open(implicitfilename, 'r')
        for line in implicitfile:
            line = line.strip()
            if re.match('X=', line):
                line = re.split('\s+',line)
                dipole[i] = np.array([float(line[1]),float(line[3]),float(line[5])])
            elif re.match('Traceless Quadrupole', line):
                multipole = True
            elif multipole:
                if re.match('XX=',line):
                    line = re.split('\s+',line)
                    quadrupole[i][0] = float(line[1])
                    quadrupole[i][1] = float(line[3])
                    quadrupole[i][2] = float(line[5])
                elif re.match('XY=',line):
                    line = re.split('\s+',line)
                    quadrupole[i][3] = float(line[1])
                    quadrupole[i][4] = float(line[3])
                    quadrupole[i][5] = float(line[5])
                    break
    return dipole, quadrupole

def readcluster(fileprefix,nfiles,nsolatoms):
    forces = []
    energy = []
    for i in xrange(nfiles):
        count = 0
        nforces = []
        forcematrix = 0
        firstcharge = False
        clusterfilename = ('gauss_' + fileprefix + '/' + fileprefix + '_cluster_' + 
                str(i) + '.log')
        clusterfile = open(clusterfilename, 'r')
        for line in clusterfile:
            line = line.strip()
            if not forcematrix and re.match('Counterpoise: corrected', line):
                edata = re.split('\s+', line)
            elif not forcematrix and re.search('Hartrees/Bohr', line):
                forcematrix = 1
            elif forcematrix:
                if re.match('---.', line):
                    count += 1
                elif count == 1:
                    tmpforce = re.split('\s+', line)
                    nforces.append([int(tmpforce[1]), float(tmpforce[2]), 
                        float(tmpforce[3]), float(tmpforce[4])])
                elif count == 2: break
        clusterfile.close()
        energy.append(float(edata[4]))
        forces.append(nforces)
    return np.array(energy), forces
    

def readwater(fileprefix,nfiles):
    energy = []
    for i in xrange(nfiles):
        waterfilename = ('gauss_' + fileprefix + '/' + fileprefix + '_water_' + str(i)
        + '.log')
        waterfile = open(waterfilename, 'r')
        for line in waterfile:
            line = line.strip()
            if re.match('SCF Done', line):
                edata = re.split('\s+', line)
        waterfile.close()
        energy.append(float(edata[4]))
    return np.array(energy)

