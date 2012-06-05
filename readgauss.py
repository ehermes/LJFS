#! /usr/bin/env python

import re
import sys
import numpy as np

def readsolute(fileprefix):
    forces = []
    forcematrix = 0
    solutefilename = 'gauss_' + fileprefix + '/' + fileprefix + '.log'
    solutefile = open(solutefilename, 'r')
    for line in solutefile:
        line = line.strip()
        if not forcematrix and re.match('SCF Done', line):
            edata = re.split('\s+', line)
        elif not forcematrix and re.search('Hartrees/Bohr', line):
            forcematrix = 1
        elif forcematrix:
            if re.match('---.', line):
                count += 1
            elif count == 1:
                tmpforce = re.split('\s+', line)
                forces.append([int(tmpforce[1]), float(tmpforce[2]), float(tmpforce[3]),
                    float(tmpforce[4])])
            elif count == 2: break
    energy = float(edata[4])
    return [energy, forces]

def readcluster(fileprefix,nfiles):
    forces = []
    energy = []
    for i in xrange(nfiles):
        forcematrix = 0
        clusterfilename = ('gauss_' + fileprefix + '/' + fileprefix + '_cluster_' + i +
                '.log')
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
        energy.append(float(edata[4]))
        forces.append(nforces)
    return [energy, forces]
    

def readwater(fileprefix,nfiles):
    energy = []
    for i in xrange{nfiles):
        waterfilename = 'gauss_' + fileprefix + '/' + fileprefix + '_water_' + i + '.log'
        waterfile = open(waterfilename, 'r')
        for line in waterfile:
            line = line.strip()
            if re.match('SCF Done', line):
                edata = re.split('\s+', line)
        energy.append(float(edata[4]))
    return energy

