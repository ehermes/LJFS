#! /usr/bin/env python

import numpy as np

def elj12(ep1,ep2,sig1,sig2,rvec):
    epsilon = np.sqrt(ep1*ep2)
    sigma = np.sqrt(sig1*sig2)
    rsq = np.dot(rvec,rvec)
    e_lj12 = 4*epsilon*(sigma**12)*(rsq**-6)
    return e_lj12

def elj6(ep1,ep2,sig1,sig2,rvec):
    epsilon = np.sqrt(ep1*ep2)
    sigma = np.sqrt(sig1*sig2)
    rsq = np.dot(rvec,rvec)
    e_lj6 = -4*epsilon*(sigma**6)*(rsq**-3)
    return e_lj6

def ecoul(conv,q1,q2,rvec):
    rsq = np.dot(rvec,rvec)
    rmag = np.sqrt(rsq)
    e_coul = conv * q1 * q2 / rmag
    return e_coul

def flj12(ep1,ep2,sig1,sig2,rvec):
    epsilon = np.sqrt(ep1*ep2)
    sigma = np.sqrt(sig1*sig2)
    f_lj12 = np.array([0.0,0.0,0.0])
    rsq = np.dot(rvec,rvec)
    f_lj12 = (48.0*epsilon*rvec)*(sigma**12)*(rsq**-7)
    return f_lj12

def flj6(ep1,ep2,sig1,sig2,rvec):
    epsilon = np.sqrt(ep1*ep2)
    sigma = np.sqrt(sig1*sig2)
    f_lj6 = np.array([0.0,0.0,0.0])
    rsq = np.dot(rvec,rvec)
    f_lj6 = (-24*epsilon*rvec)*(sigma**6)*(rsq**-4)
    return f_lj6

def fcoul(conv,q1,q2,rvec):
    f_coul = np.array([0.0,0.0,0.0])
    rsq = np.dot(rvec,rvec)
    rmag = np.sqrt(rsq)
    f_coul= conv * q1 * q2 * (rvec/rmag) / rsq
    return f_coul

def ddepelj12(ep1,ep2,sig1,sig2,rvec):
    depsilon = 0.5*np.sqrt(ep2/ep1)
    sigma = np.sqrt(sig1*sig2)
    rsq = np.dot(rvec,rvec)
    dedep_lj12 = 4*depsilon*(sigma**12)*(rsq**-6)
    return dedep_lj12

def ddsigelj12(ep1,ep2,sig1,sig2,rvec):
    epsilon = np.sqrt(ep1*ep2)
    sigma = np.sqrt(sig1*sig2)
    dsigma = 0.5*np.sqrt(sig2/sig1)
    rsq = np.dot(rvec,rvec)
    dedsig_lj12 = 48*epsilon*(sigma**11)*(rsq**-6)*dsigma
    return dedsig_lj12

def ddepelj6(ep1,ep2,sig1,sig2,rvec):
    depsilon = 0.5*np.sqrt(ep2/ep1)
    sigma = np.sqrt(sig1*sig2)
    rsq = np.dot(rvec,rvec)
    dedep_lj6 = -4*depsilon*(sigma**6)*(rsq**-3)
    return dedep_lj6

def ddsigelj6(ep1,ep2,sig1,sig2,rvec):
    epsilon = np.sqrt(ep1*ep2)
    sigma = np.sqrt(sig1*sig2)
    dsigma = 0.5*np.sqrt(sig2/sig1)
    rsq = np.dot(rvec,rvec)
    dedsig_lj6 = -24*epsilon*(sigma**5)*(rsq**-3)*dsigma
    return dedsig_lj6

def ddqecoul(conv,q1,q2,rvec):
    rsq = np.dot(rvec,rvec)
    rmag = np.sqrt(rsq)
    dedq_coul = conv * q2 / rmag
    return dedq_coul

def ddepflj12(ep1,ep2,sig1,sig2,rvec):
    depsilon = 0.5*np.sqrt(ep2/ep1)
    sigma = np.sqrt(sig1*sig2)
    dfdep_lj12 = np.array([0.0,0.0,0.0])
    rsq = np.dot(rvec,rvec)
    dfdep_lj12 = (48.0*depsilon*rvec)*(sigma**12)*(rsq**-7)
    return dfdep_lj12

def ddsigflj12(ep1,ep2,sig1,sig2,rvec):
    epsilon = np.sqrt(ep1*ep2)
    sigma = np.sqrt(sig1*sig2)
    dsigma = 0.5*np.sqrt(sig2/sig1)
    dfdsig_lj12 = np.array([0.0,0.0,0.0])
    rsq = np.dot(rvec,rvec)
    dfdsig_lj12 = (576.0*epsilon*rvec)*(sigma**11)*(rsq**-7)*dsigma
    return dfdsig_lj12

def ddepflj6(ep1,ep2,sig1,sig2,rvec):
    depsilon = 0.5*np.sqrt(ep2/ep1)
    sigma = np.sqrt(sig1*sig2)
    dfdep_lj6 = np.array([0.0,0.0,0.0])
    rsq = np.dot(rvec,rvec)
    dfdep_lj6 = (-24*depsilon*rvec)*(sigma**6)*(rsq**-4)
    return dfdep_lj6

def ddsigflj6(ep1,ep2,sig1,sig2,rvec):
    epsilon = np.sqrt(ep1*ep2)
    sigma = np.sqrt(sig1*sig2)
    dsigma = 0.5*np.sqrt(sig2/sig1)
    dfdsig_lj6 = np.array([0.0,0.0,0.0])
    rsq = np.dot(rvec,rvec)
    dfdsig_lj6 = (-144*epsilon*rvec)*(sigma**5)*(rsq**-4)*dsigma
    return dfdsig_lj6

def ddqfcoul(conv,q1,q2,rvec):
    dfdq_coul = np.array([0.0,0.0,0.0])
    rsq = np.dot(rvec,rvec)
    rmag = np.sqrt(rsq)
    dfdq_coul= conv * q2 * (rvec/rmag) / rsq
    return dfdq_coul
