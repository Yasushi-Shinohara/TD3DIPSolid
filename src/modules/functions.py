# coding: UTF-8
#Relevant functions are written
# This is created 2020/04/17 by Y. Shinohara
# This is lastly modified 2020/05/20 by Y. Shinohara
import os
import sys
import math
import numpy as np
from modules.constants import *
#
def get_vxvGvGGvGGk(param):
    alpha = 5.0e-2
    beta = 5.0e-2
    gamma = 1.0e-1
    vx = -param.v0*(1.0 - np.cos(tpi*param.x/param.a))
    if (param.flat_length > 0.0):
        if (param.flat_length > param.a):
            print('The parameter, flat_length = '+str(param.flat_length)+' [a.u.].')
            print('ERROR: flat_length should be longer than the lattice constant, a.')
            sys.exit()
        for ig in range(param.NG):
            if (param.x[ig] < (param.a - param.flat_length)):
                vx[ig] = -param.v0*(1.0 - np.cos(tpi*param.x[ig]/(param.a - param.flat_length)))
            else:
                vx[ig] = 0.0
    vG = np.fft.fft(vx)/np.float(param.NG)
    vGG = np.zeros([param.NG, param.NG], dtype='complex128') 
    for ig1 in range(param.NG):
        gind = np.remainder(ig1 - np.arange(param.NG), param.NG)
        for ig2 in range(param.NG):
            igloc = gind[ig2]
            vGG[ig1,ig2] = vG[igloc] #This definition should be carefully checked 
    vGGk = np.zeros([param.NG, param.NG, param.Nk], dtype='complex128')
    for ik in range(param.Nk):
        vGGk[:,:,ik] = 1.0*vGG[:,:]
    return vx, vG, vGG, vGGk
#
def get_tGGk(param, A):
    tGGk = np.zeros([param.NG, param.NG, param.Nk],dtype='complex128')
    for ik in range(param.Nk):
        kpA = param.k[ik] + A
        for ig in range(param.NG):
            tGGk[ig, ig ,ik] = tGGk[ig, ig, ik] + 0.5*(param.G[ig] + kpA)**2
    return tGGk

#Relevant functions
def occ_u2dns(param,occbk,uGbk):
    dns = np.zeros(param.NG,dtype='float64')
    work = np.empty_like(uGbk[:,0,0])
    NBact = np.shape(uGbk)[1]
    for ik in range(param.Nk):
        for ib in range(NBact):
            work = np.fft.ifft(uGbk[:,ib,ik])
            dns = dns + occbk[ib,ik]*(np.abs(work))**2
    return dns

def occ_u_A2J(param,occbk,uGbk,A): #Exact formula should be checked=========
    NBact = np.shape(uGbk)[1]
    J = 0.0
    for ik in range(param.Nk):
        kpA = param.k[ik] + A
        for ib in range(NBact):
            J = J + occbk[ib,ik]*(np.sum(param.G[:]*(np.abs(uGbk[:,ib,ik]))**2)*param.a/float(param.NG**2) + kpA)
    return J/param.a

def occ_u_h2Ene(param,occbk,uGbk, hGGk): #Exact formula should be checked=========
    NBact = np.shape(uGbk)[1]
    Ene = 0.0
    hk = 1.0*hGGk
    for ik in range(param.Nk):
        hubG = np.dot(hk[:,:,ik], uGbk[:,:,ik])
        for ib in range(NBact):
            Ene = Ene + occbk[ib,ik]*np.real(np.vdot(uGbk[:,ib,ik],hubG[:,ib]))
    return Ene*param.a/float(param.NG**2) #orbital function is normalized to give correct number of particle in the cell.
#

def Make_Efield(param):
    t = np.zeros([param.Nt],dtype=np.float64)
    E = np.zeros([param.Nt],dtype=np.float64)
    A = np.zeros([param.Nt,param.Ncolor],dtype=np.float64)
    for it in range(param.Nt):
        t[it] = param.dt*it
    if (param.Ncolor == 1):
        icolor = 0
        for it in range(param.Nt):
            if (t[it] < param.Tpulse):
                A[it,icolor] = (param.E0/param.omegac)*(np.sin(pi*t[it]/param.Tpulse))**param.nenvelope*np.cos(param.omegac*(t[it] - 0.5*param.Tpulse) + param.phi_CEP)
    elif (param.Ncolor > 1):
        for icolor in range(param.Ncolor):
            for it in range(param.Nt):
                if (t[it] < param.Tpulse[icolor]):
                    A[it,icolor] = (param.E0[icolor]/param.omegac[icolor])*(np.sin(pi*t[it]/param.Tpulse[icolor]))**param.nenvelope[icolor]*np.cos(param.omegac[icolor]*(t[it] - 0.5*param.Tpulse[icolor]) + param.phi_CEP[icolor])
        A = np.sum(A,axis=1)
    else :
        print('ERROR: The parameter '+str(param.Ncolor)+' is improper.')
        sys.exit()

    for it in range(1,param.Nt-1):
        E[it] = -(A[it+1] - A[it-1])/2.0/param.dt
    E[0] = 2.0*E[1] - E[2]
    E[param.Nt-1] = 2.0*E[param.Nt-2] - E[param.Nt-3]
    print('# ') 
    print('# Actual field max =', np.amax(np.abs(E)), ' [a.u.] =', np.amax(np.abs(E))*Atomfield, ' [V/nm]') 
    print('# Corresponding intensity: ', np.amax(E**2), ' [a.u.] =', np.amax(E**2)*halfepsc/1.0e9, ' [GW/cm^2]') 
    print('# Fluence of the pulse: ', np.sum(E**2)*param.dt, ' [a.u.] =', np.sum(E**2)*param.dt*Atomfluence, ' [J/cm^2]') 
    print('# ') 
    return t, A, E

#
def ft2ftave(ft):
    ftave = 0.0*ft
    Nt = len(ft)
    for it in range(Nt-1):
        ftave[it] = 0.5*(ft[it] + ft[it+1])
    ftave[Nt - 1] = 2.0*ftave[Nt - 2] - ftave[Nt-3]
    return ftave
