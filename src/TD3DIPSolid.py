#!/usr/bin/python
# coding: UTF-8
# This is created 2022/01/16 by Y. Shinohara
# This is lastly modified xxxx/xx/xx by Y. Shinohara #This part is highly doubtable because of my lazyness
import time
ts = time.time()
ver = '0.0.0'
code_name = 'TD3DIPSolid-'+ver
from modules.print_funcs import print_header, print_footer, print_midtime, print_endtime
print_header(code_name)
import sys
import numpy as np
import math
import ctypes as ct
from modules.constants import *
from modules.parameters import parameter_class
from modules.functions import * #Caution!! This should be modified 
from modules.plot_funcs import plot_potential, plot_band, plot_AE, plot_RT
from modules.RT_propagation import RT_propagation_class
from modules.parameters import parameter_class
#
param = parameter_class()
param.read_parameters()    #Initialization of the parameters and the replacement from the standard input#
param.grid_constructions() #
param.get_Nocc()           #
#
RTc = RT_propagation_class()
uGbk_forward = RTc.uGbk_forward(param.propagator_option, param.Fortlib_option)
RTc.Prep4Fortlib(param)

if (param.plot_figure_option): #Matplotlib is activated for the cluster_mode == True
    import matplotlib.pyplot as plt
    from matplotlib import cm #To include color map

#############################Prep. for the system########################
uGbk = np.zeros([np.prod(param.NG), np.prod(param.NG), np.prod(param.Nk)],dtype='complex128') #Wave function in reciprocal space
epsbk = np.zeros([np.prod(param.NG), np.prod(param.Nk)],dtype='float64') #Eigenvalue of the Hamiltonian
occbk = np.zeros([np.prod(param.NG), np.prod(param.Nk)],dtype='float64') #Occupation number

vx, vG, vGG, vGGk = get_vxvGvGGvGGk(param)
if(param.plot_figure_option):
    plot_potential(plt,cm, param, vx)
tGGk = get_tGGk(param,0.0)
hGGk = tGGk + vGGk

#Band calculation 
for ik in range(np.prod(param.Nk)):
    epsbk[:,ik], uGbk[:,:,ik] = np.linalg.eigh(hGGk[:,:,ik])
uGbk = uGbk/np.sqrt(param.a)*float(param.NG) #Normalization
Eg = np.amin(epsbk[param.Nocc,:])-np.amax(epsbk[param.Nocc - 1,:])
print('# Eg = '+str(Eg)+' a.u. = '+str(Hartree*Eg)+' eV')
if (param.plot_figure_option):
    plot_band(plt,cm, param, epsbk)

print('# Band calculation is done properly.    ')
print('######################################')

sys.exit()
if (param.temperature < 0.0):
    occbk[0:param.Nocc,:] = 2.0/float(param.Nk)
    uGbk = 1.0*uGbk[:,0:param.Nocc,:] #This "1.0*" is mandatory for calling propelry Fortran subroutin via ctypes, MAKING NEW ARRAY specification
    occbk = 1.0*occbk[0:param.Nocc,:]
    print('# The system is assumed to be zero-temperature insulator. ')
else :
    print('# ERROR: Currenty, finite-temperature occupation distribution is not supported.')
    sys.exit()


dns = occ_u2dns(param,occbk,uGbk)
print('## Check for dns at initial, '+str(np.sum(dns)*param.H))
J = occ_u_A2J(param,occbk,uGbk,0.0) #Matter current, namely negative sign need for the charge current
print('## Check for current at initial, '+str(J))
Ene = occ_u_h2Ene(param,occbk,uGbk,hGGk)
print('## Check for Ene, '+str(Ene))
print('# System energy at initial:',Ene, '[a.u.] =',Ene*Hartree, ' [eV]')

#sys.exit()

#############################Prep. for RT################################
t, A, E = Make_Efield(param)
if (param.PC_option):
    Eave = ft2ftave(E)
    Aave = ft2ftave(A)
nv = np.zeros([param.Nt],dtype=np.float64)
nc = np.zeros([param.Nt],dtype=np.float64)
Ene = np.zeros([param.Nt],dtype=np.float64)
J = np.zeros([param.Nt],dtype=np.float64)

if (np.amax(t) < np.amax(param.Tpulse)):
    print('# WARNING: max(t) is shorter than Tpulse')
        
if (param.plot_figure_option):
    plot_AE(plt,cm, param,t,A,E) #Plot shape of the electric field

tt = time.time()
print_midtime(ts,tt)
#sys.exit()

#############################RT calculation##############################
#Time-propagation
for it in range(param.Nt):
    J[it] = occ_u_A2J(param,occbk,uGbk,A[it])
    Ene[it] = occ_u_h2Ene(param,occbk,uGbk,hGGk)
    if (param.PC_option):
        tGGk = get_tGGk(param,Aave[it])
    else:
        tGGk = get_tGGk(param,A[it])
    hGGk = tGGk + vGGk
    uGbk = uGbk_forward(param, uGbk, hGGk, tGGk, vx)
    if (it%1000 == 0):
        dns = occ_u2dns(param,occbk,uGbk)
        print(it,np.sum(dns)*param.H, J[it], Ene[it])

print('# System energy at end:',Ene[param.Nt-1], '[a.u.] =',Ene[param.Nt-1]*Hartree, ' [eV]')
print('# Absorbed energy:',Ene[param.Nt-1]-Ene[0], '[a.u.] =',(Ene[param.Nt-1]-Ene[0])*Hartree, ' [eV]')

if (param.plot_figure_option):
    plot_RT(plt,cm, t,J,Ene) #Plot data obtained in real-time evolution, J, Ene

te = time.time()
print_endtime(ts,tt,te,param.Nt)

if (param.minimal_output):
    np.savez(param.sys_name+'_tJEne.npz', t=t, J=J, Ene=Ene)
else:
    np.savez(param.sys_name+'_RTall.npz', t=t, E=E, A=A, J=J, Ene=Ene)

print('# Whole calculation is done properly.    ')
print_footer() 
sys.exit()

