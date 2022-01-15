#!/usr/bin/python
# coding: UTF-8
# This is created 2020/04/20 by Y. Shinohara
# This is lastly modified 2020/04/20 by Y. Shinohara
from modules.constants import *

def plot_potential(plt,cm, param, vx):
    plt.figure()
    plt.xlim(0.0, param.a)
    plt.ylabel('Energy [a.u.]')
    plt.plot(param.x,vx,label='The local potential')
    plt.grid()
    plt.legend()
    plt.show()
#
def plot_band(plt,cm, param, epsbk, Nbandmax = 4):
    plt.figure()
    plt.xlim(-0.5*param.b, 0.5*param.b)
    plt.xlabel('$k$ [a.u.]')
    plt.ylabel('$\epsilon_{bk}$ [eV]')
    for ib in range(Nbandmax):
        plt.plot(param.k,epsbk[ib,:]*Hartree)
    plt.grid()
    plt.show()
#
def plot_AE(plt,cm, param, t, A, E):
    plt.figure()
    plt.xlabel('Time [fs]')
    plt.ylabel('Vector potential [/nm]')
    plt.xlim(0.0,np.amax(t)*Atomtime)
    plt.plot(t*Atomtime,A/aB)
    plt.plot(t*Atomtime,param.b*np.ones(param.Nt)/2.0/aB)
    plt.plot(t*Atomtime,-param.b*np.ones(param.Nt)/2.0/aB)
    plt.grid()
    plt.show()
#
    plt.figure()
    plt.title('Electric field')
    plt.xlabel('Time [fs]')
    plt.ylabel('Field strength [V/nm]')
    plt.xlim(0.0,np.amax(t)*Atomtime)
    plt.plot(t*Atomtime,E*Atomfield)
    plt.grid()
    plt.show()
#
def plot_RT(plt,cm, t,J,Ene):
    plt.figure()
    plt.xlabel('Time [fs]')
    plt.ylabel('Current density [a.u.]')
    plt.xlim(0.0,np.amax(t)*Atomtime)
    plt.plot(t*Atomtime,J)
    plt.grid()
    plt.show()
#
    plt.figure()
    plt.title('Energy of the system')
    plt.xlabel('Time [fs]')
    plt.ylabel('Energy [eV]')
    plt.xlim(0.0,np.amax(t)*Atomtime)
    plt.plot(t*Atomtime,Ene*Hartree)
    plt.grid()
    plt.show()
