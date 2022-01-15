# coding: UTF-8
# This is created 2020/04/17 by Y. Shinohara
# This is lastly modified 2020/05/20 by Y. Shinohara
import sys
from modules.constants import pi,tpi, Atomtime, Hartree, Atomfield, aB, halfepsc
import numpy as np

class parameter_class:
    def __init__(self):
        # Default values for the parameters
        ## system
        self.sys_name = ''            #Prefix of the output files
        self.plot_figure_option = False
        self.PC_option = True         #Predictor-corrector option
        self.minimal_output = True    #A flag to write-out minimal data or not
        self.Fortlib_option = False   #A flag to use Fortran acceleration by using "Fortlib"
        ## System info.
        self.a = np.array([[8.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 8.0]])     #The lattice constant
        self.alen = None
        self.b = None                 #The reciprocal lattice vector, 
        self.blen = None   
        self.vol = None               #Volume of the cell
        self.v0 = 0.37                #The depth of the local potential
        self.temperature = -1.0       #Electron temperature, negative value leading to the zero-temperature
        self.Nave = 4.0               #Number of electron in the cell, doubly degenerated due to the spin
        self.Nocc = None              #Number of occupied levels, int(Nave/2)
        ## Numerical discretization
        self.NG = np.array([12, 12, 12])    #Number of spatial grid
        self.dr = None                 #Spatial grid sizes
        self.r = None                 #Spatial grid
        self.Nk = np.array([4, 4, 4]) #Number of Brillouin zone sampling
        self.k = None                 #k-grid
        ## Time propagation
        self.propagator_option = 'exp'        #The size of the time-step
        self.dt = 5.0e-1              #The size of the time-step
        self.Nt = 4000                #The number of time steps
        self.NKS = 8                  #Number of the subspae in KS construction with "KS" propagator
        ## Field parameters
        self.Ncolor = 1               #Number of color for the field
        self.omegac = 0.3875/Hartree  #Photon energy  [a.u.]
        self.phi_CEP = 0.25           #Carrier envelope phase [2 pi]
        self.Tpulse = 40.0/Atomtime   #A parameter for pulse duration [a.u.]
        self.nenvelope = 4            #Power for the sing envelope function 
        self.E0 = 0.9*1.843/Atomfield #Field strength [a.u.]

    def read_parameters(self):
        argv = sys.argv
        argc = len(argv)
        #Reading data from the standard input
        if (argc == 1):
            print('# The default parameters are chosen.')
        elif (argc == 2):
            print('# Name of input file is "'+argv[1]+'".')
            f = open(argv[1],'r')
            lines = f.readlines()
            f.close
            Nlen = len(lines)
            text = [0]*Nlen
            for i in range(Nlen):
                text[i] = lines[i].strip()
            for i in range(Nlen):
                if (str(text[i]) == 'sys_name'):
                    self.sys_name = text[i+1].split()[0]
                if (str(text[i]) == 'plot_figure_option'):
                    if (str(text[i+1]) == 'True'):
                        self.plot_figure_option = True
                    else:
                        self.plot_figure_option = False
                if (str(text[i]) == 'PC_option'):
                    if (str(text[i+1]) == 'True'):
                        self.PC_option = True
                    else:
                        self.PC_option = False
                if (str(text[i]) == 'minimal_output'):
                    if (str(text[i+1]) == 'True'):
                        self.minimal_output = True
                    else:
                        self.minimal_output = False
                if (str(text[i]) == 'Fortlib_option'):
                    if (str(text[i+1]) == 'True'):
                        self.Fortlib_option = True
                    else:
                        self.Fortlib_option = False
                if (str(text[i]) == 'a'):
                    temp = text[i+1].split()
                    self.a[0] = float(str(temp[0]))
                    self.a[1] = float(str(temp[1]))
                    self.a[2] = float(str(temp[2]))
                if (str(text[i]) == 'v0'):
                    self.v0 = float(str(text[i+1]))
                if (str(text[i]) == 'temperature'):
                    self.temperature = float(str(text[i+1]))
                if (str(text[i]) == 'Nave'):
                    self.Nave = float(str(text[i+1]))
                    
                if (str(text[i]) == 'NG'):
                    temp = text[i+1].split()
                    self.NG[0] = int(str(temp[0]))
                    self.NG[1] = int(str(temp[1]))
                    self.NG[2] = int(str(temp[2]))
                if (str(text[i]) == 'Nk'):
                    temp = text[i+1].split()
                    self.Nk[0] = int(str(temp[0]))
                    self.Nk[1] = int(str(temp[1]))
                    self.Nk[2] = int(str(temp[2]))
                if (str(text[i]) == 'propagator_option'):
                    self.propagator_option = text[i+1].split()[0]
                if (str(text[i]) == 'dt'):
                    self.dt = float(str(text[i+1]))
                if (str(text[i]) == 'Nt'):
                    self.Nt = int(str(text[i+1]))
                if (str(text[i]) == 'NKS'):
                    self.NKS = int(str(text[i+1]))

                if (str(text[i]) == 'Ncolor'):
                    self.Ncolor = int(str(text[i+1]))
                if (str(text[i]) == 'Tpulse'):
                    self.Tpulse = float(str(text[i+1].split()[0]))
                if (str(text[i]) == 'omegac'):
                    self.omegac = float(str(text[i+1].split()[0]))
                if (str(text[i]) == 'phi_CEP'):
                    self.phi_CEP = float(str(text[i+1].split()[0]))
                if (str(text[i]) == 'nenvelope'):
                    self.nenvelope = int(str(text[i+1].split()[0]))
                if (str(text[i]) == 'E0'):
                    self.E0 = float(str(text[i+1].split()[0]))
            if (self.Ncolor >= 2):
                self.omegac = self.omegac*np.ones([self.Ncolor],dtype='float64')
                self.phi_CEP = self.phi_CEP*np.ones([self.Ncolor],dtype='float64')
                self.Tpulse = self.Tpulse*np.ones([self.Ncolor],dtype='float64')
                self.nenvelope = self.nenvelope*np.ones([self.Ncolor],dtype='int32')
                self.E0 = self.E0*np.ones([self.Ncolor],dtype='float64')
                for i in range(Nlen):
                    if (str(text[i]) == 'Tpulse'):
                        temp = text[i+1].split()
                        self.Tpulse = np.array(temp,dtype='float64')
                        if (len(self.Tpulse) != self.Ncolor):
                            print('Error: Number of argmeunt in Tpulse is wrong.')
                            sys.exit()
                    if (str(text[i]) == 'omegac'):
                        temp = text[i+1].split()
                        self.omegac = np.array(temp,dtype='float64')
                        if (len(self.omegac) != self.Ncolor):
                            print('Error: Number of argmeunt in omegac is wrong.')
                            sys.exit()
                    if (str(text[i]) == 'phi_CEP'):
                        temp = text[i+1].split()
                        self.phi_CEP = np.array(temp,dtype='float64')
                        if (len(self.phi_CEP) != self.Ncolor):
                            print('Error: Number of argmeunt in phi_CEP is wrong.')
                            sys.exit()
                    if (str(text[i]) == 'nenvelope'):
                        temp = text[i+1].split()
                        self.nenvelope = np.array(temp,dtype='int32')
                        if (len(self.nenvelope) != self.Ncolor):
                            print('Error: Number of argmeunt in nenvelope is wrong.')
                            sys.exit()
                    if (str(text[i]) == 'E0'):
                        temp = text[i+1].split()
                        self.E0 = np.array(temp,dtype='float64')
                        if (len(self.E0) != self.Ncolor):
                            print('Error: Number of argmeunt in E0 is wrong.')
                            sys.exit()
        else:
            print('Error: Number of argmeunt is wrong.')
            sys.exit()
        
        if (self.NKS > np.prod(self.NG)):
            print('# NOTICE: NKS is larger than the size of the Hilbert space, NG. ')
            print('# NKS is replaced by NG ')
            self.NKS = 1*np.prod(self.NG)

        self.alen = np.zeros(3, dtype='float64')
        self.alen[0] = np.linalg.norm(self.a[:,0])
        self.alen[1] = np.linalg.norm(self.a[:,1])
        self.alen[2] = np.linalg.norm(self.a[:,2])
        print('#=====Print the parmeters')
        print('# sys_name =', self.sys_name)
        print('# plot_figure_option =', self.plot_figure_option)
        print('# PC_option =', self.PC_option)
        print('# minimal_output =', self.minimal_output)
        print('# Fortlib_option =', self.Fortlib_option)
        print('# ')
        print('# a =', self.a, ' [a.u.] =',self.a*aB, ' [nm]')
        print('# alen =', self.alen, ' [a.u.] =',self.alen*aB, ' [nm]')
        print('# v0 =', self.v0, ' [a.u.] =', self.v0*Hartree, ' [eV]') 
        print('# Nave =', self.Nave) 
        print('# ')
        print('# NG =', self.NG) 
        print('# Nk =', self.Nk) 
        print('# ')
        print('# propagator_option =', self.propagator_option) 
        print('# dt =', self.dt, ' [a.u.] =', self.dt*Atomtime, ' [fs]') 
        print('# 2pi/dt =', tpi/self.dt, ' [a.u.] =', tpi/self.dt*Hartree, ' [eV]') 
        print('# Nt =', self.Nt) 
        print('# Nt*dt =', self.Nt*self.dt, '[a.u.] =', self.Nt*self.dt*Atomtime, '[fs]')
        print('# 2pi/(Nt*dt) =', tpi/(self.Nt*self.dt), '[a.u.] =', tpi/(self.Nt*self.dt)*Hartree, '[eV]')
        print('# NKS =', self.NKS)
        print('# ')
        print('# Number of color: Ncolor = ', self.Ncolor)
        if (self.Ncolor == 1):
            print('# Tpulse =', self.Tpulse, ' [a.u.] =', self.Tpulse*Atomtime, ' [fs]') 
            print('# 2pi/Tpulse =', tpi/self.Tpulse, ' [a.u.] =', tpi/self.Tpulse*Hartree, ' [eV]') 
            print('# nenvelope =', self.nenvelope) 
            print('# omegac =', self.omegac, ' [a.u.] =', self.omegac*Hartree, ' [eV]') 
            print('# tpi/omegac =', tpi/self.omegac, ' [a.u.] =', tpi/self.omegac*Atomtime, ' [fs]') 
            print('# phi_CEP =', self.phi_CEP, ' [2 pi]')
            self.phi_CEP = tpi*self.phi_CEP
            print('# E0 =', self.E0, ' [a.u.] =', self.E0*Atomfield, ' [V/nm]') 
            print('# e*alen*E0 =', self.alen*self.E0, ' [a.u.] =', self.alen*self.E0*Hartree, ' [eV]')
            print('# E0/omegac =', self.E0/self.omegac, ' [a.u.] =', self.E0/self.omegac/aB, ' [/nm]') 
            print('# Peak intensity of the envelope: Imax =', (self.E0)**2, ' [a.u.] =', (self.E0)**2*halfepsc/1.0e9, ' [GW/cm^2]') 
        else:
            for icolor in range(self.Ncolor):
                print('# =====', icolor,'th color ========')
                print('# Tpulse =', self.Tpulse[icolor], ' [a.u.] =', self.Tpulse[icolor]*Atomtime, ' [fs]') 
                print('# 2pi/Tpulse =', tpi/self.Tpulse[icolor], ' [a.u.] =', tpi/self.Tpulse[icolor]*Hartree, ' [eV]') 
                print('# nenvelope =', self.nenvelope[icolor]) 
                print('# omegac =', self.omegac[icolor], ' [a.u.] =', self.omegac[icolor]*Hartree, ' [eV]') 
                print('# tpi/omegac =', tpi/self.omegac[icolor], ' [a.u.] =', tpi/self.omegac[icolor]*Atomtime, ' [fs]') 
                print('# phi_CEP =', self.phi_CEP[icolor], ' [2 pi]')
                self.phi_CEP[icolor] = tpi*self.phi_CEP[icolor]
                print('# E0 =', self.E0[icolor], ' [a.u.] =', self.E0[icolor]*Atomfield, ' [V/nm]') 
                print('# e*alen*E0 =', self.alen*self.E0[icolor], ' [a.u.] =', self.alen*self.E0[icolor]*Hartree, ' [eV]')

    def grid_constructions(self):
        self.vol = np.dot(self.a[:,0], np.cross(self.a[:,1], self.a[:,2]))
        self.b = np.zeros([3,3], dtype='float64')
        self.b[0] = tpi/self.vol*np.cross(self.a[:,1], self.a[:,2])
        self.b[1] = tpi/self.vol*np.cross(self.a[:,2], self.a[:,0])
        self.b[2] = tpi/self.vol*np.cross(self.a[:,0], self.a[:,1])
        self.blen = np.zeros(3, dtype='float64')
        self.blen[0] = np.linalg.norm(self.b[:,0])
        self.blen[1] = np.linalg.norm(self.b[:,1])
        self.blen[2] = np.linalg.norm(self.b[:,2])
        self.dr = self.alen/self.NG
        self.r = np.zeros([np.prod(self.NG),3],dtype='float64')
        l=0
        for i in range(self.NG[0]):
            for j in range(self.NG[1]):
                for k in range(self.NG[2]):
                    self.r[l,0] = i*self.a[0,0]/self.NG[0] + j*self.a[0,1]/self.NG[1] + k*self.a[0,2]/self.NG[2]
                    self.r[l,1] = i*self.a[1,0]/self.NG[0] + j*self.a[1,1]/self.NG[1] + k*self.a[1,2]/self.NG[2]
                    self.r[l,2] = i*self.a[2,0]/self.NG[0] + j*self.a[2,1]/self.NG[1] + k*self.a[2,2]/self.NG[2]
                    l += 1
        self.G = np.zeros([np.prod(self.NG),3],dtype='float64')
        l=0
        fa1 = np.fft.fftfreq(self.NG[0])*(np.float(self.NG[0])) #FFT freq array
        fa2 = np.fft.fftfreq(self.NG[1])*(np.float(self.NG[1]))
        fa3 = np.fft.fftfreq(self.NG[2])*(np.float(self.NG[2]))
        for i in range(self.NG[0]):
            for j in range(self.NG[1]):
                for k in range(self.NG[2]):
                    self.G[l,0] = fa1[i]*self.b[0,0] + fa1[j]*self.b[0,1] + fa3[k]*self.b[0,2]
                    self.G[l,1] = fa1[i]*self.b[1,0] + fa1[j]*self.b[1,1] + fa3[k]*self.b[1,2]
                    self.G[l,2] = fa1[i]*self.b[2,0] + fa1[j]*self.b[2,1] + fa3[k]*self.b[2,2]
                    l += 1
        #Brillouin zone construction
        self.k = np.zeros([np.prod(self.Nk),3],dtype='float64')
        l=0
        for i in range(self.Nk[0]):
            for j in range(self.Nk[1]):
                for k in range(self.Nk[2]):
                    self.k[l,0] = i*self.b[0,0]/self.Nk[0] + j*self.b[0,1]/self.Nk[1] + k*self.b[0,2]/self.Nk[2]
                    self.k[l,1] = i*self.b[1,0]/self.Nk[0] + j*self.b[1,1]/self.Nk[1] + k*self.b[1,2]/self.Nk[2]
                    self.k[l,2] = i*self.b[2,0]/self.Nk[0] + j*self.b[2,1]/self.Nk[1] + k*self.b[2,2]/self.Nk[2]
                    l += 1
        for i in range(3):
            self.k[:,i] = self.k[:,i] - np.average(self.k[:,i])
        print('# ')
        print('# b = 2pi/a =', self.b, ' [a.u.] =', self.b/aB, '[/nm]')
        print('# blen = 2pi/a =', self.blen, ' [a.u.] =', self.blen/aB, '[/nm]')
        print('# dr =', self.dr, ' [a.u.] =',self.dr*aB, ' [nm]')
        print('# (pi/dr)^2 =', (pi/self.dr)**2, ' [a.u.] =', (pi/self.dr)**2*Hartree, ' [eV]')
        print('# Nk*alen =', self.Nk*self.alen, ' [a.u.] =', (self.Nk*self.alen)*aB, ' [nm]')
        print('# ')

    def get_Nocc(self):
        self.Nocc = int(self.Nave/2.0)
        if (abs(self.Nave - 2.0*self.Nocc) > 1.0e-8):
            print('# The average number of particle in a cell: ', self.Nave)
            print('# ERROR: Currenty, metallic occupation is not supported.')
            print('# The average number should be even.')
            sys.exit()
            
