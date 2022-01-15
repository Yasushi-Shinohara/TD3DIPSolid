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
        self.a = 8.0                  #The lattice constant
        self.b = None                 #The reciprocal lattice vector, tpi/a
        self.v0 = 0.37                #The depth of the local potential
        self.flat_length = -1.0       #The length that potential is flat with zero-value
        self.temperature = -1.0       #Electron temperature, negative value leading to the zero-temperature
        self.Nave = 4.0               #Number of electron in the cell, doubly degenerated due to the spin
        self.Nocc = None              #Number of occupied levels, int(Nave/2)
        ## Numerical discretization
        self.NG = 12                  #Number of spatial grid
        self.H = None                 #Spatial grid size, a/NG
        self.x = None                 #Spatial grid
        self.Nk = 20                  #Number of Brillouin zone sampling
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
                    self.a = float(str(text[i+1]))
                if (str(text[i]) == 'flat_length'):
                    self.flat_length = float(str(text[i+1]))
                if (str(text[i]) == 'v0'):
                    self.v0 = float(str(text[i+1]))
                if (str(text[i]) == 'temperature'):
                    self.temperature = float(str(text[i+1]))
                if (str(text[i]) == 'Nave'):
                    self.Nave = float(str(text[i+1]))
                    
                if (str(text[i]) == 'NG'):
                    self.NG = int(str(text[i+1]))
                if (str(text[i]) == 'Nk'):
                    self.Nk = int(str(text[i+1]))

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
        
        if (self.NKS > self.NG):
            print('# NOTICE: NKS is larger than the size of the Hilbert space, NG. ')
            print('# NKS is replaced by NG ')
            self.NKS = 1*self.NG

        print('#=====Print the parmeters')
        print('# sys_name =', self.sys_name)
        print('# plot_figure_option =', self.plot_figure_option)
        print('# PC_option =', self.PC_option)
        print('# minimal_output =', self.minimal_output)
        print('# Fortlib_option =', self.Fortlib_option)
        print('# ')
        print('# a =', self.a, ' [a.u.] =',self.a*aB, ' [nm]')
        print('# flat_length =', self.flat_length, ' [a.u.] =',self.flat_length*aB, ' [nm]: IGNORED if negative value')
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
            print('# e*a*E0 =', self.a*self.E0, ' [a.u.] =', self.a*self.E0*Hartree, ' [eV]')
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
                print('# e*a*E0 =', self.a*self.E0[icolor], ' [a.u.] =', self.a*self.E0[icolor]*Hartree, ' [eV]')

    def grid_constructions(self):
        self.b = tpi/self.a
        self.H = self.a/np.float(self.NG)
        self.x = np.linspace(0.0, self.a, num=self.NG, endpoint=False, dtype='float64')
        self.G = np.fft.fftfreq(self.NG)*(self.b*np.float(self.NG))
        #Brillouin zone construction
        self.k = np.linspace(-0.5*self.b, 0.5*self.b, num=self.Nk, endpoint=False, dtype='float64')
        self.k = self.k + (0.5*self.b)/np.float(self.Nk)
        print('# ')
        print('# b = 2pi/a =', self.b, ' [a.u.] =', self.b/aB, '[/nm]')
        print('# H =', self.H, ' [a.u.] =',self.H*aB, ' [nm]')
        print('# (pi/H)^2 =', (pi/self.H)**2, ' [a.u.] =', (pi/self.H)**2*Hartree, ' [eV]')
        print('# Nk*a =', self.Nk*self.a, ' [a.u.] =', (self.Nk*self.a)*aB, ' [nm]')
        print('# ')

    def get_Nocc(self):
        self.Nocc = int(self.Nave/2.0)
        if (abs(self.Nave - 2.0*self.Nocc) > 1.0e-8):
            print('# The average number of particle in a cell: ', self.Nave)
            print('# ERROR: Currenty, metallic occupation is not supported.')
            print('# The average number should be even.')
            sys.exit()
            
