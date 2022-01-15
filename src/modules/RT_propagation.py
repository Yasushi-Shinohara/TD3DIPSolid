# coding: UTF-8
# This is created 2020/11/17 by Y. Shinohara
# This is lastly modified 2020/11/17 by Y. Shinohara
import sys
import os
from modules.constants import tpi, zI
import numpy as np
import ctypes as ct

class RT_propagation_class():
    def __init__(self):
        self.something = None
        self.FL = None
        self.ref_NG = None
        self.ref_Nocc = None
        self.ref_Nk = None
        self.ref_dt = None

    def uGbk_forward(self, propagator_option, Fortlib_option):
        if (propagator_option.lower() == 'exp'):
            uGbk_forward = self.uGbk_forward_exp
            if (Fortlib_option):
                uGbk_forward = self.uGbk_forward_exp_Fortran
            print('# The exponential expression for the temporal propagator is chosen.')
        elif (propagator_option.upper() == 'RK4'):
            uGbk_forward = self.uGbk_forward_RK4
            if (Fortlib_option):
                uGbk_forward = self.uGbk_forward_RK4_Fortran
            print('# The Runge-Kutta 4th for the temporal propagator is chosen.')
        elif ((propagator_option.upper() == 'RK4FFT') or (propagator_option.upper() == 'RK4_FFT')):
            uGbk_forward = self.uGbk_forward_RK4FFT
            if (Fortlib_option):
                uGbk_forward = self.uGbk_forward_RK4FFT_Fortran
            print('# The Runge-Kutta 4th with FFT for the temporal propagator is chosen.')
        elif (propagator_option.upper() == 'KS'):
            uGbk_forward = self.uGbk_forward_KS
            if (Fortlib_option):
                uGbk_forward = self.uGbk_forward_KS_Fortran
        else :
            print('# ERROR: undefined propagator_option is called.')
            sys.exit()
        return uGbk_forward

    def hdt2U(self, hdt):
        eigs, coef = np.linalg.eigh(hdt)
        U = np.exp(-zI*eigs[0])*np.outer(coef[:,0],np.conj(coef[:,0]))
        for ib in range(1,len(eigs)):
            U = U + np.exp(-zI*eigs[ib])*np.outer(coef[:,ib],np.conj(coef[:,ib]))
        return U
        
    def uGbk_forward_exp(self, param, uGbk, hGGk, tGGk, vx):
        for ik in range(param.Nk):
            U = self.hdt2U(hGGk[:,:,ik]*param.dt)
            uGbk[:,:,ik] = np.dot(U, uGbk[:,:,ik])
        return uGbk

    def u_h2hu(self, uGb, hGG):
        return np.dot(hGG,uGb)

    def uGbk_forward_RK4(self, param, uGbk, hGGk, tGGk, vx):
        for ik in range(param.Nk):
            k1 = self.u_h2hu(uGbk[:,:,ik], hGGk[:,:,ik])/zI
            k2 = self.u_h2hu(uGbk[:,:,ik] + 0.5*param.dt*k1, hGGk[:,:,ik])/zI
            k3 = self.u_h2hu(uGbk[:,:,ik] + 0.5*param.dt*k2, hGGk[:,:,ik])/zI
            k4 = self.u_h2hu(uGbk[:,:,ik] + param.dt*k3, hGGk[:,:,ik])/zI
            uGbk[:,:,ik] = uGbk[:,:,ik] + (k1 + 2.0*k2 + 2.0*k3 + k4)*param.dt/6.0 
        return uGbk

    def u_t_v2hu_FFT(self, uGb, tGGdiag, vx):
        NBact = np.shape(uGb)[1]
        vuxb = np.empty_like(uGb)
        tuGb = np.empty_like(uGb)
        uxb = np.fft.ifft(uGb, axis=0)
        for ib in range(NBact):
            vuxb[:,ib] = vx[:]*uxb[:,ib]
            tuGb[:,ib] = tGGdiag[:]*uGb[:,ib]
        vuGb = np.fft.fft(vuxb, axis=0)
        return tuGb+vuGb

    def uGbk_forward_RK4FFT(self, param, uGbk, hGGk, tGGk, vx):
        tGGdiagk = np.diagonal(tGGk, axis1 = 0, axis2 = 1).T
        for ik in range(param.Nk):
            k1 = self.u_t_v2hu_FFT(uGbk[:,:,ik], tGGdiagk[:,ik], vx)/zI
            k2 = self.u_t_v2hu_FFT(uGbk[:,:,ik] + 0.5*param.dt*k1, tGGdiagk[:,ik], vx)/zI
            k3 = self.u_t_v2hu_FFT(uGbk[:,:,ik] + 0.5*param.dt*k2, tGGdiagk[:,ik], vx)/zI
            k4 = self.u_t_v2hu_FFT(uGbk[:,:,ik] + param.dt*k3, tGGdiagk[:,ik], vx)/zI
            uGbk[:,:,ik] = uGbk[:,:,ik] + (k1 + 2.0*k2 + 2.0*k3 + k4)*param.dt/6.0 
        return uGbk

    def uGbk_forward_KS(self, param, uGbk, hGGk, tGGk, vx):
        tGGdiagk = np.diagonal(tGGk, axis1 = 0, axis2 = 1).T
        O = np.diag(param.G)
        Odiag = 1.0*param.G
        for ik in range(param.Nk):
            #U = self.u_hdt2U_KS(uGbk[:,:,ik], hGGk[:,:,ik]*param.dt, NKS = param.NKS)                    #An option
            #U = self.u_tt_vt2U_KS(uGbk[:,:,ik], tGGdiagk[:,ik]*param.dt, vx*param.dt, NKS = param.NKS)   #An option
            U = self.u_hdt_Odiag2U_KS(uGbk[:,:,ik], hGGk[:,:,ik]*param.dt, Odiag, NKS = param.NKS)         #An option
            #U = self.u_hdt_O2U_KS(uGbk[:,:,ik], hGGk[:,:,ik]*param.dt, O, NKS = param.NKS)         #An option
            uGbk[:,:,ik] = np.dot(U, uGbk[:,:,ik])
        return uGbk

    def u_hdt_O2U_KS(self, u, hdt, O, NKS = 4):
        k = 1.0*u
        ktemp = 1.0*u
        Nocc = np.shape(u)[1]
        nKS = np.shape(ktemp)[1]
        while (nKS < NKS):  #Making Krylov subspace (KS) up to a dimension of NKS via G-matrix.
            ktemp = np.dot(O,ktemp)
            k = np.hstack((k, ktemp))
            nKS = np.shape(k)[1]
        q, r = np.linalg.qr(k)

        hdtr = np.dot(np.conj(q).T, self.u_h2hu(q, hdt))
        eigs, coef = np.linalg.eigh(hdtr)
        Ur = np.exp(-zI*eigs[0])*np.outer(coef[:,0],np.conj(coef[:,0]))  #Constructing a unitary operator in the KS
        for i in range(1,len(eigs)):
            Ur = Ur + np.exp(-zI*eigs[i])*np.outer(coef[:,i],np.conj(coef[:,i]))

        U = np.dot(q, np.dot(Ur , np.conj(q).T)) #Converting the subspace representation to the original space
        return U

    def u_tt_vt2U_KS(self, u, tdiagt, vt, NKS = 4):
        k = 1.0*u
        ktemp = 1.0*u
        nKS = np.shape(ktemp)[1]
        while (nKS < NKS):
            ktemp = self.u_t_v2hu_FFT(ktemp, tdiagt, vt)
            k = np.hstack((k, ktemp))
            nKS = np.shape(k)[1]
        q, r = np.linalg.qr(k)

        hdtr = np.dot(np.conj(q).T, self.u_t_v2hu_FFT(q, tdiagt, vt))
        eigs, coef = np.linalg.eigh(hdtr)
        Ur = np.exp(-zI*eigs[0])*np.outer(coef[:,0],np.conj(coef[:,0]))
        for i in range(1,len(eigs)):
            Ur = Ur + np.exp(-zI*eigs[i])*np.outer(coef[:,i],np.conj(coef[:,i]))

        U = np.dot(q, np.dot(Ur , np.conj(q).T))
        return U

 
    def u_hdt_Odiag2U_KS(self, u, hdt, Odiag, NKS = 4):
        k = 1.0*u
        ktemp = 1.0*u
        Nocc = np.shape(u)[1]
        nKS = np.shape(ktemp)[1]
        while (nKS < NKS):  #Making Krylov subspace (KS) up to a dimension of NKS via G-matrix.
            for ib in range(Nocc):
                ktemp[:,ib] = Odiag[:]*ktemp[:,ib]
            k = np.hstack((k, ktemp))
            nKS = np.shape(k)[1]
        q, r = np.linalg.qr(k)

        hdtr = np.dot(np.conj(q).T, self.u_h2hu(q, hdt))
        eigs, coef = np.linalg.eigh(hdtr)
        Ur = np.exp(-zI*eigs[0])*np.outer(coef[:,0],np.conj(coef[:,0]))  #Constructing a unitary operator in the KS
        for i in range(1,len(eigs)):
            Ur = Ur + np.exp(-zI*eigs[i])*np.outer(coef[:,i],np.conj(coef[:,i]))

        U = np.dot(q, np.dot(Ur , np.conj(q).T)) #Converting the subspace representation to the original space
        return U

    def u_hdt2U_KS(self, u, hdt, NKS = 4):
        k = 1.0*u
        ktemp = 1.0*u
        nKS = np.shape(ktemp)[1]
        while (nKS < NKS):
            ktemp = self.u_h2hu(ktemp, hdt)
            k = np.hstack((k, ktemp))
            nKS = np.shape(k)[1]
        q, r = np.linalg.qr(k)

        hdtr = np.dot(np.conj(q).T, self.u_h2hu(q, hdt))
        eigs, coef = np.linalg.eigh(hdtr)
        Ur = np.exp(-zI*eigs[0])*np.outer(coef[:,0],np.conj(coef[:,0]))
        for i in range(1,len(eigs)):
            Ur = Ur + np.exp(-zI*eigs[i])*np.outer(coef[:,i],np.conj(coef[:,i]))

        #NG = np.shape(u)[0]
        #U = 0.0*hdt
        #for iG in range(NG):
        #    for jG in range(NG):
        #        for i in range(len(eigs)):
        #            for j in range(len(eigs)):
        #                U[iG, jG] = U[iG, jG] + q[iG,i]*Ur[i,j]*np.conj(q[jG,j])
        #for iG in range(NG):
        #    for jG in range(NG):
        #        U[iG, jG] = np.dot(q[iG,:], np.dot(Ur, np.conj(q[jG,:])))

        U = np.dot(q, np.dot(Ur , np.conj(q).T))
        return U

    def Prep4Fortlib(self, param):
#        self.ref_NG   = ct.byref(ct.c_int32(param.NG)  )
        self.ref_Nocc = ct.byref(ct.c_int32(param.Nocc))
#        self.ref_Nk   = ct.byref(ct.c_int32(param.Nk)  )
        self.ref_dt   = ct.byref(ct.c_double(param.dt) )
        dir_name = os.path.dirname(os.path.abspath(__file__)).strip('modules')
        print('# Fortlib.so: ',dir_name+"Fortlib.so")
        self.FL = np.ctypeslib.load_library(dir_name+"Fortlib.so",".")
        self.FL.ugbk_forward_rk4_.argtypes = [
            np.ctypeslib.ndpointer(dtype='complex128'),  #ubk
            np.ctypeslib.ndpointer(dtype='complex128'),  #hGGk
#            ct.POINTER(ct.c_int32),                      #NG
            np.ctypeslib.ndpointer(dtype='int32'),       #NG
            ct.POINTER(ct.c_int32),                      #Nocc
#            ct.POINTER(ct.c_int32),                      #Nk
            np.ctypeslib.ndpointer(dtype='int32'),       #Nk
            ct.POINTER(ct.c_double),]                    #dt
        self.FL.ugbk_forward_rk4_.restype = ct.c_void_p
        self.FL.ugbk_forward_rk4fft_.argtypes = [
            np.ctypeslib.ndpointer(dtype='complex128'),  #ubk
            np.ctypeslib.ndpointer(dtype='complex128'),  #tGGk
            np.ctypeslib.ndpointer(dtype='float64'),     #vx
#            ct.POINTER(ct.c_int32),                      #NG
            np.ctypeslib.ndpointer(dtype='int32'),       #NG
            ct.POINTER(ct.c_int32),                      #Nocc
#            ct.POINTER(ct.c_int32),                      #Nk
            np.ctypeslib.ndpointer(dtype='int32'),       #Nk
            ct.POINTER(ct.c_double),]                    #dt
        self.FL.ugbk_forward_rk4fft_.restype = ct.c_void_p
        self.FL.ugbk_forward_exp_.argtypes = [
            np.ctypeslib.ndpointer(dtype='complex128'),  #ubk
            np.ctypeslib.ndpointer(dtype='complex128'),  #hGGk
#            ct.POINTER(ct.c_int32),                      #NG
            np.ctypeslib.ndpointer(dtype='int32'),       #NG
            ct.POINTER(ct.c_int32),                      #Nocc
#            ct.POINTER(ct.c_int32),                      #Nk
            np.ctypeslib.ndpointer(dtype='int32'),       #Nk
            ct.POINTER(ct.c_double),]                    #dt
        self.FL.ugbk_forward_exp_.restype = ct.c_void_p
        #self.FL.ugbk_forward_ks_.argtypes = [
        #    np.ctypeslib.ndpointer(dtype='complex128'),  #ubk
        #    np.ctypeslib.ndpointer(dtype='complex128'),  #tGGk
        #    np.ctypeslib.ndpointer(dtype='float64'),     #vx
        #    np.ctypeslib.ndpointer(dtype='complex128'),  #O
        #    ct.POINTER(ct.c_int32),                      #NG
        #    ct.POINTER(ct.c_int32),                      #Nocc
        #    ct.POINTER(ct.c_int32),                      #Nk
        #    ct.POINTER(ct.c_double),                     #dt
        #    ct.POINTER(ct.c_int32),]                     #NKS
        #self.FL.ugbk_forward_ks_.restype = ct.c_void_p
        self.FL.writeout_ompinfo_.argtypes = [
            ct.POINTER(ct.c_int32),]                     #Nk
        self.FL.writeout_ompinfo_.restype = ct.c_void_p
        self.FL.writeout_ompinfo_(self.ref_Nk)

    def uGbk_forward_RK4_Fortran(self, param, uGbk, hGGk, tGGk, vx):
        self.FL.ugbk_forward_rk4_(uGbk, hGGk, self.ref_NG, self.ref_Nocc, self.ref_Nk, self.ref_dt)
        return uGbk

    def uGbk_forward_RK4FFT_Fortran(self, param, uGbk, hGGk, tGGk, vx):
        self.FL.ugbk_forward_rk4fft_(uGbk, tGGk, vx, self.ref_NG, self.ref_Nocc, self.ref_Nk, self.ref_dt)
        return uGbk

    def uGbk_forward_exp_Fortran(self, param, uGbk, hGGk, tGGk, vx):
        self.FL.ugbk_forward_exp_(uGbk, hGGk, self.ref_NG, self.ref_Nocc, self.ref_Nk, self.ref_dt)
        return uGbk

    #def uGbk_forward_KS_Fortran(self, param, uGbk, hGGk, tGGk, vx):
        #O = np.diag(param.G + 0.0j)
        #self.FL.ugbk_forward_ks_(uGbk, tGGk, vx, O, self.ref_NG, self.ref_Nocc, self.ref_Nk, self.ref_dt, self.ref_NKS)
        #return uGbk
