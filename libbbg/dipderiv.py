# --------------------------------------------------------------- #
#                       DIPOLE DERIVATIVES                        #
# --------------------------------------------------------------- #

__all__=["DipDeriv"]

import os
from numpy import zeros, float64, array, transpose, dot
from .gaussfreq import FREQ
from .units import UNITS, Atom

class DipDeriv(FREQ):
    """dipole derivatives wrt normal mode class in helico representation"""
    
    def __init__(self,file="",L=0,step=0.006):
        FREQ.__init__(self,file)
        # step of differentiation
        self.h = step
        # derive first derivatives (helico!)
        self.fdip = self.FDeriv(Print=0,Debye=0,divide=0)
        # calculate second derivatives (helico)
        self.sdip = self.SDeriv()
        
    def SDeriv(self,symm=0):
        """ computes second derivatives of dipole moment wrt normal modes """
        A = os.listdir('.')
        A.sort()
        B = A[:]
        for file in B:
            if (not file.endswith('_.log') ) : A.remove(file)

        a = len(A)
        b = self.Nmodes ; N = self.Natoms
        D = zeros((a  ,3*N,3  ))  # first derivatives wrt cartesian coordinate
        for i in range(a):
            D[i] = self.DipoleDeriv(A[i])

        # make second and third derivatives wrt cartesian coordinates
        S = zeros((3*N,3*N,3  ))   # second ij      # tensor 9x9x3 dla wody
        for i in range(3*N):
            K = 4*i + 1
            S[i] = (1./12.) * ( 8.*(D[K+1] - D[K+2]) + 
                                   (D[K+3] - D[K+0]) )  / (self.h * self.AngstromToBohr )
        LT = transpose(self.L)

        # transform to normal mode space
        XS = dot(dot(LT,S[:,:,0]),self.L)
        YS = dot(dot(LT,S[:,:,1]),self.L)
        ZS = dot(dot(LT,S[:,:,2]),self.L)

        E = tensordot(S,self.L,(1,0))
        E = tensordot(LT,S,(1,0))
        if symm:        
           for M in [XS,YS,ZS]:
               for i in range(len(XS)):
                   for j in range(len(XS)):
                       if i>j:
                          M[i][j] = 0.5 * ( M[i][j] + M[j][i] )
                          M[j][i] = M[i][j]

        s_diag = zeros((b,3),dtype=float64)
        for i in range(b):
            s_diag[i,0] = XS[i,i]
            s_diag[i,1] = YS[i,i]
            s_diag[i,2] = ZS[i,i]
        
            
        return XS,YS,ZS,s_diag     
       
    
    
    def __repr__(self):
        log = "\n"
        return str(log)
    
