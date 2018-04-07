# ------------------------------------------ #
#                UTILITIES                   #
# ------------------------------------------ #

__all__=['SVDSuperimposer','ParseDMA','RotationMatrix',
         'Periodic','Atomn','DMAMadrixMultiply',
         'PRINT','PRINTV','PRINTL','ParseDmatFromFchk',
         'Emtp','FrequencyShift','func_to_method',
         'Trapezoidal','Simpson','GaussLegendre2','integrate',
         'Read_xyz_file','Vr_dma','Allign','ElectricField',
         'FrequencyShiftPol','Newton','Parse_EDS_InteractionEnergies',
         'CalcStep','ModifyStruct','ParseUnitedAtoms',
         'MakeSoluteAndSolventFiles','GROUPS','DistanceRelationMatrix',
         'status','ROTATE','get_tcf','choose','get_pmloca',
         'ParseVecFromFchk','interchange','Peak','PUPA','VIB',
         'ParseFCFromFchk','ParseDipoleDerivFromFchk',
         'ParseFockFromGamessLog','lind','order','check_sim','MakeMol',
         'ParseDistributedPolarizabilitiesFromGamessEfpFile','reorder',
         'ParseEFPInteractionEnergies','secant','RungeKutta',
         'numerov1','numerov2','simpson','simpson_nonuniform','fder5pt',
         'ParseDMAFromGamessEfpFile','dihedral','Peak2DIR',
         'text_to_list','QMFile','Emtp_charges','MDOut',
         'ParseLmocFromGamessEfpFile','resol','ft_1d','FF_scheme','diff',
         'calc_tcf','autocorr','crosscorr','ParseEnergyFromFchk',
         'make_bqc_inp','bcolors','ParseDipoleMomentFromFchk',
         'ParseGradFromFchk','distribute','ParseAlphaOrbitalEnergiesFromFchk','TIMER',
         'ParseElectronsFromFchk','ParseDistributedPolarizabilitiesWrtImFreqFromGamessEfpFile',
         'ParseChargesFromFchk','glue','UnitaryOptimizer','UnitaryOptimizer_4_2'] #'gen_camm'
         
__version__ = '3.3.4'

import os
import re
import qm
import copy
import time
import string
import operator
import math
import numpy
import numpy.linalg
import scipy.optimize
import scipy.interpolate
import PyQuante
import dma
import units
import letters
import fourier
import re_templates

import_matplotlib = int(os.environ['LIBBBG_IMPORT_MATPLOTLIB'])

if import_matplotlib:
   import matplotlib.font_manager
   import pylab
   import scitools.numpyutils
   __all__+= ['QMOscillator','circles','PotentialContourMap','Grid2D']

   #import re, qm, PyQuante,  \
   #       scipy.optimize, scipy.integrate, numpy,\
   #       math, numpy.linalg, dma, units, re_templates,\
   #       copy, os, math, matplotlib.font_manager,\
   #       pylab, scitools.numpyutils, scipy.interpolate,\
   #       letters, fourier, string, time #, coulomb.multip
#else:
#   import re, qm, PyQuante,  \
#          scipy.optimize, scipy.integrate, numpy,\
#          math, numpy.linalg, dma, units, re_templates,\
#          copy, os, math,\
#          scipy.interpolate,\
#          letters, fourier, string, time #, coulomb.multip

uAtom = units.Atom
uUNITS= units.UNITS

class UnitaryOptimizer_4_2(object):
  """
 ---------------------------------------------------------------------------------

 Finds the unitary matrix X that optimizes the following function:
 
 Z(X) = \sum_{ijklmn} X_{ki} X_{lj} X_{mi} X_{nj} R_{ijklmn} 
      + \sum_{ijk}    X_{ji} X_{ki}               P_{ijk}
 
 where 
   * X is a square unitary matrix of size N x N
   * R is a general real 6-th rank tensor of size N^6
   * P is a general real 3-rd rank tensor of size N^3
 
 Usage:
   optimizer = UnitaryOptimizer(R, P, conv=1.0e-8, maxiter=100, verbose=True)
   optimizer.maximize() #or minimize()
   X = optimizer.X
   Z = optimizer.Z
 
 ---------------------------------------------------------------------------------
                                                       Last Revision: 07.04.2018
 """
  def __init__(self, R, P, conv=1.0e-8, maxiter=100, verbose=True):
      """Initialize with R and P matrix, as well as optimization options"""
      self._R   = R.copy()
      self._P   = P.copy()
      self._R0  = R.copy()
      self._P0  = P.copy()
      self.X    = None
      self._d   = len(P)
      # optimization options
      self.conv   = conv
      self.maxiter=maxiter
      self.verbose=verbose
      # initial Z value
      self._Zinit = self._eval_Z(numpy.identity(self._d), self._R0, self._P0)

  def maximize(self):
      """Maximize the Z function under unitary constraint for X""" 
      self.run('max')

  def minimize(self):
      """Minimize the Z function under unitary constraint for X"""
      self.run('min')

  def run(self, opt='minimize'):
      """Perform the optimization"""
      assert (opt.lower().startswith('min') or opt.lower().startswith('max')), 'Unrecognized optimization mode < %s >' % opt
      self._refresh()
      self._run(opt.lower())

  @property
  def Z(self):
      """Return the current value of objective function"""
      z = self._eval_Z(self.X, self._R0, self._P0)
      return z

  # -- protected
  def _run(self, opt):
      """Perform the optimization (protected interface)"""
      conv = 1e8
      #Xold = numpy.identity(self._d)
      Zold = self._Zinit
      Xacc = numpy.identity(self._d)
      success = False
      if self.verbose: 
         print     " Start  : Z[1] = %15.6f" % Zold
      niter = 0
      while (conv > self.conv):
         i, j, gamma = self._find_next(opt)
         #print " Chosen x = %14.4f" % gamma
         Xnew = self._form_X(i, j, gamma)
         self._update_RP(Xnew)
         Znew = self._eval_Z(numpy.identity(self._d), self._R, self._P)
         conv = abs(Znew-Zold)
         #Xold = Xnew.copy()
         Zold = Znew
         niter += 1
         Xacc = numpy.dot(Xacc, Xnew)
         if self.verbose:
            print  " Iter %2d: Z[X] = %15.6f  Conv= %15.6f" % (niter, Znew, conv)
         if niter > self.maxiter: 
            print " Optimization unsuccesfull! Maximum iteration number exceeded!"
            success = False
            break
      success = True if niter <= self.maxiter else False
      self.X = Xacc.copy()
      if (self.verbose and success):
         print " Optimization succesfull!\n"
         print " Optimized Z[X] value: %15.6f" % self.Z

  def _update_RP(self, X):
      """Update R and P tensors by transforming them by using X matrix"""
      P = numpy.tensordot(self._P, X, (1,0))
      P = P.transpose(0,2,1)
      P = numpy.tensordot(P, X, (2,0))

      R = numpy.tensordot(self._R, X, (2,0)) 
      R = R.transpose(0,1,5,2,3,4)
      R = numpy.tensordot(R  , X, (3,0))
      R = R.transpose(0,1,2,5,3,4)
      R = numpy.tensordot(R  , X, (4,0))
      R = R.transpose(0,1,2,3,5,4)
      R = numpy.tensordot(R  , X, (5,0))

      #P = self._P.copy()
      #R = self._R.copy()
      #P.fill(0); R.fill(0)
      #N = self._d

      #for i in range(N):
      #    for J in range(N): 
      #        for K in range(N):
      #            v = 0.0
      #            for j in range(N):
      #                for k in range(N):
      #                    v += X[j,J] * X[k,K] * self._P[i,j,k]
      #            P[i,J,K] = v


      #for i in range(N):
      #    for j in range(N):
      #        for K in range(N):
      #            for L in range(N):
      #                for M in range(N):
      #                    for NN in range(N):
      #                        v = 0.0
      #                        for k in range(N):
      #                            for l in range(N):
      #                                for m in range(N):
      #                                    for n in range(N):
      #                                        v += X[k,K] * X[l,L] * X[m,M] * X[n,NN] * self._R[i,j,k,l,m,n]
      #                        R[i,j,K,L,M,NN] = v
      self._P = P.copy()
      self._R = R.copy()

  def _refresh(self):
      """Restore the initial state of the optimizer"""
      self._R = self._R0.copy()
      self._P = self._P0.copy()
      self.X  = None

  def _find_next(self, opt):
      """Determine next pair of degrees of freedom for 2D rotation"""
      optfunc = operator.lt if opt.startswith('min') else operator.gt
      I, J = 0, 1
      Gamma = None
      dZold = 1e8 if opt.startswith('min') else -1e8
      for j in range(self._d):
          for i in range(j):
              a0, a1, a2, a3, a4, b1, b2, b3, b4 = self._get_Fourier_coeffs(i, j)
              gamma   = self._find_x(a0, a1, a2, a3, a4, b1, b2, b3, b4, i, j, opt)
              dZ = self._eval_dZ(gamma, self._P, self._R, i, j)
              if optfunc(dZ, dZold):
                 Gamma = gamma
                 I = i
                 J = j
                 dZold = dZ
      return I, J, Gamma

  def _find_x(self, a0, a1, a2, a3, a4, b1, b2, b3, b4, i, j, opt):
      """Find the optimal 2D rotation angle"""
      # Boyd's method in 4 dimensions: Boyd, J.P.; J. Eng. Math. (2006) 56:203-219
      d = a4 - 1.0j * b4
      K = numpy.zeros((8,8), numpy.complex64)
      K[0,1] = 1.0
      K[1,2] = 1.0
      K[2,3] = 1.0
      K[3,4] = 1.0
      K[4,5] = 1.0
      K[5,6] = 1.0
      K[6,7] = 1.0
      K[7,0] =-(a4 + 1.0j * b4) / d
      K[7,1] =-(a3 + 1.0j * b3) / d
      K[7,2] =-(a2 + 1.0j * b2) / d
      K[7,3] =-(a1 + 1.0j * b1) / d
      K[7,4] =- a0 * 2.0        / d
      K[7,5] =-(a1 - 1.0j * b1) / d
      K[7,6] =-(a2 - 1.0j * b2) / d
      K[7,7] =-(a3 - 1.0j * b3) / d
      #
      E, X = numpy.linalg.eig(K)
      X   = -1.0j * numpy.log(E)
      #print "Imaginary part of X: "
      #libbbg.utilities.PRINT(X.imag)
      X = X.real
      X[numpy.where(X<0.0)] += numpy.pi
      #print "Real      part of X: "
      #libbbg.utilities.PRINT(X)

      # Find optimal gamma 
      gamma = None
      if opt.startswith('min'):
         Zold = 1.0e8
         for x in X:
             Z = self._eval_Z(self._form_X(i, j, x), self._R, self._P)
             if Z < Zold:
                gamma = x
                Zold = Z
      else:
         Zold = -1e8
         for x in X:
             Z = self._eval_Z(self._form_X(i, j, x), self._R, self._P)
             if Z > Zold:
                gamma = x
                Zold = Z
      assert gamma is not None, "Error while searching for optimum!"
      return gamma

  def _get_Fourier_coeffs(self, I, J):
      """Retrieve ABCD parameters for root search"""
      a0, a1, a2, a3, a4, b1, b2, b3, b4 = 0, 0, 0, 0, 0, 0, 0, 0, 0
      d  = lambda i, j: 0.0 if (i!=j) else 1.0
      a_ = lambda i, k: (-d(I,i)*d(J,k)+d(I,k)*d(J,i))*(1.0-d(i,k))
      b_ = lambda i, k: d(i,k) * (d(I,k)+d(J,i))
      c_ = lambda i, k: d(i,k) * (1.0-d(I,k))*(1.0-d(J,i))
      N = self._d

      # P-contribution
      for i in range(N):
          for j in range(N):
              for k in range(N):
                  p   = self._P[i,j,k]
                  A =   (d(I,j)*d(J,i) - d(I,i)*d(J,j)) * (1-d(i,j))
                  B =   -d(i,j)*(d(I,j)+d(J,i))
                  C =    d(i,k)*(1-d(I,k))*(1-d(J,i))
                  D =    d(i,k)*(d(I,k)+d(J,i))
                  E =   -d(I,i)*d(J,k)*(1-d(i,k))
                  F =    d(I,k)*d(J,i)*(1-d(i,k))
                  G =    (1-d(i,k))*(d(I,k)*d(J,i) - d(I,i)*d(J,k))
                  H =   -d(i,k)*(d(I,k)+d(J,i))
                  II=    d(i,j)*(1-d(I,j))*(1-d(J,i))
                  JJ=    d(i,j)*(d(I,j)+d(J,i))
                  K =   -d(I,i)*d(J,j)*(1-d(i,j))
                  L =    d(I,j)*d(J,i)*(1-d(i,j))
                  a0 += p * (A*D+G*JJ+B*(E+F)+H*(K+L)) / 2.0
                  a1 += p * (A*C+G*II)
                  a2 += p * (A*D+G*JJ-B*(E+F)-H*(K+L)) / 2.0
                  b1 += p * (B*C+H*II)
                  b2 += p * (H*JJ+B*D+G*(K+L)+A*(E+F)) / 2.0

      # R-contribution
      for i in range(N):
          for j in range(N):
              for k in range(N):
                  for l in range(N):
                      for m in range(N):
                          for n in range(N):
                              r = self._R[i,j,k,l,m,n] 
                              # First R batch
                              A = a_(i,k)
                              B =-b_(i,k)
                              C = a_(i,m)
                              D = c_(i,m)
                              E = b_(i,m)
                              F = a_(j,l)
                              G = c_(j,l)
                              H = b_(j,l)
                              II= a_(j,n)
                              JJ= c_(j,n)
                              K = b_(j,n)
                              a0 += r * (A*C*F*K+A*C*H*II+4*A*D*G*K+4*A*D*H*JJ+A*E*F*II+4*A*E*G*JJ+3*A*E*H*K\
                                        +3*B*C*F*II+4*B*C*G*JJ+B*C*H*K+4*B*D*F*JJ+4*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0
                              a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4*A*D*G*JJ+3*A*D*H*K+3*A*E*G*K+3*A*E*H*JJ\
                                        +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0
                              a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0
                              a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ\
                                         -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0
                              a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0
                              b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3*B*C*F*JJ+3*B*C*G*II\
                                         +3*B*D*F*II+4*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b2 += r * (A*C*F*II+2*A*C*G*JJ+A*C*H*K+2*A*D*F*JJ+2*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K\
                                        +B*C*H*II+2*B*D*G*K+2*B*D*H*JJ+B*E*F*II+2*B*E*G*JJ+B*E*H*K) / 4.0
                              b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II\
                                        -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0
                              # Second R batch
                              A = a_(i,m)
                              B =-b_(i,m)
                              C = a_(i,k)
                              D = c_(i,k)
                              E = b_(i,k)
                              F = a_(j,l)
                              G = c_(j,l)
                              H = b_(j,l)
                              II= a_(j,n)
                              JJ= c_(j,n)
                              K = b_(j,n)
                              a0 += r * (A*C*F*K+A*C*H*II+4*A*D*G*K+4*A*D*H*JJ+A*E*F*II+4*A*E*G*JJ+3*A*E*H*K\
                                        +3*B*C*F*II+4*B*C*G*JJ+B*C*H*K+4*B*D*F*JJ+4*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0
                              a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4*A*D*G*JJ+3*A*D*H*K+3*A*E*G*K+3*A*E*H*JJ\
                                        +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0
                              a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0
                              a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ\
                                         -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0
                              a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0
                              b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3*B*C*F*JJ+3*B*C*G*II\
                                         +3*B*D*F*II+4*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b2 += r * (A*C*F*II+2*A*C*G*JJ+A*C*H*K+2*A*D*F*JJ+2*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K\
                                        +B*C*H*II+2*B*D*G*K+2*B*D*H*JJ+B*E*F*II+2*B*E*G*JJ+B*E*H*K) / 4.0
                              b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II\
                                        -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0
                              # Third R batch
                              A = a_(j,l)
                              B =-b_(j,l)
                              C = a_(i,k)
                              D = c_(i,k)
                              E = b_(i,k)
                              F = a_(i,m)
                              G = c_(i,m)
                              H = b_(i,m)
                              II= a_(j,n)
                              JJ= c_(j,n)
                              K = b_(j,n)
                              a0 += r * (A*C*F*K+A*C*H*II+4*A*D*G*K+4*A*D*H*JJ+A*E*F*II+4*A*E*G*JJ+3*A*E*H*K\
                                        +3*B*C*F*II+4*B*C*G*JJ+B*C*H*K+4*B*D*F*JJ+4*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0
                              a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4*A*D*G*JJ+3*A*D*H*K+3*A*E*G*K+3*A*E*H*JJ\
                                        +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0
                              a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0
                              a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ\
                                         -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0
                              a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0
                              b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3*B*C*F*JJ+3*B*C*G*II\
                                         +3*B*D*F*II+4*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b2 += r * (A*C*F*II+2*A*C*G*JJ+A*C*H*K+2*A*D*F*JJ+2*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K\
                                        +B*C*H*II+2*B*D*G*K+2*B*D*H*JJ+B*E*F*II+2*B*E*G*JJ+B*E*H*K) / 4.0
                              b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II\
                                        -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0
                              # Fourth R batch
                              A = a_(j,n)
                              B =-b_(j,n)
                              C = a_(i,k)
                              D = c_(i,k)
                              E = b_(i,k)
                              F = a_(i,m)
                              G = c_(i,m)
                              H = b_(i,m)
                              II= a_(j,l)
                              JJ= c_(j,l)
                              K = b_(j,l)
                              a0 += r * (A*C*F*K+A*C*H*II+4*A*D*G*K+4*A*D*H*JJ+A*E*F*II+4*A*E*G*JJ+3*A*E*H*K\
                                        +3*B*C*F*II+4*B*C*G*JJ+B*C*H*K+4*B*D*F*JJ+4*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0
                              a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4*A*D*G*JJ+3*A*D*H*K+3*A*E*G*K+3*A*E*H*JJ\
                                        +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0
                              a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0
                              a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ\
                                         -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0
                              a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0
                              b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3*B*C*F*JJ+3*B*C*G*II\
                                         +3*B*D*F*II+4*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b2 += r * (A*C*F*II+2*A*C*G*JJ+A*C*H*K+2*A*D*F*JJ+2*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K\
                                        +B*C*H*II+2*B*D*G*K+2*B*D*H*JJ+B*E*F*II+2*B*E*G*JJ+B*E*H*K) / 4.0
                              b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II\
                                        -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                              b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0
      
      return a0, a1, a2, a3, a4, b1, b2, b3, b4

  def _eval_Z(self, X, R, P):
      """Evaluate the objective Z function"""            
      Z = 0.0
      N = self._d
      p = numpy.tensordot(P, X, (1,0))
      p = numpy.tensordot(p, X, (1,0))
      Z = p.diagonal().diagonal().sum()
     
      r = numpy.tensordot(R, X, (2,0)) 
      r = numpy.tensordot(r, X, (2,0)) 
      r = numpy.tensordot(r, X, (2,0)) 
      r = numpy.tensordot(r, X, (2,0)) 
      for i in range(N): 
          for j in range(N): 
              Z+=r[i,j,i,j,i,j]

      #for i in range(N):
      #    for j in range(N):
      #        for k in range(N):
      #            Z += P[i,j,k] * X[j,i] * X[k,i]
      #            for l in range(N):
      #                for m in range(N):
      #                    for n in range(N):
      #                        pass
      #                        Z += R[i,j,k,l,m,n] * X[k,i] * X[l,j] * X[m,i] * X[n,j]
      
      return Z

  def _eval_dZ(self, g, P, R, I, J):
      """Compute the change in Z"""
      X = self._form_X(I, J, g)
      E = numpy.identity(self._d)
      Z_new = self._eval_Z(X, R, P)
      Z_old = self._eval_Z(E, R, P)
      dZ = Z_new - Z_old
      return dZ

  def _form_X(self, i, j, gamma):
      """Form unitary matrix X"""                                        
      X = numpy.identity(self._d)
      T = numpy.cos(gamma)
      g = numpy.sin(gamma)
      X[i,i] = T
      X[j,j] = T
      X[i,j] = g
      X[j,i] =-g
      return X


class UnitaryOptimizer(object):
  """
 ---------------------------------------------------------------------------------

 Finds the unitary matrix X that optimizes the following function:
 
 Z(X) = \sum_{ijkl} X_{ij}X_{kl} R_{jl} - \sum_{ij} X_{ij}P_{j} 
 
 where 
   * X is a square unitary matrix of size N x N
   * R is a square, in general non-symmetric matrix of size N x N
   * P is a vector of length N   
 
 Usage:
   optimizer = UnitaryOptimizer(R, P, conv=1.0e-8, maxiter=100, verbose=True)
   optimizer.maximize() #or minimize()
   X = optimizer.X
   Z = optimizer.Z
 
 ---------------------------------------------------------------------------------
                                                       Last Revision: 25.03.2018
 """
  def __init__(self, R, P, conv=1.0e-8, maxiter=100, verbose=True):
      """Initialize with R and P matrix, as well as optimization options"""
      self._R   = R.copy()
      self._P   = P.copy()
      self._R0  = R.copy()
      self._P0  = P.copy()
      self.X    = None
      self._d   = P.size
      # optimization options
      self.conv   = conv
      self.maxiter=maxiter
      self.verbose=verbose
      # initial Z value
      self._Zinit = self._eval_Z(numpy.identity(self._d), self._R0, self._P0)
      # functions do find roots in 2D unitary optimization step
      self._f = lambda x, A, B, C, D:   A*numpy.sin(x)+B*numpy.cos(x)+  C*numpy.sin(2*x)+  D*numpy.cos(2*x)
      self._fg= lambda x, A, B, C, D:   A*numpy.cos(x)-B*numpy.sin(x)+2*C*numpy.cos(2*x)-2*D*numpy.sin(2*x)
      self._fh= lambda x, A, B, C, D: -(A*numpy.sin(x)+B*numpy.cos(x)+4*C*numpy.sin(2*x)+4*D*numpy.cos(2*x))

  def maximize(self):
      """Maximize the Z function under unitary constraint for X""" 
      self.run('max')

  def minimize(self):
      """Minimize the Z function under unitary constraint for X"""
      self.run('min')

  def run(self, opt='minimize'):
      """Perform the optimization"""
      assert (opt.lower().startswith('min') or opt.lower().startswith('max')), 'Unrecognized optimization mode < %s >' % opt
      self._refresh()
      self._run(opt.lower())

  @property
  def Z(self):
      """Return the current value of objective function"""
      z = self._eval_Z(self.X, self._R0, self._P0)
      return z

  # -- protected
  def _run(self, opt):
      """Perform the optimization (protected interface)"""
      conv = 1e8
      Xold = numpy.identity(self._d)
      Zold = self._Zinit
      Xacc = numpy.identity(self._d)
      success = False
      if self.verbose: 
         print     " Start  : Z[1] = %15.6f" % Zold
      niter = 0
      while (conv > self.conv):
         i, j, gamma = self._find_next(opt)
         Xnew = self._form_X(i, j, gamma)
         self._uptade_RP(Xnew)
         Znew = self._eval_Z(numpy.identity(self._d), self._R, self._P)
         conv = abs(Znew-Zold)
         Xold = Xnew.copy()
         Zold = Znew
         niter += 1
         Xacc = numpy.dot(Xnew, Xacc)
         if self.verbose:
            print  " Iter %2d: Z[X] = %15.6f  Conv= %15.6f" % (niter, Znew, conv)
         if niter > self.maxiter: 
            print " Optimization unsuccesfull! Maximum iteration number exceeded!"
            success = False
            break
      success = True if niter <= self.maxiter else False
      self.X = Xacc
      if (self.verbose and success):
         print " Optimization succesfull!\n"
         print " Optimized Z[X] value: %15.6f" % self.Z

  def _uptade_RP(self, X):
      self._P = numpy.dot(X, self._P)
      self._R = numpy.dot(X, numpy.dot(self._R, X.T))

  def _refresh(self):
      """Restore the initial state of the optimizer"""
      self._R = self._R0.copy()
      self._P = self._P0.copy()
      self.X  = None

  def _find_next(self, opt):
      """Determine next pair of degrees of freedom for 2D rotation"""
      optfunc = operator.lt if opt.startswith('min') else operator.gt
      I, J = 0, 1
      Gamma = 0.0
      dZold = 1e8 if opt.startswith('min') else -1e8
      for j in range(self._d):
          for i in range(j):
              A,B,C,D = self._get_ABCD(i, j)
              gamma   = self._find_x(A, B, C, D, i, j, opt)
              dZ = self._eval_dZ(gamma, self._P, self._R, i, j)
              if optfunc(dZ, dZold):
                 Gamma = gamma
                 I = i
                 J = j
                 dZold = dZ
      return I, J, Gamma

  def _find_x(self, A, B, C, D, i, j, opt):
      """Find the optimal 2D rotation angle"""
      #f = lambda x, A, B, C, D:   A*numpy.sin(x)+B*numpy.cos(x)+  C*numpy.sin(2*x)+  D*numpy.cos(2*x)
      #fg= lambda x, A, B, C, D:   A*numpy.cos(x)-B*numpy.sin(x)+2*C*numpy.cos(2*x)-2*D*numpy.sin(2*x)
      #fh= lambda x, A, B, C, D: -(A*numpy.sin(x)+B*numpy.cos(x)+4*C*numpy.sin(2*x)+4*D*numpy.cos(2*x))

      # Boyd's method in 4 dimensions: Boyd, J.P.; J. Eng. Math. (2006) 56:203-219
      d = D - 1.0j * C                        
      K = numpy.zeros((4,4), numpy.complex64)
      K[0,1] = 1.0
      K[1,2] = 1.0
      K[2,3] = 1.0
      K[3,0] = -(D + 1.0j * C)/d
      K[3,1] = -(B + 1.0j * A)/d
      K[3,2] = 0.0
      K[3,3] = -(B - 1.0j * A)/d
      E, X = numpy.linalg.eig(K)
      X   = -1.0j * numpy.log(E)
      X[numpy.where(X<0.0)] += 2.0 * numpy.pi
      X = X.real

      # discriminate between minima and maxima
      Xmin = list()
      Xmax = list()
      for x in X:
          F = self._f (x,A,B,C,D)
          g = self._fg(x,A,B,C,D)
          if   g> 0.0: Xmin.append(x)
          elif g< 0.0: Xmax.append(x)
          #else: raise ValueError, "The Hessian of objective function is zero at X=%15.5f" % x
      Xmin = numpy.array(Xmin)
      Xmax = numpy.array(Xmax)
    
      # Find optimal gamma 
      gamma = None
      if opt.startswith('min'):
         Zold = 1.0e8
         for x in Xmin:
             Z = self._eval_Z(self._form_X(i, j, x), self._R, self._P)
             if Z < Zold:
                gamma = x
             Zold = Z
      else:
         Zold = -1e8
         for x in Xmax:
             Z = self._eval_Z(self._form_X(i, j, x), self._R, self._P)
             if Z > Zold:
                gamma = x
             Zold = Z
      assert gamma is not None, "Error while searching for optimum!"
      return gamma

  def _get_ABCD(self, i, j):
      """Retrieve ABCD parameters for root search"""
      A = self._P[i]+self._P[j]
      B = self._P[i]-self._P[j] 
      C =-2.*(self._R[i,j]+self._R[j,i])
      D = 2.*(self._R[j,j]-self._R[i,i])
      ii = numpy.setdiff1d(numpy.arange(self._P.size), [i,j])
      A -= self._R[j,ii].sum() + self._R[i,ii].sum() + self._R[ii,j].sum() + self._R[ii,i].sum()
      B += self._R[j,ii].sum() - self._R[i,ii].sum() + self._R[ii,j].sum() - self._R[ii,i].sum()
      return A, B, C, D

  def _eval_Z(self, X, R, P):
      """Evaluate the objective Z function"""            
      z1 = numpy.dot(X, numpy.dot(R,X.T))
      z2 = numpy.dot(X, P)
      return z1.sum() - z2.sum()

  def _eval_dZ(self, g, P, R, i, j):
      """Compute the change in Z"""
      dZ = (1.0 - numpy.cos(g)) * (P[i] + P[j]) \
                + numpy.sin(g)  * (P[i] - P[j])       \
                + numpy.sin(2.*g) * (R[j,j] - R[i,i]) \
             - 2.*numpy.sin(g)**2 * (R[i,j] + R[j,i])
      ii = numpy.setdiff1d(numpy.arange(P.size), [i,j])
      dZ-= (1-numpy.cos(g))*(R[j,ii].sum() + R[i,ii].sum() + R[ii,j].sum() + R[ii,i].sum())
      dZ+=    numpy.sin(g) *(R[j,ii].sum() - R[i,ii].sum() + R[ii,j].sum() - R[ii,i].sum())
      return dZ

  def _form_X(self, i, j, gamma):
      """Form unitary matrix X"""                                        
      X = numpy.identity(self._d)
      T = numpy.cos(gamma)
      g = numpy.sin(gamma)
      X[i,i] = T
      X[j,j] = T
      X[i,j] = g
      X[j,i] =-g
      return X
 
  # private
  def __same(self, x, X, tol=0.0001):
      """Check if the solution was already obtained"""
      #pipi = 2*numpy.pi
      result = False
      for xp in X:
          dx = x-xp
          if (numpy.abs(dx) < tol) or (numpy.abs(dx+2*numpy.pi) < tol) or (numpy.abs(dx-2*numpy.pi) < tol): 
          #if not (numpy.abs(dx)%pipi):
             result = True
             break
      return result

def glue(a, b, axis=None, buff=10, update_time_axis=None, lprint=False):
    """
 --------------------------------------------------------------------------------------------
 Concatenate overlapping numpy.ndarrays. 
 --------------------------------------------------------------------------------------------
 Assumes that 'a' ends with slice that overlaps with the beginning of 'b'. 'a' and 'b' can 
 have arbitrary but equal shape. If overlap is found returns the concatenated numpy.ndarray,
 otherwise returns the None object.

 Usage:

    c = glue(a, b, axis=None, buff=10, update_time_axis=None, lprint=False)
 
 where:

    o a, b             - input numpy.ndarrays containing numbers
    o c                - the resulting concatenated numpy.ndarray                      
    o axis             - column id to compare the overlap. If index=None it is assumed
                         'a' and 'b' have only one axis (row vectors).
    o buff             - the buffer size. Overlap length should be less than buffer
    o update_time_axis - column which is to be updated based on the constant increment
                         assumed in 'a[:,update_time_axis]'. If update_time_axis=None 
                         no update is done.
    o lprint           - print the square differences of trial slices or not

 --------------------------------------------------------------------------------------------
                                                                 Last Revision: 31 Aug 2017
"""
    assert len(a.shape)==len(b.shape), " Error: The shapes of input arrays are not the same!"
    assert buff <= min(len(a),len(b)), " Error: Buffer size too large!"
    # parse columns
    if len(a.shape)>1: 
       if axis is None: 
          raise ValueError, " Error: Column index has not ben specified. Specify the index of column to analyze by setting index=<id> where <id> is in Python convention!"
       a2 = a[:,axis]; b2 = b[:,axis]
    else: 
       a2 = a.copy(); b2 = b.copy()
    # search for overlap
    i_ov = None
    for i in range(1,buff+1):
        a2_slice = a2[-i:  ]
        b2_slice = b2[  : i]
        dsq = sum( (a2_slice - b2_slice)**2 )
        if lprint: print " Trial overlap length= %4d Square difference= %14.6f" % (i,dsq)
        if dsq==0.0000000: i_ov = i
    # check if overlap found
    if i_ov is None: 
        print " Warning: The sequences are not mergeable or you must increase the buffer > %d" % buff
        return None
    b_copy = b.copy()
    # update time axis
    if len(a.shape)>1 and update_time_axis is not None:
        assert update_time_axis!=index, " Error: Glue column cannot be the time column!"
        dt_a = a[:,update_time_axis][1] - a[:,update_time_axis][0]
        t_a_last = a[:,update_time_axis][-1]
        delta_t_a = dt_a * (len(a)+0)
        b_copy[:,update_time_axis] += delta_t_a
    # concatenate
    c = numpy.concatenate([a,b_copy[i_ov:]]) 
    return c


class TIMER:
    """process timing statistics utility class""" 
    
    def __init__(self, name=None):
        self.name = name
        self.occurence = 'start'
        self.t0 = time.time()
        self.tp = self.t0
        self.occurence_list = []
        self.log = "\n"
        self.log+= " ------------------ \n"
        self.log+= " Timing information \n"
        self.log+= " ------------------ \n"
        self.log+= "\n" 
        if name: self.log+= " %s\n\n" % name#.capitalize()
        self.total_time = 0
        
    def actualize(self,new_occurence):
        """actualizes new occurence in the clock history"""
        
        self.occurence = new_occurence
        self.measure()
        self.total_time = sum(self.occurence_list)
        return
        
    def measure(self):
        """measures length of occurence"""
        
        self.tn = time.time()
        self.length = self.tn - self.tp
        self.occurence_list.append(self.length)
        self.log += " - %60s %14.5f sec\n" % (self.occurence.ljust(60),self.length)
        self.tp = self.tn
        return
        
    def __repr__(self):
        """Print the timing of entire process"""
        
        suma = sum(self.occurence_list)
        t = self.total_time
        log = self.log.split('\n')
        N = 6 if self.name is None else 8
        for i in range(len(log) - N):
            log[i+N-1] += " (%4.1f%%)" % (self.occurence_list[i]/suma*100.0)
        self.log = '\n'.join(log)
        self.log+= '\n'
        self.log+= " =========================================================================================\n"
        self.log+= " TOTAL TIME:  %d days %d hours %d min %d sec \n" % ( t/86400,t/3600,t/60,int(t) )
        self.log+= "\n\n"
        return str(self.log)


def distribute(box_size, mol_size, max_n_mol, box_type='cubic'):
    """
Distribute molecules into a regular lattice in a cubic box. 
Finds optimial number of molecules per each box dimension to fill the 
space. Pay caution: to many molecules per box may cause clashes and the algorithm
will produce an error message.

Usage:

   nxyz = distribute(box_size, mol_size, max_n_mol, box_type='cubic')
 
Arguments:

   box_size = [Lx, Ly, Lz] - list or numpy.ndarray of length 3
                             specifying the box dimensions in unit U
   mol_size = [hy, hy, hz] - list or numpy.ndarray of length 3
                             specifying the molecule dimensions in unit U    
   max_n_mol               - maximal number of molecules to be packed
   box_type                - type of a box. Now only cubic box is implemented

Returns:

   nxyz = [nx, ny, nz]     - numpy.ndarray of length 3 giving the optimal
                             fractional number of molecules per axis
"""
    assert box_type=='cubic', 'Box type %s is not implemented yet! Quitting...' % box_type
    n = max_n_mol**(1./3)
    nxyz_0 = numpy.array([n,n,n], numpy.float64)
    Lx, Ly, Lz = box_size
    hx, hy, hz = mol_size
    fmin_slsqp = scipy.optimize.fmin_slsqp

    # minimizing function                                                                                            
    I = lambda nxyz, n, Lx, Ly, Lz, hx, hy, hz: \
        (Lx-nxyz[0]*hx)/(nxyz[0]-1) + \
        (Ly-n/(nxyz[0]*nxyz[2])*hy)/(n/(nxyz[0]*nxyz[2])-1) + \
        (Lz-n/(nxyz[0]*nxyz[1])*hz)/(n/(nxyz[0]*nxyz[1])-1)

    # constraints
    dx = lambda nxyz, n, Lx, Ly, Lz, hx, hy, hz: (Lx-nxyz[0]*hx)/(nxyz[0]-1)
    dy = lambda nxyz, n, Lx, Ly, Lz, hx, hy, hz: (Ly-n/(nxyz[0]*nxyz[2])*hy)/(n/(nxyz[0]*nxyz[2])-1)
    dz = lambda nxyz, n, Lx, Ly, Lz, hx, hy, hz: (Lz-n/(nxyz[0]*nxyz[1])*hz)/(n/(nxyz[0]*nxyz[1])-1)

    # go! 
    x = fmin_slsqp(I, nxyz_0, eqcons=[dx, dy, dz], f_eqcons=None, ieqcons=(), f_ieqcons=None,
                      bounds=(), fprime=None, fprime_eqcons=None, fprime_ieqcons=None,
                      args=(max_n_mol, Lx, Ly, Lz, hx, hy, hz), iter=100, acc=1e-06, iprint=1, disp=None, full_output=0,
                      epsilon=1.4901161193847656e-08)
    return x

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def make_bqc_inp(mol):
    nelec = mol.get_atno().sum()
    bfs   = mol.get_bfs()

    ntypes = {'[0 0 0]':1 ,
              '[1 0 0]':2 , '[0 1 0]':3 , '[0 0 1]':4 ,
              '[2 0 0]':5 , '[0 2 0]':6 , '[0 0 2]':7 , '[1 1 0]':8 , '[1 0 1]':9 , '[0 1 1]':10,
              '[3 0 0]':11, '[0 3 0]':12, '[0 0 3]':13, '[2 1 0]':14, '[2 0 1]':15, '[1 2 0]':16, 
                            '[0 2 1]':17, '[1 0 2]':18, '[0 1 2]':19, '[1 1 1]':20}
    
    # gather the operational arguments
    natoms = len(mol)
    nbasis = len(bfs)
    nt     = bfs.get_bfst()
    
    # collect the basis set
    exps  = list()
    coefs = list()
    origs = list()
    nfirst= list()
    nlast = list()
    
    ngmx = 0
    for i in range(nbasis):
        bf = bfs.bfs[i]
        n_cont = len(bf.prims)
        nfirst.append(ngmx+1)
        nlast.append(ngmx+n_cont)
        for j in range(n_cont):
            origs += list(bf.prims[j].origin) 
        orig  = bf.origin
        exps +=bf.pexps
        coefs+=bf.pcoefs
        #
        ngmx += n_cont
    
    # build the final data structures
    vlist = numpy.zeros((natoms,4),numpy.float64)
    eta   = numpy.zeros((ngmx  ,5),numpy.float64)
    
    # make vlist
    vlist[:,:3] = mol.get_pos()
    vlist[:,3] = mol.get_atno()
    # make eta
    eta[:,:3] = numpy.array(origs,numpy.float64).reshape(ngmx,3)
    eta[:,3] = numpy.array(exps ,numpy.float64)
    eta[:,4] = numpy.array(coefs,numpy.float64)
    # make other
    ncntr = numpy.array(bfs.LIST1, int) + 1
    ncmx  = natoms
    nbfns = nbasis
    # make ntype
    ntype = [ntypes[str(i)] for i in nt]
    task = 90
    crit = 1.00E-07
    damp = 0.33
    interp = 30
    read_eri = 0 # NO. For YES set 1
    charge = 0
    multiplicity = 1
    ntmx = 20
    return vlist, eta, ncntr, ntype, nfirst, nlast, ncmx, ngmx, ntmx, nbfns

if import_matplotlib:
   def circles(x, y, s, c='b', ax=None, vmin=None, vmax=None, **kwargs):
       """
       Make a scatter of circles plot of x vs y, where x and y are sequence 
       like objects of the same lengths. The size of circles are in data scale.
   
       Parameters
       ----------
       x,y : scalar or array_like, shape (n, )
           Input data
       s : scalar or array_like, shape (n, ) 
           Radius of circle in data scale (ie. in data unit)
       c : color or sequence of color, optional, default : 'b'
           `c` can be a single color format string, or a sequence of color
           specifications of length `N`, or a sequence of `N` numbers to be
           mapped to colors using the `cmap` and `norm` specified via kwargs.
           Note that `c` should not be a single numeric RGB or
           RGBA sequence because that is indistinguishable from an array of
           values to be colormapped.  `c` can be a 2-D array in which the
           rows are RGB or RGBA, however.
       ax : Axes object, optional, default: None
           Parent axes of the plot. It uses gca() if not specified.
       vmin, vmax : scalar, optional, default: None
           `vmin` and `vmax` are used in conjunction with `norm` to normalize
           luminance data.  If either are `None`, the min and max of the
           color array is used.  (Note if you pass a `norm` instance, your
           settings for `vmin` and `vmax` will be ignored.)
   
       Returns
       -------
       paths : `~matplotlib.collections.PathCollection`
   
       Other parameters
       ----------------
       kwargs : `~matplotlib.collections.Collection` properties
           eg. alpha, edgecolors, facecolors, linewidths, linestyles, norm, cmap
   
       Examples
       --------
       a = np.arange(11)
       circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')
   
       License
       --------
       This code is under [The BSD 3-Clause License]
       (http://opensource.org/licenses/BSD-3-Clause)
       """
       from matplotlib.patches import Circle
       from matplotlib.collections import PatchCollection
       import pylab as plt
       #import matplotlib.colors as colors
   
       if ax is None:
           ax = plt.gca()    
   
       if isinstance(c,basestring):
           color = c     # ie. use colors.colorConverter.to_rgba_array(c)
       else:
           color = None  # use cmap, norm after collection is created
       kwargs.update(color=color)
   
       if isinstance(x, (int, long, float)):
           patches = [Circle((x, y), s),]
       elif isinstance(s, (int, long, float)):
           patches = [Circle((x_,y_), s) for x_,y_ in zip(x,y)]
       else:
           patches = [Circle((x_,y_), s_) for x_,y_,s_ in zip(x,y,s)]
       collection = PatchCollection(patches, **kwargs)
   
       if color is None:
           collection.set_array(np.asarray(c))
           if vmin is not None or vmax is not None:
               collection.set_clim(vmin, vmax)
   
       ax.add_collection(collection)
       return collection
   
   class PotentialContourMap:
      """
    ---------------------------------------------------------------------------------------------------------------------------------------
    Represents the Potential Contour Map for a molecule. 
   
                                                                                                                   Author: Bartosz Błasiak
    ---------------------------------------------------------------------------------------------------------------------------------------
   
    DESCRIPTION:
   
       Makes the contour plot of generalized potential (normal elextrostatic, solvatochromic transition etc) for a given molecular data.
    It masks areas where multipole expansion diverges.  
   
   
    USAGE:
   
               from libbbg.utilities import PotentialContourMap as PCM                  
                                                                                        
               map = PCM(dma, atoms, bonds, allign_atid, allign_axes=(1,2,0), 
                         pad_x=4.0, pad_y=4.0, dx=0.5, dy=0.5, levels=None,
                         pad_left=None, pad_right=None, pad_down=None, pad_up=None, 
                         radii=None, dmat=None, vec=None, basis=None, 
                         dma_pot=None, qm_mask_thr=0.01, 
                         colors =[(0.0, 0.1, 1.0),
                                  (1.0, 1.0, 1.0),
                                  (1.0, 0.1, 0.0)], levs=60, linthresh = 0.0014,
                         label=False, fmt='%2.1f', font_size=6, block=False, name=None,
                         atom_colors=None, bond_width=4.0, bond_fmt='k-',
                         )
                                                                                        
               map.make()
   
   
    ARGUMENTS:
   
    Basic inputs:
    -------------
    dma         - libbbg.dma.DMA object
    
    atoms       - list of atoms to be drawn on the plane of a contour plot. (normal numbers)
                  Example: [1,2,3,4]
    
    bonds       - list of lists of bonds to be drawn on the plane of a contour plot. (normal numbers)
                  Example: [[1,2], [3,2]] 
   
    allign_atid - allignment specification for the molecule. Describes the new axes of coordinate system. 
                  See libbbg.utilities.Allign class for more information.
   
    allign_axes - permutation indices for axes. 
                  See libbbg.utilities.Allign class for more information.
                  Example: allign_atid=[3,1,4] 
                           allign_axes=[1,2,0]
                           it places 3-rd atom at the origin, 1-st atom specifies y-axis, 
                           4-th atom lies in the plane of contour plot, z-axis is perpendicular to the plane of contour plot.
   
    max_moment  - maximum rank of multipole moment included in map generation (for dma object only). Default is 3 (octupole),
                  2 - quadrupole, 1 - dipole and 0 - charge).
   
   
    Masking options:
    ----------------
   
    Choose either *radii* or the remaining options. 
   
    radii       - van der Waals radii for atoms. It is used to mask the regions where multipole expansion diverges. 
                  If radii=None you must specify density matrix, basis set and the appropriate DMA object computed at the same level of theory. 
                  Then the region with deviations 
                  larger than qm_mask_thr (in absolute value in AU units of potential) from the exact QM potential are masked.
   
                                - OR -
   
    basis       - basis set (eg. 6-311++G**)
   
    dmat        - density matrix for a molecule in alligned geometry! To obtain this geometry first make a map using van der Waals radii. The
                  program will print the alligned geometry. For THIS geometry compute the density matrix. Otherwise the result will be incorrect.
   
    vec         - optional to dmat, specify the eigenvectors (wavefunction LCAO coefficients). These eigenvectors will be rotated so
                  they should be in the geometry identical to the DMA supplied (dma). From the rotated eigenvectors density matrix is computed.
                  This option is still buggy and provides poor masking since rotation of wavefunction is either buggy or just not accurate enough.
   
    dma_pot     - dma distribution in exactly the same orientation as DMA suppplied (the one which is to be plotted). dma_pot represents
                  the elecrostatic potential and is assumed to have the same divergence spheres as DMA to be plotted.
   
    atnos       - atomic numbers of all atoms in DMA distribution. 
   
   
    Plot options:
    -------------
   
    pad_x/pad_y - padding in x/y directions (in Bohr)
   
    dx/dy       - grid spacing in x/y directions (in Bohr)
   
    levels      - Levels of contour plot. If None, default ones will be used by matplotlib.contour
   
    label       - if True, the isovalues of potential are printed within isobars. Default is False.
   
    fmt         - format of labels. Relevant if label=True.
   
    font_size   - Size of labels. Relevant if label=True.
   
    block       - show the plot before ending the execution of the main script. Default is False.
   
    name        - if not None, the map will be saved to the file=name. Default is None.
   
    atom_colors - list of size len(atoms) specifying matplotlib.colors for atoms. Default all black.
   
    bond_width  - number specifying width of stick representation of bonds. Default 4.0
   
    bond_fmt    - format of bond line. Default 'k-'.
   
    colors      - RGB values for color map. Default is Red(+)-White(0)-Blue(-)
   
    levs        - number of levels in the colorbar
   
    linthresh   - Number specifying the maximum value (for positive and negative potential regions) where the normalization is linear. 
                  In other points the normalization is logarithmic. This is necessary parameter because otherwise the zero-valued regions
                  would blow up to infinity during logarithmic normalization.
   
   
    ---------------------------------------------------------------------------------------------------------------------------------------
                                                                                                         Last Revision: 14 Jan 2015
   """
      def __init__(self, dma, atoms, bonds, allign_atid, allign_axes=(1,2,0), max_moment=3,
                         pad_x=4.0, pad_y=4.0, dx=0.5, dy=0.5, levels=None,
                         pad_left=None, pad_right=None, pad_down=None, pad_up=None, 
                         radii=None, vec=None, dmat=None, dma_pot=None, qm_mask_thr=0.10,
                         atnos=None, basis=None,
                         colors =[(0.0, 0.1, 1.0),
                                  (1.0, 1.0, 1.0),
                                  (1.0, 0.1, 0.0)], levs=60, linthresh = 0.0014,
                         label=False, fmt='%2.1f', font_size=6, block=False, name=None, 
                         atom_colors=None, bond_width=4.0, bond_fmt='k-',
                         ):
          self.__dma = dma.copy()
          self.__atoms = numpy.array(atoms,int)-1
          self.__bonds = numpy.array(bonds,int)-1
          self.__levels = levels
          self.__allign = (allign_atid, allign_axes)
          self.__max_moment = max_moment
          if radii is None: self.__radii = None
          else: self.__radii = numpy.array(radii)
          self.__atnos = atnos
          self.__basis = basis
          self.__vec = vec
          self.__dmat = dmat
          self.__dma_pot = dma_pot; self.__qm_mask_thr = qm_mask_thr
          if pad_left is None and pad_right is not None or pad_right is None and pad_left is not None:
             print " You must specify padding in both right and left directions!"; exit()
          if pad_down is None and pad_up is not None or pad_up is None and pad_down is not None:
             print " You must specify padding in both up and down directions!"; exit()
          if pad_left is None and pad_right is None:
             self.__pad_left = pad_x ; self.__pad_right= pad_x
          else:
             self.__pad_left= pad_left ; self.__pad_right= pad_right
          if pad_down is None and pad_up is None:
             self.__pad_down= pad_y ; self.__pad_up= pad_y
          else:
             self.__pad_down= pad_down ; self.__pad_up= pad_up
          self.__dx = dx       ; self.__dy = dy
          self.__colors = colors; self.__levs = levs; self.__linthresh = linthresh
          self.__label = label ; self.__fmt = fmt; self.__font_size = font_size
          self.__block = block ; self.__name = name
          self.__atom_colors = atom_colors
          self.__bond_width = bond_width; self.__bond_fmt = bond_fmt
          if atom_colors is None: 
             self.__atom_colors = ['black'] * len(atoms)
          if radii is None:
             error = " Van der Waals radii not specified so you must provide basis set (basis) and density matrix (dmat) for alligned molecule!"
             assert (dmat is not None or vec is not None) , error
             self.__mask_with_qm = True
          else: self.__mask_with_qm = False
   
          self._create_cmap()
          self._prepare()
   
      def make(self):
          """Main routine for map generation"""
          self._allign_molecule()
          self._calc_potential()
          self._mask()
          self._make_plot()
          return
   
      def _allign_molecule(self):
          """Alligns the molecule (structure and DMA)"""
          atid, axes = self.__allign
          alligner = Allign(self.__xyz, atid=atid, axes=axes, dma=self.__dma)
          if self.__mask_with_qm:
             alligner2 = Allign(self.__xyz, atid=atid, axes=axes, dma=self.__dma_pot)
          self.__dma, self.__xyz = alligner.get_transformed()
          self.__dma.MAKE_FULL()
          self.__dma.MakeTraceless()
          self.__dma.makeDMAfromFULL()
          rot = alligner.rot
   
          # rotate (allign) the additional electrostatic distribution
          if self.__mask_with_qm:
             self.__dma_pot, smiec = alligner2.get_transformed()
             self.__dma_pot.MAKE_FULL()
             self.__dma_pot.MakeTraceless()
             self.__dma_pot.makeDMAfromFULL()
             self._create_mol()
          print " The alligned coordinates [Angstrom]:\n"
          PRINTL(self.__xyz*units.UNITS.BohrToAngstrom,'','')
          # generate the rotated density matrix
          if self.__mask_with_qm:
             if self.__dmat is None:
                typs= self.__bfs.get_bfst().sum(axis=1)                       
                nbasis = len(self.__bfs)
                self.__vec = qm.efprot.vecrot(self.__vec, rot.T, typs)
                self.__dmat = numpy.zeros((nbasis,nbasis),numpy.float64)
                for i in range(len(self.__vec)):
                    v = self.__vec[i]
                    self.__dmat += numpy.outer(v,v)
                self.__dmat = self.__dmat.ravel() * 2.0   
             else:
                self.__dmat = self.__dmat.ravel()
          return
   
      def _calc_potential(self):
          """Calculate potential from DMA and WFN if needed"""
          ndma = numpy.array([len(self.__dma),],int)
          rdma = self.__dma.get_origin().ravel()
          chg  = self.__dma.get_charges()
          dip  = self.__dma.get_dipoles().ravel()
          qad  = self.__dma.get_quadrupoles().ravel()
          oct  = self.__dma.get_octupoles().ravel()
          points = numpy.zeros((self.__nx*self.__ny*8), numpy.float64)
          points = qm.make_points.make_points(self.__nx,self.__ny,self.__dx,self.__dy,self.__x_min,self.__y_min,points)
          if self.__mask_with_qm: 
             points2 = points.copy()
             points3 = points.copy()
          points = qm.clemtp.potdma(points,rdma,ndma,chg,dip,qad,oct)
          points = points.reshape(self.__nx,self.__ny,8)
          if self.__max_moment == 3:
             self.__z = points[:,:,-1]
          elif self.__max_moment == 2:
             self.__z = points[:,:,-5] + points[:,:,-4] + points[:,:,-3]
          elif self.__max_moment == 1:
             self.__z = points[:,:,-5] + points[:,:,-4]
          elif self.__max_moment == 0:
             self.__z = points[:,:,-5]
          if self.__mask_with_qm:
             # exact QM potential
             vlist, eta, ncntr, ntype, nfirst, nlast, ncmx, ngmx, ntmx, nbfns = make_bqc_inp(self.__mol)
             points2 = qm.wfn.scawf2(points2,self.__dmat,eta,nfirst,nlast,ntype,vlist)
             points2 = points2.reshape(self.__nx,self.__ny,8)
             self.__z_exact = points2[:,:,4]
             # potential from DMA
             ndma = numpy.array([len(self.__dma_pot),],int)  
             rdma = self.__dma_pot.get_origin().ravel()
             chg  = self.__dma_pot.get_charges()
             dip  = self.__dma_pot.get_dipoles().ravel()
             qad  = self.__dma_pot.get_quadrupoles().ravel()
             oct  = self.__dma_pot.get_octupoles().ravel()
             points3 = qm.clemtp.potdma(points3,rdma,ndma,chg,dip,qad,oct)
             points3 = points3.reshape(self.__nx,self.__ny,8)
             self.__z_test = points3[:,:,-1]
             del points2, points3
          del points
          return
   
      def _mask(self):
          """masks the region where multipole expansion diverges"""
          # masking according to QM calculation
          if self.__mask_with_qm:
             ke = self.__z_exact.reshape(self.__nx*self.__ny)
             kt = self.__z_test .reshape(self.__nx*self.__ny)
             #for i in range(len(ke)):
             #    ff = numpy.abs(ke[i]-kt[i])/numpy.abs(ke[i]) * 100.0
             #    print " %16.5E %16.5E %16.5f %16.5f" % (ke[i], kt[i], numpy.abs(ke[i]-kt[i]), ff)
             #interior = numpy.abs(self.__z_exact - self.__z_test)/numpy.abs(self.__z_exact) > self.__qm_mask_thr
             interior = numpy.abs(self.__z_exact - self.__z_test) > self.__qm_mask_thr
             self.__z[interior] = numpy.NaN #numpy.ma.masked
          # masking according to vdW radii
          else:
              X,Y = numpy.meshgrid(self.__y,self.__x)
              atoms = self.__xyz[self.__atoms]
              for i in range(len(atoms)):
                  interior = numpy.sqrt(((X-atoms[i,1])**2) + ((Y-atoms[i,0])**2)) < self.__radii[i]
                  self.__z[interior] = numpy.ma.masked
          return
   
      def _make_plot(self):
          """Makes the contour plot"""
          # [1] mask
          X,Y = numpy.meshgrid(self.__y,self.__x)
          Z = self.__z 
   
          # [2] make contour plot
          CP1 = pylab.contour(X,Y,Z,levels=self.__levels,colors='k')
          if self.__label: pylab.clabel(CP1, colors='k', fmt=self.__fmt, nline=True, fontsize=self.__font_size)
          pylab.gca().patch.set_color('1.0')
          CP2 = pylab.contourf(X,Y,Z,levels=self.__levels,cmap=self.__cmap,norm=matplotlib.colors.SymLogNorm(self.__linthresh))
          pylab.colorbar(CP2)
          pylab.title('Contour plot')
          pylab.xlabel('x')
          pylab.ylabel('y')
   
          x = self.__xyz[:,0]
          y = self.__xyz[:,1]
   
          # [3] plot circles on mapped areas
          if not self.__mask_with_qm:
             s = self.__radii[self.__atoms]
   
             circles(y[self.__atoms], x[self.__atoms], s=s, alpha=1.0, c='white',facecolor='white')
   
          # [4] plot bonds
          for bond in self.__bonds:
              pylab.plot(y[bond], x[bond], self.__bond_fmt, lw=self.__bond_width)
   
          # [5] plot atoms
          for i in range(len(self.__atoms)):
              circles(y[self.__atoms[i]], x[self.__atoms[i]], s=0.3, alpha=1.0, c='black',facecolor=self.__atom_colors[i])
             
          
          pylab.axes().set_aspect('equal', 'datalim')
          if self.__name is not None: pylab.savefig(self.__name)
          pylab.show(block=self.__block)
          return
   
      # helper methods
      def _create_mol(self):
          """Make PyQuante.Molecule object in the alligned geometry"""
          xyz = self.__xyz.copy()
          Coords = list()
          for i in range(len(xyz)):
              atom  = (self.__atnos[i], (xyz[i,0],
                                         xyz[i,1],
                                         xyz[i,2]) )
              Coords.append(atom)
          Mol = PyQuante.Molecule('buddy',Coords,units='Bohr',
                                  multiplicity=1,charge=0,
                                  basis=self.__basis)
          self.__mol = Mol
          self.__bfs = Mol.get_bfs()
          return
    
      def _prepare(self):
          """Generate the X and Y axis and initialize points' values Z"""
          xyz = self.__dma.get_pos()
          x_min = xyz[:,0].min()-self.__pad_left; x_max = xyz[:,0].max()+self.__pad_right
          y_min = xyz[:,1].min()-self.__pad_down; y_max = xyz[:,1].max()+self.__pad_up
          nx = numpy.int64((x_max-x_min)/self.__dx + 1)
          ny = numpy.int64((y_max-y_min)/self.__dy + 1)
          x = numpy.linspace(x_min, x_max, nx)
          y = numpy.linspace(y_min, y_max, ny)
          self.__xyz = xyz
          self.__x_min = x_min; self.__x_max = x_max
          self.__y_min = y_min; self.__y_max = y_max
          self.__nx = nx
          self.__ny = ny
          self.__x = x
          self.__y = y
          self.__z = numpy.zeros((nx,ny),numpy.float64)
          print " This map will contain %10d points" % (nx*ny)
          return 
   
      def _create_cmap(self):
          levs = range(self.__levs)
          assert len(levs) % 2 == 0, 'N levels must be even.'
          cmap = matplotlib.colors.LinearSegmentedColormap.from_list(name='red_white_blue', 
                                                                     colors =self.__colors, 
                                                                     N=len(levs)-1, )
          self.__cmap = cmap
          return

def autocorr(y, normalize=True, t_min=0.0, t_max=1.0 ):
    """Auto-correlation"""
    result = numpy.correlate(y, y, mode='full')
    r = result[result.size/2:]
    if normalize:
       M = lambda L, t, dt: L - t/dt
       L = len(y); dt = (t_max - t_min)/(L-1)
       for i in range(L):
           t = t_min + i*dt
           r[i]/= M(L, t, dt)
       r/=r[0]
    return r

def crosscorr(y1, y2, normalize=True, t_min=0.0, t_max=1.0 ):
    """Cross-correlation"""
    assert y1.shape==y2.shape, 'Shapes of input arrays are not the same!'
    result = numpy.correlate(y1, y2, mode='full')
    r = result[result.size/2:]
    if normalize:
       M = lambda L, t, dt: L - t/dt
       L = len(y); dt = (t_max - t_min)/(L-1)
       for i in range(L):
           t = t_min + i*dt
           r[i]/= M(L, t, dt)
       r/=r[0]
    return r

class FF_scheme(object):
   """
 ------------------------------------------------------------------------------
 Finite differentiation scheme

 Contains various displacement schemes for numerical differentiation
 (generates the displacement coordinates)
 ------------------------------------------------------------------------------
                                                 Last revision: 5 Nov 2014
"""
   def __init__(self, step, scheme='3pt', DIM=3):
       self._step   = step
       self._scheme = scheme.lower()
       self._DIM    = DIM
       self._disp1 , self._disp2 = self._make()

   # PUBLIC

   def get(self):
       return self._disp1, self._disp2

   # PROTECTED

   def _get_Ki(self, I):
       """Get the principal index of first derivative displacement matrix"""
       if self._scheme=='3pt':
          return I*2

   def _get_Kij(self, I, J):
       """Get the principal index of second derivative displacement matrix"""
       if self._scheme=='3pt':
          S = I*(2*self._DIM-I-1)/2 + J - I - 1
          return S*2

   def _make(self):
       """Generate the displacement vectors"""
       k   = self._step
       DIM = self._DIM
       if self._scheme == '3pt':
          disp1 = numpy.zeros((DIM* 2     , DIM), numpy.float64)
          disp2 = numpy.zeros((DIM*(DIM-1), DIM), numpy.float64)
          for i in range(DIM):
              K = 2*i
              disp1[K  ,i  ] =+k
              disp1[K+1,i  ] =-k
          I = 0
          for i in range(DIM):
              for j in range(i+1,DIM):
                  K = 2*I
                  disp2[K  ,i] =+k
                  disp2[K  ,j] =+k
                  disp2[K+1,i] =-k
                  disp2[K+1,j] =-k
                  I+=1
          return disp1, disp2
       else:
          raise Exception('The scheme %s is not implemented!' % self._scheme)


class diff(FF_scheme):
   """
 ------------------------------------------------------------------------------
 Numerical Finite Difference Method

 Calculates 1-st and 2-nd derivatives of a function f(x) at point x_0.
 ------------------------------------------------------------------------------

 Usage:
   c = diff(func=None, step=0.001, scheme='3pt', DIM=3)
   fder, sder = c.eval(x_0, symm=True)

 where:
   o func       function defining f(x). It should take x as a first argument
                If function is None then you must provide function evaluations
                for every displaced structure in eval (according to scheme pointity)
   o step       differentiation step
   o scheme     numerical scheme: 3pt - 3-point central method
   o DIM        dimension of x_0 (number of variables)

 Return:
   o fder       numpy.ndarray of shape DIM
   o sder       numpy.ndarray of shape DIM x DIM
 ------------------------------------------------------------------------------
                                                 Last revision: 5 Nov 2014
"""
   def __init__(self, func=None, step=0.001, scheme='3pt', DIM=3):
       super(diff, self).__init__(step=step, scheme=scheme.lower(), DIM=DIM)
       self.__func = func

   def eval(self, x, symm=True):
       """
If func is not None:
  Evaluate the derivatives of f(x). Usage: eval(x, symm=True) where x are arguments of f(x)
  and if symm=True then do the symmetrization of the offdiagonal derivatives (default).

If func is None:
  then x means the set of function evaluations for an appropriate displaced
  arguments (they must be consistet with this package pointity scheme conventions!).
"""
       if self.__func is not None: fder, sder = self._der(x, symm)
       else:                       fder, sder = self._der_ext(x[0],x[1], symm)
       return fder, sder

   def _make_disp(self, x, disp):
       """create the displaced arguments"""
       n = len(disp)
       x_disp = numpy.zeros( (n, self._DIM), numpy.float64 )
       for i in range(n):
           x_disp[i] = x.copy() + disp[i]
       return x_disp

   def _der(self, x, symm):
       """calculate the derivatives - now only 3pt scheme is implemented"""
       # make displaced coordinates
       x_d1 = self._make_disp(x, self._disp1)
       x_d2 = self._make_disp(x, self._disp2)
       # initialize the derivatives
       fder = numpy.zeros( self._DIM,             numpy.float64)
       sder = numpy.zeros((self._DIM, self._DIM), numpy.float64)
       # calculate first and diagonal second derivatives
       for i in range(self._DIM):
           Ki  = self._get_Ki(i)
           f1 = self.__func( x_d1[Ki+0] ) - self.__func( x_d1[Ki+1] )
           f1/= 2.*self._step
           fder[i] = f1
           f2 = self.__func( x_d1[Ki+0] ) + self.__func( x_d1[Ki+1] ) - 2.*self.__func( x )
           f2/= self._step**2.0
           sder[i,i] = f2
           # calculate offdiagonal second derivatives
           for j in range(self._DIM):
             if i!=j:
               Kj = self._get_Ki(j)
               if i>j: 
                  I = j
                  J = i
               else:
                  I = i
                  J = j
               Kij = self._get_Kij(I,J)
               f2 = self.__func( x_d2[Kij ]      ) + self.__func( x_d2[Kij+1]     ) \
                  - self.__func( x_d1[Ki+0]      ) - self.__func( x_d1[Kj+0]      ) \
                  - self.__func( x_d1[Ki+1]      ) - self.__func( x_d1[Kj+1]      ) \
                  + 2.*self.__func( x )
               f2/= 2.*self._step**2.
               sder[i,j] = f2
       # symmetrize the second derivatives
       if symm:
          tl = numpy.tril(sder,k=-1)
          tu = numpy.triu(sder,k=+1).transpose()
          av =(tl + tu)/2.0
          sder = numpy.diag(sder.diagonal()) + av + (av.copy()).transpose()
       return fder, sder

   def _der_ext(self, f_d1, f_d2, symm):
       """calculate the derivatives - now only 3pt scheme is implemented"""
       # initialize the derivatives
       fder = numpy.zeros( self._DIM,             numpy.float64)
       sder = numpy.zeros((self._DIM, self._DIM), numpy.float64)
       # calculate first and diagonal second derivatives
       for i in range(self._DIM):
           Ki  = self._get_Ki(i) + 1
           f1 = f_d1[Ki+0] - f_d1[Ki+1]
           f1/= 2.*self._step
           fder[i] = f1
           f2 = f_d1[Ki+0] + f_d1[Ki+1] - 2.*f_d1[0]
           f2/= self._step**2.0
           sder[i,i] = f2
           # calculate offdiagonal second derivatives
           for j in range(self._DIM):
             if i!=j:
               Kj = self._get_Ki(j) + 1
               if i>j: 
                  I = j
                  J = i
               else:
                  I = i
                  J = j
               Kij = self._get_Kij(I,J)
               f2 = f_d2[Kij ] + f_d2[Kij+1]     \
                  - f_d1[Ki+0] - f_d1[Kj+0]      \
                  - f_d1[Ki+1] - f_d1[Kj+1]      \
                  + 2.*f_d1[0]
               f2/= 2.*self._step**2.
               sder[i,j] = f2
       # symmetrize the second derivatives
       if symm:
          tl = numpy.tril(sder,k=-1)
          tu = numpy.triu(sder,k=+1).transpose()
          av =(tl + tu)/2.0
          sder = numpy.diag(sder.diagonal()) + av + (av.copy()).transpose()
       return fder, sder

def calc_tcf(f,t_max,dt,n=None):
    """Calculate the time autocorrelation function by using the Wiener-Khintchine theorem"""
    # evaluate FFT of f(t)
    v, gr, gi, v_max, v_res = ft_1d(f,t_max,dt=dt,n=n,algorithm='fft',cunit=None)
    # evaluate the |g(w)|**2
    g2 = gr*gr + gi*gi
    # evaluate IFFT of g2(w)
    N = numpy.where(v>v_max)[0][0]
    #t, cr, ci, v_max, v_res = ft_1d(g2[:N],v[N],dt=dt,n=n,algorithm='ifft',cunit=None)
    t, cr, ci, v_max, v_res = ft_1d(g2,v[-1],dt=dt,n=n,algorithm='ifft',cunit=None)
    return t, cr, ci

def ft_1d(f,t,dt,n=None,algorithm='fft',cunit=None):
    """\
------------------------------------------------------------------------------
Compute Discrete Fourier Transform of the time-domain signal f(t) 
measured through t seconds and sampled every dt seconds.

                      or

Compute Discrete Inverse Fourier Transform of the frequency-domain signal f(w) 
measured up to w Hz and sampled every dw Hz.
------------------------------------------------------------------------------

Usage:
v, gr, gi, v_max, v_res = ft(f,t,dt,n=None,algorithm='fft',cunit=None)

Input:
 - f    - 1d ndarray of time-domain signal points (real or complex)
 - t    - time of measurement [s]
 - dt   - sampling time [s]
 - n    - request extended number of points
 - algorithm - default: FFT of Cooley and Tukey (scales as np*np)
               other:   DFT - explicit implementation (scales as np*log_2(np))
 - cunit- return frequencies converted from Hz to another unit (detault is Hz)
         'ang' - angular frequency (multiply by 2pi)
         'cm-1'- wavenumber
          any number can be also provided. Then frequencies will be multiplied
          by it.
Output:
 - v    - ndarray of frequencies [Hz]
 - gr   - real part of Fourier spectrum
 - gi   - imaginary part of Fourier spectrum
 - v_max- maximum reasonable frequency [Hz]
 - v_res- resolution of Fourier spectrum [Hz]
 
Notes:
 - zero padding is automatic if <n> is switched and is larger 
   than len(<f>). 
 - If FFT algorithm is switched, len(<f>) is not a power of 2 and <n=None>
   the ft_1d function will zero-padd to reach the total number of points which
   is equal to the nearest power of 2. However, if the user requests specific
   number of points for FFT then <n> has to be a power of two.
 - the output Fourier spectrum in frequency domain should be taken up to the
   <v_max>. Note that the spectrum for higher frequencies should be rejected!
   To increase <v_max> decrease sampling time.
------------------------------------------------------------------------------
                                                 Last revision:  9 Oct 2014
"""
    if   algorithm.lower()=='fft' : ftfunc = fourier.ft.fft
    elif algorithm.lower()=='ifft': ftfunc = fourier.ft.fftinv
    elif algorithm.lower()=='dft' : ftfunc = fourier.ft.dft
    elif algorithm.lower()=='idft': ftfunc = fourier.ft.dftinv
    else: raise Exception('The algorighm %s is not implemented' % algorithm.lower())
    # determine the type of Fourier Transform
    if algorithm.lower().startswith('i'): inverse = True
    else:                                 inverse = False
    # check if the data type is complex or purely real
    if 'complex' in str(type(f[0])):
       f_real = f.real.copy()
       f_imag = f.imag.copy()
    else:
       f_real = f.copy()
       f_imag = None
    # check changing units
    uconv = 1.000
    if cunit is not None:
       try:
         if not inverse:
            if   cunit.lower() == 'ang' : uconv = 2.00*math.pi          
            elif cunit.lower() == 'cm-1': uconv = units.UNITS.HzToCmRec*2.00*math.pi
            elif cunit.lower() == 'hz'  : pass
         else:
            if   cunit.lower() == 'sec' : pass
       except TypeError:
         uconv = numpy.float64(cunit)
    #
    nf = len(f)
    ht = numpy.float(t)/(nf-1)
    # Fast Fourier Transform (Cooley-Tukey)
    if algorithm.lower()[-3:]=='fft':
       # prepare the data points
       N  = numpy.array([ 2**i for i in range(4,30) ], int)
       if n is None:                         
          if not nf in N:
             np = N[numpy.where(N>nf)][0]
             fr = numpy.zeros(np,numpy.float64)
             fi = numpy.zeros(np,numpy.float64)
             fr[:nf] = f_real.copy()
             if f_imag is not None: 
                fi[:nf] = f_imag.copy()
          else:
             np = nf
             fr = f_real.copy()
             if f_imag is not None:
                fi = f_imag.copy()
             else: 
                fi = numpy.zeros(np,numpy.float64)
          #
       else:
          nf = len(f)
          merror = " Invalid number of points requested!"
          assert n in N, merror
          merror = " The number of points <n=%i> is smaller than the size of input data: %i!" % (n,nf)
          assert n >= nf, merror
          np = n
          fr = numpy.zeros(np,numpy.float64)
          fi = numpy.zeros(np,numpy.float64)
          fr[:nf] = f_real.copy()
          if f_imag is not None:
             fi[:nf] = f_imag.copy()
       #
       m = int(numpy.log2(np))
       # do the FFT
       v_max = 1./(2.*ht)  # Hz
       v_res = 1./ (np*ht) # Hz
       v     = numpy.linspace(0,np,np) / (ht*np)  # Hz

       gr,gi = ftfunc(fr,fi,m) 
       gr/= math.sqrt(np)
       gi/= math.sqrt(np)
          
    # Direct Discrete Fourier Transform
    elif algorithm.lower()[-3:]=='dft':
       if n is not None: 
          message = " Requested number of points <n=%d> is smaller than data size <%d> so n is ignored" % (n,nf)
          if n<nf: 
             print message
             np = len(nf)
             fr = f_real.copy()
             if f_imag is not None:
                fi = f_imag.copy()
          else: 
             np = n
             fr = numpy.zeros(np,numpy.float64)
             fi = numpy.zeros(np,numpy.float64)
             fr[:nf] = f_real.copy()
             if f_imag is not None:
                fi[:nf] = f_imag.copy()
       else:
          np = nf
          fr = f_real.copy()
          fi = numpy.zeros(np,numpy.float64)
          if f_imag is not None:
             fi = f_imag.copy()
       #
       v_max = 1./(2.*ht)  # Hz
       v_res = 1./ (np*ht) # Hz
       v     = numpy.linspace(0,np,np) / (ht*np)  # Hz

       gr = numpy.zeros(np,numpy.float64)
       gi = numpy.zeros(np,numpy.float64)
       gr,gi = ftfunc(fr,fi,gr,gi) 
       gr/= math.sqrt(np)
       gi/= math.sqrt(np)

    # bad algorithm request
    else:
       error = " Not known algorithm=%s requested!" % algorithm.lower()
       raise ValueError, error
    #
    v *= uconv
    gr*= uconv   ; v_max*=uconv
    gi*= uconv   ; v_res*=uconv
    #
    return v, gr, gi, v_max, v_res

def resol(t,dt):
    """
 Provide the maximum frequency and resolution (in [cm-1])
 Usage:
    resol(t,dt), where t is sampling time and dt is spacing between samples 
                 (everything is vigen in [ps])
"""
    m = int(t/dt)
    df= 1./(dt*m) *1.0e+12  # Hz
    df*= units.UNITS.HzToCmRec * 2.*math.pi
    #
    fm= 1./(2.*dt) * 1.0e+12 * units.UNITS.HzToCmRec * 2.*math.pi
    print " Max: %10.1f [cm-1], Resolution: %10.1f [cm-1] " % (fm,df)
    return

def ParseLmocFromGamessEfpFile(efp_file):
    """Parse LMO centroinds from GAMESS *.efp file. Returns them in A.U."""
    # open the *.efp file and read the contents
    f = open(efp_file)
    text = f.read()
    f.close()

    # search for section with LMOCs
    sec = re.compile(' STOP', re.DOTALL)
    s = re.split(sec,text)

    lmoc_text = s[5]

    # extract LMOCs
    k = re_templates.re_real + '.*'
    templ = 'CT.* '+ 3*k + '\n'

    l = re.findall(templ, lmoc_text)

    LMOC = []
    for i in l:
        LMOC.append(i.split()[1:])
    LMOC = numpy.array(LMOC,numpy.float64)

    return LMOC

class MDOut(units.UNITS):
   """
 Represents MD output file containing basic information
 about the trajectory, i.e. time, temperatures and so on.
 
 Usage:
 
 d = MDOut(file,pkg='amber')    # amber is default
 d.read()                       # read the file
 d.write(time=True)             # write the reports. 
                                # time=True is default
                                # time=True: write in time (ps)
                                # time=False: write in time steps
 d.get_obs()                    # return MD observable dictionary

 Notes:

   1) Now only AMBER out files are supported
   2) Units saved and stored are in the units identical as in output
      file!
"""

   def __init__(self,file=None,pkg='amber'):
       self._init(pkg)
       if file is not None:
          self.__call__(file)
       return

   # --- P U B L I C ---
   
   def __call__(self,file):
       """open the file"""
       self._open(file)
       return

   def read(self):
       """read the MD output file"""
       if self.__pkg == 'amber':
          for obs,val in self.__obs.items():
              new_val = self._findall(obs)
              self.__obs[obs] = new_val
       else: 
          raise NotImplementedError
       return

   def write(self,time=True):
       """write the reports"""
       if time: x= self.__obs['time']
       else   : x= self.__obs['step']
       for obs,val in self.__obs.items():
           if obs!='time' and obs!='step':
              self._write_report(obs,x,val)
       return
 
   def get_obs(self):
       """return the dictionary with observables"""
       return self.__obs
  
   # --- P R O T E C T E D ---
   
   def _open(self,file):
       """read the text from MD output file"""
       if type(file) == '''<type 'str'>''':
          print 'fff'
          p = open(file)
          self.__text = p.read()
          p.close()
       else: 
          text = ''
          for f in file:
              p = open(f)
              text+= p.read()
              p.close()
          self.__text = text
       return

   def _findall(self,obs):
       """find all mathces for an observable"""
       p = self.__pattern[obs]
       m = re.findall(p, self.__text)
       if obs!='step':
          m = numpy.array(m,numpy.float64)
       else:
          m = numpy.array(m,int)
       return m

   def _write_report(self,obs,x,y):
       """write a report for a particular observable"""
       report = open('summary.'+obs.upper(),'w')
       title = '#\n' ; log = ''
       report.write(title)
       if len(x)!=len(y): 
          print " Not equal number of X and Y points for observable %s" % obs.upper()
       else:   
          for i in range(len(x)):
              log+= '%16.3f %16.6f\n' % (x[i],y[i])
          report.write(log)
       report.close()
       return
 
   def _update(self):
       """update the object"""
       # update the observable dictionary
       self.__obs = {'time'  : self.__time,
                     'step'  : self.__step,
                     'volume': self.__volume,
                     'pres'  : self.__pres,
                     'dens'  : self.__dens,
                     'temp'  : self.__temp,
                     'ekin'  : self.__ekin,
                     'epot'  : self.__epot,
                     'etot'  : self.__etot,
                     'escf'  : self.__escf,
                     'bond'  : self.__bond,
                     'angle' : self.__angle,
                     'dihed' : self.__dihed,
                     'ehbond': self.__ehbond,
                     '14el'  : self.__14el,
                     'elec'  : self.__elec,
                     'vdw'   : self.__vdw,
                     'virial': self.__virial, }
       return

   def _init(self,pkg):
       """create the namespace of memorials"""
       self.__pkg = pkg
       # MD observables
       self.__time    = None ;  self.__bond    = None
       self.__volume  = None ;  self.__angle   = None
       self.__pres    = None ;  self.__dihed   = None
       self.__dens    = None ;  self.__ehbond  = None
       self.__temp    = None ;  self.__14el    = None
       self.__ekin    = None ;  self.__elec    = None
       self.__epot    = None ;  self.__vdw     = None
       self.__etot    = None ;  self.__virial  = None
       self.__escf    = None ;  self.__step    = None
       self.__restr   = None
       # search patterns
       _real = r'\s*(-?\d*\.\d*)'
       _int  = r'\s*(\d*)'
       self.__pattern = { 'time'  : r'TIME\(PS\) ='  +_real,
                          'step'  : r'NSTEP ='     +_int,
                          'volume': r'VOLUME     ='+_real,
                          'pres'  : r'PRESS ='     +_real,
                          'dens'  : r'Density    ='+_real,
                          'temp'  : r'TEMP\(K\) ='   +_real,
                          'ekin'  : r'EKtot   ='   +_real,
                          'epot'  : r'EPtot      ='+_real,
                          'etot'  : r'Etot   ='    +_real,
                          'escf'  : r'ESCF='       +_real,
                          'bond'  : r'BOND   ='    +_real,
                          'angle' : r'ANGLE   ='   +_real,
                          'dihed' : r'DIHED      ='+_real,
                          'ehbond': r'EHBOND  ='   +_real,
                          '14el'  : r'1-4 EEL ='   +_real,
                          'elec'  : r'EELEC  ='    +_real,
                          'vdw'   : r'VDWAALS    ='+_real,
                          'virial': r'VIRIAL  ='   +_real, 
                          'restr' : r'RESTRAINT  ='+_real }
       #
       self._update()  
       return


def text_to_list(text,delimiter=None,dtype=int,dt=1):
    """
Transform one text line into a list. Returns list of integers.

Syntax: 
       text_to_list(text,delimiter=None)

Usage:

1) text = '     2 3 4 5 6    4    '

text_to_list(text) = [2 3 4 5 6 4]

2) text = ' 2-4'

text_to_list(text) = [2 3 4]

3) text = ' 3, 56, 6, -1, 33  '

text_to_list(text,delimiter=',') = [3 56 6 -1 33] 

4) text = ' -0.44, 0.122, 23.4544, -335.2   '
text_to_list(text,delimiter=',',dtype=float) = [-0.44, 0.122, 23.4544, -335.2]

"""
    if (('-' in text and '.' not in text) or ('-' in text and ',' not in text)):
       start,end = text.split('-')
       rang = numpy.arange(dtype(start),dtype(end)+1,dt)
    else:
       rang = text.split(delimiter)
       rang = map(dtype,rang)
    return numpy.array(rang)

def dihedral(A,unit='radian'):
    """Compute dihedral angle n1-n2-n3-n4. 
Usage: 
dihedral([4x3 array],<unit='radian' or 'deg'>).
Provide real atom numbers. The dihedral evaluated by this code gave opposite signs 
as compared with MOLDEN for a test NMA molecule"""
    P1 = A[1]-A[0]
    P2 = A[2]-A[1]
    P3 = A[3]-A[2]
    N1 = numpy.cross(P1,P2)
    N2 = numpy.cross(P2,P3)
    N1/= math.sqrt(sum(N1*N1))
    N2/= math.sqrt(sum(N2*N2))
    P2/= math.sqrt(sum(P2*P2))
    M1 = numpy.cross(N1,P2)
    x = sum(N1*N2)
    y = sum(M1*N2)
    angle = numpy.arctan2(y,x)
    conv = 180./math.pi
    if unit=='deg':angle*=conv
    return angle

def numerov1(q,s,x0,y0,y1,npoints,step,qarg={},sarg={}):
    """
Integrates the equation of the form:

        y''(x) + q(x)u(x) = s(x)

using Numerov algorithm. Here it is assumed that 
user provides explicitly the values of the two first
poitns of the solution which is equivalent to giving
y(0) and y'(0).

Usage:
q    - function object for q(x) or a float64 number
qarg - dictionary of optional arguments (other than x) for q
s    - function object for s(x) or a float64 number
sarg - dictionary of optional arguments (other than x) for s 
x0   - starting point of integration
y0,y1- first and second quess points for the solution
npoints - numper of points for the numerical solution
          (including the first two points)
step - delta x for numerical integration
"""
    # initialize outputs
    x = numpy.float64(numpy.linspace(x0,x0+(npoints-1)*step,npoints))
    y = numpy.zeros(npoints,numpy.float64)
    y[0] = y0; y[1] = y1
    #
    H = step*step/12.0
    if   isinstance(q,numpy.float64): Q = q*numpy.ones(npoints,numpy.float64)
    elif isinstance(q,numpy.ndarray): Q = q
    else:                             Q = q(x,**qarg)
    #
    if   isinstance(s,numpy.float64): S = s*numpy.ones(npoints,numpy.float64)
    elif isinstance(s,numpy.ndarray): S = s
    else:                             S = s(x,**sarg)
    # calculate all the coefficients
    ci1 = 1.0 +      H * Q
    ci  = 2.0 - 10.0*H * Q
    d   = H * (S[2:] + 10.0*S[1:-1] + S[:-2])
    # proceed the integration
    for i in xrange(2,npoints):
        y[i] = (d[i-2] + y[i-1]*ci[i-1] - y[i-2]*ci1[i-2]) / ci1[i]
    return x,y

def numerov2(q,x0,y0,y1,npoints,step,**qarg):
    """
Integrates the equation of the form:

        y''(x) + q(x)u(x) = 0

using Numerov algorithm. Here it is assumed that 
user provides explicitly the values of the two first
poitns of the solution which is equivalent to giving
y(0) and y'(0).

Usage:
q    - function object for q(x) or a float64 number
qarg - optional arguments passed to q (other than x)
x0   - starting point of integration
y0,y1- first and second quess points for the solution
npoints - numper of points for the numerical solution
          (including the first two points)
step - delta x for numerical integration
"""
    # initialize outputs
    x = numpy.float64(numpy.linspace(x0,x0+(npoints-1)*step,npoints))
    y = numpy.zeros(npoints,numpy.float64)
    y[0] = y0; y[1] = y1
    #
    H = step*step/12.0
    if   isinstance(q,numpy.float64): Q = q*numpy.ones(npoints,numpy.float64)
    elif isinstance(q,numpy.ndarray): Q = q
    else:                             Q = q(x,**qarg)
    # calculate all the coefficients
    ci1 = 1.0 +      H * Q
    ci  = 2.0 - 10.0*H * Q
    # proceed the integration
    for i in xrange(2,npoints):
        y[i] = (y[i-1]*ci[i-1] - y[i-2]*ci1[i-2]) / ci1[i]
    return x,y

def fder5pt(data,step):
    """Calculate first derivatives of f(x) wrt x using 5-point central Stencil method.
Uses numpy.slice objects to increase the speed instead of using slow 'for' python loop.
Note: The i-th derivative in <deriv> refers to i+2-th data point in <data>"""
    u = len(data) - 2
    deriv = (data[0:u-2]-data[4:u+2]) - 8.*(data[1:u-1]-data[3:u+1])
    deriv/= 12.*step
    return deriv

def secant(f,x,delta=1.e-6,max_iter=1000,**args):
    """Calculate the root of function f(x [,args]) wrt x"""
    x_old = x
    fx = f(x_old,**args)
    x_new = x_old - delta * fx / (fx-f(x-delta,**args))
    n_iter = 1
    while math.fabs(x_new-x_old)>delta:
          x_old = x_new
          fx = f(x_old,**args)
          x_new = x_old - delta * fx / (fx-f(x_old-delta,**args))
          n_iter += 1
          if n_iter > max_iter: break
    return x_new

def simpson(data,step):
    """Integrate f(x) from a to b using Simpson rule. <data> is a ndarray of function values
evaluated from x=a to x=b every dx where dx is assumed to be constant. Use numpy.slice
objects to increase the speed instead of using slow 'for' python loop"""
    u = len(data) - 2
    S = (data[0:u-2] + 4.*data[1:u-1] + data[2:u])*step/3.
    s = S[::2].sum()
    # in the case of even number of points
    if (u+1)%2: s+= (-data[-3] + 8.*data[-2] + 5.*data[-1])*step/12.
    return s

def simpson_nonuniform(data,args):
    """Integrate f(x) from a to b using Simpson rule. <data> is a ndarray of function values
evaluated from x=a to x=b every dx where dx is allowed to be not constant. Use numpy.slice
objects to increase the speed instead of using slow 'for' python loop"""
    aerror = 'The lengths of arguments and function values are not equal! (%d,%d)' % (len(args),len(data))
    assert len(args)==len(data), aerror
    u = len(args) - 1
    H = args[1:u+1] - args[:u]
    u-= 1
    hi  = H[1:u+1]
    him1= H[:u]
    hi6 = 6.*hi
    #
    a = (2.*hi*hi + hi*him1 - him1*him1)/hi6
    b = (hi+him1)**3 / (hi6*him1)
    c = (-hi*hi + hi*him1 + 2*him1*him1)/hi6
    #
    S = data[2:u+2] * a + data[1:u+1] * b + data[0:u  ] * c
    s = S[::2].sum()
    # in the case of even number of points
    if (u+1)%2:
        hnm1 = H[-1]; hnm2 = H[-2]
        d1   = data[-1]; d2 = data[-2]; d3 = data[-3]
        a = hnm1/6. * (3. - hnm1/(hnm1+hnm2))
        b = hnm1/6. * (3. - hnm1/hnm2)
        c =-hnm1/6. * hnm1*hnm1/(hnm2*(hnm1+hnm2))
        s+= d1 * a + d2 * b + d3 * c
    return s

class RungeKutta:
   """
 -----------------------------------------------------------------------------
             Numerical Solver of Ordinary Differential Equations
                        *** Runge-Kutta Method ***
                Bartosz Błasiak, blasiak.bartosz@gmail.com
 -----------------------------------------------------------------------------

 Solves numerically the set of n first-order ordinary differential equations of
 the form:
 
 
    dy_1(t)
    ------- = f(y_1(t), y_2(t), ..., f_n(t)) = f({y_i(t)})
      dt 

                    ···
    
    dy_n(t)
    ------- = f(y_1(t), y_2(t), ..., f_n(t)) = f({y_i(t)})
      dt 

 RungeKutta class is a purely-Python vectorized implementation of Kunge-Kutta 
 method.
 
 Usage:
 
    solver = RungeKutta()
    solver.set(func,tau,init)
    solver.eval(n_pass,**kwargs)
    result = solver()  or  result = solver.get()
 
 Inputs:
 
    func = function returning an ndarray of f({y_i(t)}) functions 
           for each dy_n(t)/dt (see below for an example)
    tau  = step of integration
    init = initial conditions
    **kwargs = additional arguments of func
    
 Notes:
 
    o func and init are ndarrays of length n
    
 Example:
 --------
 
 Here we consider a driven pendulum. We want to obtain the functions of angle
 theta(t) and angular frequency dtheta(t)/dt as a function of time. We denote
 y_1(t) = theta(t) and y_2(t) = dtheta(t)/dt. w_0 and b are the frequency and
 amplitude of a driving force whereas q is a damping coefficient (q and b >0).
 The set of two coupled 1-st order differential equations is as follows:
 
    dy_1(t)
    ------- = y_2(t)                                                (1)
      dt 

    dy_2(t)
    ------- = -q y_2(t) - sin(y_1(t)) + b cos(w_0 t)                (2)
      dt 

  Note that these can be written as a one second-order differential equation.
  Now, we want to construct <func>. Below there is an example of solving this
  problem using RungeKutta class method for the case when initial conditions
  are theta(0) = 2.3 radians and omega(0) = 0.3 radians/second. What is the be
  haviour of the pendulum? Try to change b value to 1.5 and 0.8 and examine 
  the differences:
  
    import utilities
    from numpy import array, linspace, sin, cos, sqrt, pi
                                          
    # define inputs
    gfc= lambda y, tau, q, b, w_0: array([   y[1]                              ,
                                          -q*y[1] - sin(y[0]) + b*cos(w_0*tau) ])
                                          
    steps = 10000; dt = 3*pi/100.; w_0 = 2./3.; q = 0.5; b = 3.5
    init  = [ 2.3, 0.3 ]
    
    # create RungeKutta instance and solve the problem
    pendulum = utilities.RungeKutta()
    pendulum.set(gfc,tau=dt,init=init)
    pendulum.eval(steps,w_0=w_0,q=q,b=b)
    t = linspace(0.,steps*dt,steps)
    r = pendulum()
    y_1 = r[:,0]; y_2 = r[:,1]
    
    # make frequency-domain power spectrum signal
    v, gr, gi, v_max, v_res = utilities.ft_1d(y_1,steps,dt,algorithm='FFT',
                                                           cunit='Hz')
    s = sqrt(gi**2. + gr**2.)
    print " Resolution of spectrum      : %10.6f" % v_res
    print " Maximal reasonable frequency: %10.2f" % v_max
    
    # plot the happy results
    import pylab as pl
    f = pl.figure()
    ax1 = f.add_subplot(311)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('$\\Theta(t)$ [rad]')
    ax1.plot(t,y_1,'b-')
    ax2 = f.add_subplot(312)
    ax2.set_xlabel('$\\Theta(t)$ [rad]')
    ax2.set_ylabel('$\\omega(t)$ [rad/s]')
    ax2.plot(y_1,y_2,'g-')
    ax3 = f.add_subplot(313)
    ax3.set_xlabel('$\\upsilon$ [Hz]')
    ax3.set_ylabel('Intensity')
    ax3.plot(v,s,'r-')
    ax3.set_xlim([0.0,v_max])
    pl.show(block=True)

    
  Note that:
    o keyword for step in <func> has to be exactly <tau>
    o the function <func> has to return ndarray of functions (NOT list!)
    o init can be tuple, list or array
    o steps is the number of steps including the initial condition
 -----------------------------------------------------------------------------
                                               Last revision: 10 Oct 2014
"""
   def __init__(self):
       """happy birthday!"""
       self.__gfc    = None
       self.__ndim   = None
       self.__result = None
       pass

   def __call__(self):
       """return the nd-array of time-dependent variables"""
       return self.get()

   def set(self,func,tau,init):
       """set an array of g-functions!"""
       self.__gfc = func
       self.__tau = tau
       self.set_init(init)
       return

   def set_init(self,init):
       """set the initial variables"""
       self.__init = numpy.array(init,numpy.float64)
       self.__ndim = len(init)
       return

   def get(self):
       """return the nd-array of time-dependent variables"""
       return self.__result

   def eval(self,n_pass,**kwargs):
       """evaluate the dynamic properties as time goes"""
       y = numpy.zeros((n_pass,self.__ndim),numpy.float64)
       # initial conditions
       y[0] = self.__init
       # proceed with time
       for i in xrange(1,n_pass):
           ti = self.__tau*i
           a = y[i-1]
           c1 = self.__gfc( a      , tau=ti               ,**kwargs) * self.__tau
           c2 = self.__gfc( a+c1/2., tau=ti+self.__tau/2. ,**kwargs) * self.__tau
           c3 = self.__gfc( a+c2/2., tau=ti+self.__tau/2. ,**kwargs) * self.__tau
           c4 = self.__gfc( a+c3   , tau=ti+self.__tau    ,**kwargs) * self.__tau
           y[i] = y[i-1] + 1./6. * (c1+2.*(c2+c3)+c4)
       self.__result = y
       return

if import_matplotlib:
   class Graph:
       "graph environment"
       fs = 20
       font0 = matplotlib.font_manager.FontProperties()
       font1 = font0.copy()
       font1.set_family('sans-serif')
       font1.set_weight('bold')
       marker = 'o'
    
class QMFile(units.UNITS):
   """
---------------------------------------------------------
           Represents the following file formats: 
                   XYZ, DMA, LOG and FCHK.

These files contain molecular structure and other quantum 
chemistry information. These formats mean:

1) XYZ - structure file format. Contains number of atoms,
         comment, list of atoms and cartesian coordinates
2) DMA - Solvshift Distributed Multipole Analysis format.
         Contains atoms, atomic cartesian coordinates, di
         stributed multipoles and their origins.
3) FCHK- Gaussian formatted checkpoint file. Contains ato
         mic structure, wavefunction of a molecule and ma
         ny other useful information.
---------------------------------------------------------
                                             May 5 2014
Usage:

f = QMFile()                                 # initialize object

atoms, xyz = f(file,format,**kwargs)         # get atomic symbols and coordinates
                                             # directly by invoking QMFile class 
                                             # instance as a function

f.open(file,format='fchk',**kwargs)          # open file using initialized object

str   = f.get_pos(**kwargs)                  # return atomic coordinates 
atoms = f.get_atoms()                        # return atom list
mol   = f.get_mol()                          # return Molecule object
dma   = f.get_dma()                          # return DMA object
misc  = f.get_misc()                         # return miscellanea

Notes:

 1) if format is not provided, it will be determined automatiacally
    from file name extension

 2) Now, only structural information is supported (coordinates and atomic lists
    only)

 3) Additional options for **kwargs:

    - units = 'angs' or 'bohr'. Default is 'bohr'. This option specify
            final units of memorials being returned (structural information). 
    
 4) all memorials are stored in A.U.!
---------------------------------------------------------------------------------
                                                                      @Globulion
"""

   def __init__(self,file=None,format=None,**kwargs):
       self._init()
       if file is not None: self.__call__(file,format,**kwargs)
       return

   def open(self,file,format=None,
                 units='Angstrom',name='Unnamed molecule',
                 mult=1,charge=0,method='HF',basis='3-21G',
                 mol=True,freq=False,anh=False,oniom=False,
                 pol=False):
       """open file"""
       self.__file_name = file
       self.__format = format
       if format is None:
          if   file.endswith('.xyz'):   self._open_xyz(file,units,name,mult,charge,method,basis,mol)
          elif file.endswith('.dma'):   self._open_dma(file)
          elif file.endswith('.fchk'):  self._open_fchk(file,units,name,method,basis)
          elif (file.endswith('.g09') or file.endswith('.log')):   self._open_g09(file,freq,anh,oniom,pol)
       else:
          self.__format_dict[format](file,units,name,mult,charge,method,basis,mol)
       return

   def close(self):
       """close the file and clean the memory"""
       # close the appropriate file object
       if self.__file_obj is not None: self.__file_obj.close()
       # clean the memory
       self._init()
       return
   
   def set_pos(self,pos,units='bohr'):
       """set the position array"""
       assert len(pos) == len(self.__atoms), "The number of atoms is not matching to the coordinates (%d,%d)" % (len(self.__atoms),len(pos))
       if (units.lower()).startswith('angs'): self.__pos = pos * self.AngstromToBohr
       else: self.__pos = pos
       return 

   def set_atoms(self,pos,atoms):
       """set the atoms and positions"""
       self.__pos = pos
       self.__atoms = atoms
       return
  
   def set_charges(self, chg):
       """set the charges array"""
       if   self.__pos.shape[1]==3:
            chg = numpy.array(chg).reshape(1,len(chg))
            self.__pos = numpy.hstack((self.__pos, chg.T))
       elif self.__pos.shape[1]==4:
            self.__pos[:,3] = chg
       return

   def set_misc(self, text):
       """set the miscellaneum: text as a description"""
       self.__misc = text
       return
 
   def write(self, name='default', pkg=None, template=None, delim='@', ext=None, overwrite=False, **kwargs):
       """
 Write the structure into the XYZ file or input file for GAMESS, Gaussian or Coulomb.

 Usage:

   QMFile_instance.write(self, name='default', pkg=None, template=None, delim='@', ext=None, **kwargs)

 Notes:

   template - is a usual string with a template argument for string.Template
   delim    - is a delimiter used by string.Template
   **kwargs are the keys for template substitutions (eg CHK when @CHK with delim='@')
   *COORD keyword needs always to be present in template string object where * is the chosen delimiter
        (eg. @COORD, $COORD etc)
   ext      - optional extension for file (default are .xyz or .inp)

 Examplary template:

         templ = '''\
         %mem=@MEM
         %nprocshared=@NCPUS
         %chk=@CHK
         
         # RHF/6-311++G** nosymm scf(conver=10,xqc) 6D 10F
         
         test
         
         0 1
         @COORD
         '''

 Warning: keyword pkg is not working yet (so if wanted to make GAMESS inputs it will not add atomic numbers
          and you will have to add them manually or using external script.
"""
       if not overwrite:
          if ext is None:                               
             if template is None: ext = '.xyz' 
             else               : ext = '.inp'
          if not name.endswith(ext):  name = name + ext
          f = open(name,'w')
       else            : f = name

       if template is None:
          f.write('%d\n' % len(self.__atoms))
          if self.__misc is None: f.write('\n')
          else:                   f.write('%s' % self.__misc)

       log = ''
       if self.__pos.shape[1]==3:
          for i in range(len(self.__atoms)):
              #log+= '%2s ' % self.__atoms[i]
              #print self.__atoms[i].symbol
              log+= '%2s ' % self.__atoms[i]#.symbol
              log+= '%14.6f %14.6f %14.6f\n' % tuple(self.__pos[i]*self.BohrToAngstrom)
       elif self.__pos.shape[1]==4:
          for i in range(len(self.__atoms)):
              log+= '%2s ' % self.__atoms[i]#.symbol
              log+= '%14.6f %14.6f %14.6f' % tuple(self.__pos[i,:3]*self.BohrToAngstrom)
              log+= '%14.6f\n' % self.__pos[i,3]

       if template is None: f.write(log)
       else:
          class LocalTemplate(string.Template): delimiter=delim
          template = LocalTemplate(template)
          txt = template.substitute(COORD=log, **kwargs)
          f.write(txt)
          
       if not overwrite: f.close()
       return
   
   def insert(self,xyz=None,id=0,atoms=None,mol=None):
       """
Insert atoms into the structure.

Usage:

xyz   - array of coordinates in Bohr (and charges, optional)
id    - insert xyz after (id)th atom. Default is 0 (insert at the beginnig)
atoms - list of atomic symbols. Default is None (dummy atoms, 'X')
mol   - molecule object that will be inserted (only coordinates and atoms)
"""
       if id==-1: id = len(self.__atoms)
       N,M = self.__pos.shape                              

       if mol is not None:
          xyz = mol.get_pos()
          atoms = mol.get_atoms()
          n,m = xyz.shape; I=id
       else:
          n,m = xyz.shape
          # case where there is mismatch in number of columns
          if m+1==M:
             xew = numpy.zeros((n,m+1),numpy.float64)
             xew[:,:3] = xyz
             xyz = xew
          I = id
       # insert the structure
       new = numpy.zeros((N+n,M),numpy.float64)
       new[:I] = self.__pos[:I]
       new[I:I+n] = xyz
       new[I+n:] = self.__pos[I:]
       # insert atom symbols
       if atoms is None:
          atoms = ['X' for i in range(n)]
       for i in range(n):
           self.__atoms.insert(id+i,atoms[i])
       # save
       self.__pos = new
       return
   
   def get_all(self):
       """return all memorials in a tuple: (atoms, pos, mol, dma, misc)"""
       return self.__atoms, self.__pos, self.__mol, self.__dma, self.__misc

   def get_pos(self):
       """return positions (and charges, if present)"""
       return self.__pos
   
   def get_atoms(self):
       """return atoms"""
       return self.__atoms

   def get_coord_list(self):
       """return coord list containing atom symbols and position vectors"""
       return self.__coord_list

   def get_mol(self):
       """return Molecule object"""
       return self.__mol

   def get_bfs(self):
       """return Basis Set object"""
       return self.__mol.get_bfs()

   def get_dma(self):
       """return DMA object"""
       return self.__dma

   def get_dmat(self):
       """Return RHF density matrix"""
       run = PyQuante.SCF(self.__mol,bfs=self.__mol.get_bfs())
       run.iterate()
       return run.dmat

   def get_misc(self):
       """return misc (e.g.: comment second line in XYZ file)"""
       return self.__misc
   
   def get_file_name(self):
       """return file name"""
       return self.__file_name
   
   def get_freq(self):
       """return frequencies, reduced masses, force constants, IR intensities and eigenvectors"""
       return self.__freq, self.__redmss, self.__forcec, self.__irints, self.__eigvec
   
   def translate(self,transl):
       """translate tensors by <transl> cartesian displacement vector"""
       if self.__pos   is not None: self.__pos   += transl
       return
   
   def rotate(self,rot):
       """rotate the tensors by <rot> unitary matrix"""
       # transform the atomic position static and dynamic information
       if self.__pos   is not None:
          self.__pos    = numpy.dot(self.__pos , rot)
       if self.__eigvec  is not None:
          print "Eigenvectors were not rotated!"
       #   self.__eigvec   = dot(self.__lvec, rot)
       if self.__pol is not None:
          self.__pol = numpy.dot(numpy.transpose(rot),numpy.dot(self.__pol,rot))
       return
   
   def scale_freq(self,s):
       """scale the frequencies using a scaling factor"""
       self.__freq*= s
       return
   
   def get_spectrum(self,func='g',w=10.0,n=10000, sc=1.00):
       """return the spectrum"""
       x, s = self._func(func,w,n,sc)
       return x, s
   
   def get_oniom_perc(self):
       """return percentages of vibration between model and real systems in ONIOM calculation"""
       return self.__modelP, self.__realP
   
   def get_pol(self):
       """return molecular total polarizability"""
       return self.__pol
   
   def set_pol(self,pol):
       """set the polarizability"""
       self.__pol = pol
       return
   
   def _func(self,func,w,n,sc):
       """phenomenological line shape function for FTIR spectra"""
       x = numpy.linspace(self.__freq.min(),self.__freq.max(),n) * sc
       s = numpy.zeros(len(x),numpy.float64)
       # Gaussian function
       if func=='g':
          for i in range(self.__nmodes):
              s+= self.__irints[i] * exp( -4.*log(2.) * (x-self.__freq[i])**2/w**2 )
          s *= w * numpy.sqrt(math.pi/4. * log(2.))
       # Lorenzian function
       if func=='l':
          for i in range(self.__nmodes):
              s+= self.__irints[i] / ( 4.0 * (x-self.__freq[i])**2 + w**2 )
          s *= w * 2. / math.pi
       return x, s
   
   def __call__(self,file,format=None,**kwargs):
       """open file and return atom list (as symbols) and atomic coordinates (as ndarray)"""
       self.open(file,format,**kwargs)
       return self.__atoms, self.__pos

   def __repr__(self):
       """print me!"""
       pos_q, mol_q, dma_q, pol_q = 'No','No','No','No'
       if self.__pos is not None: pos_q = 'Yes'
       if self.__mol is not None: mol_q = 'Yes'
       if self.__dma is not None: dma_q = 'Yes'
       if self.__pol is not None: pol_q = 'Yes'
       log = '\n'
       log+= ' File: %s\n' % self.__file_name
       log+= '   * Atoms  : %i\n' % len(self.__atoms)
       log+= '   * Coord  : %s\n' % pos_q
       log+= '   * Mol    : %s\n' % mol_q
       log+= '   * DMA    : %s\n' % dma_q
       log+= '   * POL    : %s\n' % pol_q
       log+= '   * Misc   : %s\n' % str(self.__misc)
       return str(log)

   # protected

   def _init(self):
       """create the namespace of variables"""
       self.__pos       = None
       self.__atoms     = None
       self.__mol       = None
       self.__dma       = None
       self.__misc      = None
       self.__file_name = None
       self.__file_obj  = None
       self.__format    = None
       self.__pol       = None
       self.__eigvec    = None
       self.__units     = 'au'
       self.__defaults  = {'name'  :'Unnamed molecule',
                           'mult'  :1,
                           'charge':0,
                           'method':'HF',
                           'basis' :'3-21G'}
       self.__format_dict = {'xyz' :self._open_xyz,
                             'dma' :self._open_dma,
                             'fchk':self._open_fchk,
                             'log' :self._open_g09,
                             'gms' :self._open_gms,}
       return

   def _make_mol(self,name,coord,mult,charge,method,basis):
       """make the PyQuante.Molecule object"""
       Coords = []
       for i in range(len(coord)):
           atom  = (uAtom(coord[i][0]).atno, (coord[i][1], 
                                              coord[i][2],
                                              coord[i][3]) )
           Coords.append(atom)
       Mol = PyQuante.Molecule(name,Coords,units='Bohr',
                               multiplicity=mult,charge=charge,
                               basis=basis,method=method)
       return Mol

   def pop(self, idx, mol=True):
       """Remove atom from a molecule and create molecule object if specified - this method is NOT FINISHED !!! AND WORKS ONLY FOR GENERATING SHRINKED MOL OBJECT!!!!!!!!!!!!!!!!!!!!!!"""
       self.__coord_list.pop(idx) 
       if mol: self.__mol = self._make_mol('happy buddy', self.__coord_list, mult=1, charge=0, method='HF', basis='3-21G')
       return

   def _open_xyz(self,file,units,name,mult,charge,method,basis,mol=True):
       """open xyz file"""
       self.__file_obj = open(file)
       data = self.__file_obj.readlines()
       n_atoms = int(data[0])
       misc = data[1]
       data.pop(0);data.pop(0)
       coord = []
       atoms = []
       for i in range(n_atoms):
           line_spl = data[i].split() 
           coord.append(line_spl[:])   # it was [:4] instead of [:]
           atoms.append(line_spl[0])
           coord[i][1:] = map(numpy.float64,coord[i][1:])
           if units.lower()=='angstrom':
              for j in range(3):
                  coord[i][j+1]*= self.AngstromToBohr
            
       data = [map(numpy.float64,[x for x in coord[y][1:]]) \
                                    for y in range( len(coord))]
       data = numpy.array(data,dtype=numpy.float64)

       Coords = []
       for i in range(n_atoms):
           atom  = (uAtom(coord[i][0]).atno, (coord[i][1], 
                                              coord[i][2],
                                              coord[i][3]) )
           Coords.append(atom)
       Mol = None
       if mol:
          Mol = PyQuante.Molecule(name,Coords,units='Bohr',
                                  multiplicity=mult,charge=charge,
                                  basis=basis,method=method)
   
       self.__mol = Mol
       self.__pos = data
       self.__atoms = atoms
       self.__misc = misc
       self.__coord_list = coord
       return

   def _open_fchk(self,file,units,name,method,basis):
       """open fchk file"""
       file = open(file)
       self.__file_obj = file
       line = file.readline()
       g = lambda n,m: n/m+bool(n%m)

       # search for charge and multiplicity
       querry = "Charge"
       while True:
           if querry in line: break
           line = file.readline()
       charge = int(line.split()[-1])

       querry = "Multiplicity"
       while True:
           if querry in line: break
           line = file.readline()
       multiplicity = int(line.split()[-1])
       
       # search for atomic numbers
       querry = "Atomic numbers"
       while True:
           if querry in line: break
           line = file.readline()       
       n_atoms = int(line.split()[-1])
       line = file.readline()
       
       atnos = []
       for i in xrange(g(n_atoms,6)):
           atnos += line.split()
           line = file.readline()
           
       atnos = numpy.array(atnos,dtype=int)
       
       # search for atomic coordinates       
       querry = "Current cartesian coordinates"
       while True:
           if querry in line: break
           line = file.readline()
       N = int(line.split()[-1])
       line = file.readline()
       
       coord = []
       for i in xrange(g(N,5)):
           coord += line.split()
           line = file.readline()
           
       coord = numpy.array(coord,dtype=numpy.float64).reshape(n_atoms,3)

       # create Molecule object
       Coords = []
       for i in range(n_atoms):
           atom  = (atnos[i], (coord[i][0],
                               coord[i][1],
                               coord[i][2]) )
           Coords.append(atom)
       Mol = PyQuante.Molecule(name,Coords,units='Bohr',
                               multiplicity=multiplicity,charge=charge,
                               basis=basis,method=method)
   
       self.__mol = Mol                        
       self.__pos = coord
       self.__atoms = atnos
       return

   def _open_dma(self,file,**kwargs):
       """open dma file"""
       dma = dma.DMA(file=file)
       self.__dma = dma
       self.__pos = dma.get_pos()
       return
   
   def _open_gms(self,file):
       """open GAMESS log file (not working yet)"""
       file = open(file)
       self.__file_obj = file
       line = file.readline()
       raise NotImplementedError
   
   def _open_g09(self,file,freq,anh,oniom,pol):
       """open Gaussian G09RevD.01 log file. Attention: #P (extended printout) and NOSYMM keywords are necessary!"""
       file = open(file)
       self.__file_obj = file
       line = file.readline()
       ### search for atomic symbols, charge and multiplicity
       atoms = []
       querry1 = " symbolic z-matrix:"
       querry2a= " charge = "
       querry2b= " multiplicity = "
       while True:
           if (querry1 in line.lower() or querry2a in line.lower() and querry2b in line.lower()): break
           line = file.readline()
       if line.lower().startswith(' charge ='):
          chg = numpy.float64(line.split()[2])
          mult= numpy.float64(line.split()[5])
          line = file.readline()
       else:
          line = file.readline()
          chg = numpy.float64(line.split()[2])
          mult= numpy.float64(line.split()[5])
       line = file.readline()
       if oniom: 
          line = file.readline()
          line = file.readline()
       while len(line)>2:
           atom = uAtom(line.split()[0])
           atoms.append(atom)
           line = file.readline()
       ### search for atomic positions
       pos = []
       querry = "Input orientation:"
       while True:
           if querry in line: break
           line = file.readline()
       for i in range(4): line = file.readline()
       for i in range(len(atoms)):
           line = file.readline()
           xyz  = numpy.array(line.split()[3:],numpy.float64)
           pos.append(xyz)
       pos = numpy.array(pos,numpy.float64) * self.AngstromToBohr
       
       self.__nmodes = len(atoms) * 3 - 6
       self.__n3 = len(atoms) * 3
       ### search for polarizability
       if pol:
          querry = " Exact polarizability:"
          while True:
            if querry in line: break
            line = file.readline()
          pol = numpy.array(line.split()[-6:],numpy.float64)
          self.__pol = numpy.zeros((3,3),numpy.float64)
          self.__pol[0,0] = pol[0]
          self.__pol[0,1] = pol[1]
          self.__pol[1,0] = pol[1]
          self.__pol[1,1] = pol[2]
          self.__pol[0,2] = pol[3]
          self.__pol[2,0] = pol[3]
          self.__pol[1,2] = pol[4]
          self.__pol[2,1] = pol[4]
          self.__pol[2,2] = pol[5]
          
       ### search for frequencies, reduced masses, foce constants, IR intensities 
       if freq:
          nmodes = len(atoms) * 3 - 6
          querry = " Frequencies --- "
          n = nmodes
          while True: 
            if querry in line: break 
            line = file.readline()
          redmss = numpy.zeros(n,numpy.float64)
          freqs  = numpy.zeros(n,numpy.float64)
          irints = numpy.zeros(n,numpy.float64)
          forcec = numpy.zeros(n,numpy.float64)
          eigvec = numpy.zeros((self.__n3,n),numpy.float64)
          modelP = numpy.zeros(n,numpy.float64)
          realP  = numpy.zeros(n,numpy.float64)
          for j in range( n/5+bool(n%5) ):
              # frequencies
              freqs[(j*5):j*5+self.__dupa(j)] =\
              [ numpy.float64(line.replace('D','E').split()[-self.__dupa(j):][x])\
                                              for x in range(self.__dupa(j)) ]
              # reduced masses
              line = file.readline()
              redmss[(j*5):j*5+self.__dupa(j)] =\
              [ numpy.float64(line.replace('D','E').split()[-self.__dupa(j):][x])\
                                              for x in range(self.__dupa(j)) ]
              # force constants
              line = file.readline()
              forcec[(j*5):j*5+self.__dupa(j)] =\
              [ numpy.float64(line.replace('D','E').split()[-self.__dupa(j):][x])\
                                              for x in range(self.__dupa(j)) ]
              if oniom: 
                 # Percent ModelSys
                 line = file.readline()
                 modelP[(j*5):j*5+self.__dupa(j)] =\
                    [ numpy.float64(line.replace('D','E').split()[-self.__dupa(j):][x])\
                                                    for x in range(self.__dupa(j)) ]
                 # Percent RealSys
                 line = file.readline()
                 realP [(j*5):j*5+self.__dupa(j)] =\
                    [ numpy.float64(line.replace('D','E').split()[-self.__dupa(j):][x])\
                                                    for x in range(self.__dupa(j)) ]
              # IR intensities
              line = file.readline()
              irints[(j*5):j*5+self.__dupa(j)] =\
              [ numpy.float64(line.replace('D','E').split()[-self.__dupa(j):][x])\
                                              for x in range(self.__dupa(j)) ]
              # Eigenvectors
              line = file.readline()
              line = file.readline()
              for i in range(self.__n3):
                  eigvec[i][(j*5):j*5+self.__dupa(j)] =\
                  [ numpy.float64(line.replace('D','E').split()[-self.__dupa(j):][x])\
                                                  for x in range(self.__dupa(j)) ]
                  if (i+1)==self.__n3:
                     for h in range(3): line = file.readline()
                  else: line = file.readline()
              #for h in range(4): line = file.readline()
              
          freqs  = freqs [::-1]
          redmss = redmss[::-1]
          forcec = forcec[::-1]
          irints = irints[::-1]
          eigvec = eigvec[:,::-1]
          modelP = modelP[::-1]
          realP  = realP [::-1]
          
       # save
       self.__atoms = atoms
       self.__pos = pos
       if freq: 
          self.__freq = freqs
          self.__redmss = redmss
          self.__forcec = forcec
          self.__irints = irints
          self.__eigvec = eigvec
          if oniom:
             self.__modelP = modelP
             self.__realP  = realP
          else:
             self.__modelP = None
             self.__realP  = None
       return
   
   # U T I L I T I E S
   def __dupa(self,j):
       """some strange but extremely helpful utility:D"""
       if self.__nmodes-j*5 >= 5: return 5
       else                   : return self.__nmodes%5

   

if import_matplotlib:
   class QMOscillator(Graph):
       """Solver for 1D QM oscillator confined to a given potential.
   
   Usage:
   
   wfn = QMOscillator()
   wfn.set(self,V,x0,xn,np=500,y1=0.01,match_center=False,**kwargs)
   e,x,y = wfn.eval(<your guess for eigenvalue>)
   e,x,y = wfn.get()
   wfn.plot(name='fig.eps')
   
   Notes:
   
   Inputs to set():
   -----------------------------------------------------------------------
   V   - function object for the potential 
   kwargs - optional arguments for V
   x0  - starting point
   xn  - ending point
   y1  - arbitrary second small value needed to start Numerov integrations
   np  - number of intervals between x0 and xn
   match_center - if True - the matching will be done in the center 
                  of the argument range; else: in the xr point being
                  the argument for V(xr) = guess
   
   Returned by eval():
   -----------------------------------------------------------------------
   e   - eigenvalue for the solution
   x,y - argument and wave-function values for these arguments
   """
       
       def __init__(self,**kwargs):
           if kwargs: self.set(**kwargs)
           return
   
       # public
   
       def set(self,V,x0,xn,np=500,y1=0.01,match_center=False,**kwargs):
           "set the potential function and other parameters"
           self.__V = V
           self.__Vxr = lambda xr,guess: self.__V(xr,**kwargs) - guess
           self.__x0 = numpy.float64(x0)
           self.__xn = numpy.float64(xn)
           self.__y0 = 0.000
           self.__y1 = numpy.float64(y1)
           self.__np = np
           self.__step = (xn-x0)/np
           self.__x = numpy.linspace(x0,xn,np+1)
           self.__v = V(self.__x,**kwargs)
           self.__match_center = match_center
           self.__fig1 = None
           return
   
       def eval(self,guess):
           "calculate eigenvalue and eigenfunction"
           xr = self._find_xr(guess)
           self.__eps = secant(self._eigen,guess,xr=xr)
           # normalize the wavefunction
           prob = self.__wfn*self.__wfn
           Nsq  = simpson(prob,self.__step)
           self.__wfn /= math.sqrt(Nsq)
           return self.__eps, self.__x, self.__wfn
       
       def get(self):
           "return the final quantities"
           return self.__eps, self.__x, self.__wfn
   
       def plot(self,name='fig.eps'):
           "plot the wavefunction"
           #if self.__fig1 is None:
           self.__fig1 = pylab.figure(figsize=(12,9))
           axes11 = self.__fig1.add_subplot(111)
           axes12 = axes11.twinx()
           
           axes11.set_xlabel('$x$',fontsize=self.fs)
           axes11.set_ylabel('$V(x)$',fontsize=self.fs)
           axes12.set_ylabel('$\\psi(x)$',fontsize=self.fs)
           
           axes11.plot(self.__x,self.__v,'-',color='gray',
                       label='$V(x)$',lw=1)
           axes11.legend(loc=1,numpoints=1,
                         prop={'size':20}).draw_frame(False)
                         
           axes12.plot(self.__x,self.__wfn,'-',color='green',
                       label='$\\epsilon=%4.5f$'%self.__eps,lw=2)
           axes12.legend(loc=2,numpoints=1,
                         prop={'size':20}).draw_frame(False)
                         
           pylab.plt.xticks(fontsize=10, fontproperties=self.font1)
           pylab.plt.yticks(fontsize=10, fontproperties=self.font1)
           pylab.plt.draw()
           pylab.plt.show()
           pylab.plt.savefig(name)
           return
       
       # protected
   
       def _eigen(self,guess,xr):
           "calculate eigenvalue"
           # find matching point
           if self.__match_center:
              ixl = self.__np/2+6
           else:
               for i in range(self.__np):
                  if self.__x[i] > xr:
                     break
               ixl = i
           ixr = self.__np-ixl
           # left pre-solution
           q  = 2.0 * (guess - self.__v[:ixl+2])
           x,yl = numerov2(q,self.__x0,self.__y0,self.__y1,
                           ixl+2,self.__step)
   
           f2 = yl[ixl]
           # right pre-solution 
           q  = 2.0 * (guess - self.__v[ixl-1:])
           q  = q[::-1]
           x,yr = numerov2(q,self.__xn,self.__y0,self.__y1,
                           ixr+2,self.__step)
           k1 = yr[ixr-1] # r+h
           k2 = yr[ixr  ] # r
           k3 = yr[ixr+1] # r-h
           # rescale the left pre-solution
           yl *= k2/f2
           f1 = yl[ixl+1]   # l+h
           f2 = yl[ixl  ]   # l
           f3 = yl[ixl-1]   # l-h
           #
           yr = yr[::-1]
           f = ((f1-f3)-(k1-k3))/(2.0*self.__step*k2)
           # store the wave functions
           self.__yl = yl
           self.__yr = yr
           self.__wfn = numpy.concatenate([yl[:-2],yr[1:]])
           return f
   
       def _find_xr(self,guess):
           "find the matching point"
           xr = secant(self.__Vxr, 1.0 , guess=guess)
           return xr
    

def check_sim(l):
    """check the sim list"""
    log = ' --- OK ---'
    for x,y in l:
        i=0;j=0
        for a,b in l:
            if a==x: i+=1
            if b==y: j+=1
        if (i>1 or j>1): 
            log = " --- !ERROR! --- "
            break
    return log
        
def lind(file,querry):
    """
Parses line indices containing the querry
Usage:

ind = lind(file,querry)

Notes: returns list of integers. You can
go the the line of this file by:
for i in range(index_of_a_line): file_obj.readline()
line = file_obj.readline() is the desired line!
"""
    ind = []
    #if file.closed: data = open(file) # dopisz to potem zeby nie czytał pliku jescze raz od  nowa!!!!!
    data = open(file)
    line = data.readline()
    n = 0
    while line:
       if querry in line:
          ind.append(n)
       n+=1
       line = data.readline()
    data.close()
    return ind

class VIB(units.UNITS):
    """
Represents vibrational analysis tool
    """
    def __init__(self,mol,hess,weight=True,proj_grad=None):
        self.hess = hess
        self.mol = mol
        self._prepare()
        if weight: self._weight(proj_grad)

    # public
    
    def eval(self):
        """find normal modes doing tANDr-projection-out"""
        # transformation matrix
        self.__u = self._projout()
        F = numpy.dot(numpy.transpose(self.__u),numpy.dot(self.hess,self.__u))
        #F = dot(self.__u,dot(self.hess,transpose(self.__u)))
        E = numpy.linalg.eigh(F)[0]
        # frequencies
        self.__freq = numpy.where(E>0.,
            numpy.sqrt( E)*self.HartreePerHbarToCmRec,
           -numpy.sqrt(-E)*self.HartreePerHbarToCmRec)
        # reduced masses and cartesian l-matrix
        self.__u_cart, self.__redmass = self._redmass()
        return
    
    def get(self):
        """get frequencies, reduced masses and normal modes"""
        return self.__freq, self.__redmass, self.__u

    def __repr__(self):
        """print happy me!"""
        log = '\n'
        return str(log)
        
    # protected
    
    def _M(self):
        """M matrix"""
        N = len(self.mol.atoms)
        M = numpy.zeros((N*3,N*3),dtype=numpy.float64)
        for i in xrange(N):
            m = 1./numpy.sqrt(self.masses[i])
            M[3*i+0,3*i+0] = m
            M[3*i+1,3*i+1] = m
            M[3*i+2,3*i+2] = m
        return M
    def _M1(self):
        """M^-1 matrix"""
        N = len(self.mol.atoms)
        M = numpy.zeros((N*3,N*3),dtype=numpy.float64)
        for i in xrange(N):
            m = numpy.sqrt(self.masses[i])
            M[3*i+0,3*i+0] = m
            M[3*i+1,3*i+1] = m
            M[3*i+2,3*i+2] = m
        return M        

    def _redmass(self):
        """calculate reduced masses and cartesian l matrix"""
        u_cart = numpy.dot(self._M(),self.__u)
        redmass = 1./numpy.sum(u_cart**2,axis=0)*self.ElectronMassToAmu
        # normalize u_cart
        u_cart = u_cart/numpy.sqrt(numpy.sum(u_cart**2,axis=0))
        return u_cart, redmass
    
    def _prepare(self):
        """prepare coordinates and masses"""
        coords = []; masses = []
        for atom in self.mol.atoms:
            coords.append(atom.pos())
            masses.append(atom.mass())
        self.masses = numpy.array(masses,dtype=numpy.float64)*self.AmuToElectronMass
        self.coords = numpy.array(coords,dtype=numpy.float64)
        return
    def _weight(self,proj_grad=None):
        """weight hessian"""
        M = self._M()
        self.hess = numpy.dot(M,numpy.dot(self.hess,M))
        if proj_grad is not None:
           g = numpy.dot(M,proj_grad)
           g/=-numpy.linalg.norm(g)
           T = numpy.identity(len(g)) - numpy.outer(g,g)
           self.hess = numpy.dot(numpy.dot(T,self.hess),T)
        return
    def _weight_old(self):
        """weight Hessian"""
        for i in xrange(len(self.mol.atoms)):
            mi = 1./numpy.sqrt(self.masses[i])
            for j in xrange(len(self.mol.atoms)):
                mj = 1./numpy.sqrt(self.masses[j])
                mij = mi*mj
                self.hess[3*i+0,3*j+0] *= mij
                self.hess[3*i+1,3*j+1] *= mij
                self.hess[3*i+2,3*j+2] *= mij
                
                self.hess[3*i+0,3*j+1] *= mij
                self.hess[3*i+1,3*j+0] *= mij
                self.hess[3*i+0,3*j+2] *= mij
                self.hess[3*i+2,3*j+0] *= mij
                self.hess[3*i+1,3*j+2] *= mij
                self.hess[3*i+2,3*j+1] *= mij
        return
    def _rcom(self):
        """calculate center of mass"""
        return numpy.sum(self.coords*self.masses[:,numpy.newaxis],axis=0)/sum(self.masses)
    def _r(self):
        """translate all coordinates to R_COM"""
        return self.coords - self._rcom()
    def _I(self):
        """moment of inertia tensor"""
        I = numpy.zeros((3,3),numpy.float64)
        r = self._r()
        for i in xrange(len(self.coords)):
            x,y,z = r[i]
            m = self.masses[i]
            I[0,0]+= m*(y**2+z**2)
            I[1,1]+= m*(x**2+z**2)
            I[2,2]+= m*(x**2+y**2)
            I[0,1]-= m*x*y
            I[0,2]-= m*x*z
            I[1,2]-= m*y*z
        I[1,0] = I[0,1]
        I[2,0] = I[0,2]
        I[2,1] = I[1,2]
        return I
    def _vec_t(self):
        """translational vectors"""
        D1 = []; D2 = []; D3 = []
        for i in xrange(len(self.coords)):
            m = numpy.sqrt(self.masses[i])
            v1 = [m,0,0]
            v2 = [0,m,0]
            v3 = [0,0,m]
            D1+=v1; D2+=v2; D3+=v3
        D1 = self._norm(numpy.array(D1,numpy.float64))
        D2 = self._norm(numpy.array(D2,numpy.float64))
        D3 = self._norm(numpy.array(D3,numpy.float64))
        return D1,D2,D3
    def _vec_r(self):
        """rotational vectors"""
        D4 = []; D5 = []; D6 = []
        X = numpy.linalg.eigh(self._I())[1]
        X1,X2,X3 = X
        r = self._r()
        for i in xrange(len(self.coords)):
            m = numpy.sqrt(self.masses[i])
            ri = r[i]
            P1 = numpy.dot(ri,X1)
            P2 = numpy.dot(ri,X2)
            P3 = numpy.dot(ri,X3)
            d4 = [(P2*X[0,2]-P3*X[0,1])/m,
                  (P2*X[1,2]-P3*X[1,1])/m,
                  (P2*X[2,2]-P3*X[2,1])/m]
            d5 = [(P3*X[0,0]-P1*X[0,2])/m,
                  (P3*X[1,0]-P1*X[1,2])/m,
                  (P3*X[2,0]-P1*X[2,2])/m]
            d6 = [(P1*X[0,1]-P2*X[0,0])/m,
                  (P1*X[1,1]-P2*X[1,0])/m,
                  (P1*X[2,1]-P2*X[2,0])/m]
            #P = array([P1,P2,P3],float64)
            #d1 = (cross(P,X1)/m).tolist()
            #d2 = (cross(P,X2)/m).tolist()
            #d3 = (cross(P,X3)/m).tolist()
            D4+=d4; D5+=d5; D6+=d6
        D4 = self._norm(numpy.array(D4,numpy.float64))
        D5 = self._norm(numpy.array(D5,numpy.float64))
        D6 = self._norm(numpy.array(D6,numpy.float64))
        return D4,D5,D6
    def _norm(self,u):
        """normalize vector u"""
        return u/numpy.sqrt(numpy.dot(u,u))
    def _gs(self,u,v):
        """Gram-Schmidt ortogonalization of vector v wrt u"""
        puv = numpy.dot(v,u)/numpy.dot(u,u)*u
        return v-puv
    def _pre_diag(self):
        """diagonalize Hessian matrix"""
        return numpy.linalg.eigh(self.hess)[1]
    def _projout(self):
        """project out translations and rotations"""
        U = self._pre_diag()
        u = U[:,6:]
        #u = U[6:]
        #t1,t2,t3,r1,r2,r3 = transpose(U[:,:6])
        d1,d2,d3 = self._vec_t()
        d4,d5,d6 = self._vec_r()
        for i in xrange(u.shape[1]):
            vec = u[:,i]
            #vec = u[i]
            for trvec in [d1,d2,d3,d4,d5,d6]:
                vec = self._gs(trvec,vec)
            #u[:,i] = self._norm(vec)
            #u[i] = self._norm(vec)
        return U
    
class Peak2DIR(units.UNITS):
    """
Represent a 2DIR peak in the spectrum of signals

Usage:

a = Peak(x=<coordsx>,y=<coordsy>,z=<coordsz>,Tw=<waiting time>,
         t_max=1.0     ,dt=0.008,
         n_points=2**10,w_cent=1500.0,
         wx_max=2000.0 ,wx_min=0.0,
         wy_max=2000.0 ,wy_min=0.0)
a.set_peak(<n>=1,<func_name>='g')
a.fit(<parameters>,[method='slsqp',bounds=[],**kwargs])
a.get_r2()
a.get_parameters()
print peak

Notes:

 1) <x> and <y> are 1D lists or numpy ndarrays
 2) <z> is 2D ndarray or a list of 2D ndarrays
    of length N = len(Tw) (corresponding to each
    waiting time Tw)
 3) <Tw> is single waiting time or the list (ndarray)
    of waiting times of length N
 4) <n> is number of subpeaks subject to extract
 5) 'g' - Gaussian distributions (default)
 6) <parameters> is the list of lists or tuples
    of the form:

      [['par1',val1],
       ['par2',val2],
           . . .
       ['parn',valn]]

    Parameter (par) have to be strings 
    whereas initial values (valn) floats.

    <parameters> depend on a type of function
    selected. Here I list {parn} names:
    a) Gaussian 1 peak:
       A, sigma, x_o
    b) Gaussian 2 peaks:
       A_1, sigma_1, x_o1,
       A_2, sigma_2, x_o2
       and so on
"""
    def __init__(self, x, y, z, Tw,
                 t_max=1.0     , dt=0.008     ,
                 n_points=2**10, w_cent=1500.0,
                 wx_max=2000.0 , wx_min=0.0   ,
                 wy_max=2000.0 , wy_min=0.0   ,):
        
        self.ToCmRec = 1.e-12*self.SpeedOfLight*100.0
        self.FromAngToCmRec = self.ToCmRec*2.*math.pi
        dt *= 2.*math.pi
        self.__nTw = None
        
        self.__t_max_1D    = 60.0
        self.__n_points_1D = 2**10
        self.__time_1D     = numpy.linspace(0.,self.__t_max_1D,self.__n_points_1D) * 2.*math.pi
        self.__dt_1D       = self.__time_1D[1]-self.__time_1D[0]
        
        ### range of calculation
        self.freq  = numpy.fft.fftshift( numpy.fft.fftfreq(n_points,
                                               d=dt*self.ToCmRec) ) + w_cent
        self.__kx = numpy.where(numpy.logical_and(self.freq<wx_max, self.freq>wx_min))[0]
        self.__ky = numpy.where(numpy.logical_and(self.freq<wy_max, self.freq>wy_min))[0]
        self.X = self.freq.copy()[self.__kx]
        self.Y = self.freq.copy()[self.__ky]
        
        
        ### MULTI-WAITING-TIME MODE
        try:
            self.__nTw = len(Tw)
            
            ### simulation grid
            self.sim_grid = Grid3D(xmin=0.0,xmax=t_max*2.*math.pi,
                                   ymin=0.0,ymax=t_max*2.*math.pi,
                                   zmin=0  ,zmax=self.__nTw-1,
                                   dx=dt,dy=dt,dz=1.)
            self.sim_grid.newaxis(2,Tw)
            
            ### experimental data
            Z = numpy.zeros((self.Y.shape[0],self.X.shape[0],self.__nTw),numpy.float64)
            for i in xrange(self.__nTw):
                exp_grid = scipy.interpolate.RectBivariateSpline(y,x,z[:,:,i])
                Z[:,:,i] = exp_grid(self.Y,self.X)
            self.Z = Z
            self.Z/= self.Z.max()
            
        ### SINGLE-WAITING-TIME MODE
        except TypeError:
            
            ### simulation grid
            self.exp_grid = scipy.interpolate.RectBivariateSpline(y,x,z)
            self.sim_grid = Grid2D(xmin=0.0,xmax=t_max*2.*math.pi,
                                   ymin=0.0,ymax=t_max*2.*math.pi,
                                   dx=dt,dy=dt)
        
            ### experimental data
            self.Z = self.exp_grid(self.Y,self.X)
            self.Z/= self.Z.max()
        
        ### memorize important data
        self.x = x
        self.y = y
        self.z = z
        self.func = None
        self.init = None
        self.out  = None
        self.status =-1
        self.param  = None
        self.__func = None
        self.n      = None
        self.__n_points = n_points
        self.__dt = dt
        self.Tw = Tw
        self.w_cent = w_cent
        
    # public
    def set_peak(self,n=1,func_name='r'):
        """set the type of peak"""
        self.__func = func_name
        
        if self.__nTw is None:
           if func_name=='r':
              if   n==1: self.func = self._r1_no_exch
              elif n==2: self.func = self._r2_no_exch
        else:
           if func_name=='r':
              if   n==1: self.func = self._r1_no_exch_m
              elif n==2: self.func = self._r2_no_exch_m
        self.n = n

    def set_1D(self,t_max,n_points):
        """
Settings for 1D FTIR simulation

Usage:

t_max        - maximum time in [ps] 
n_points     - number of time points (should be a power of 2).
               Good value is 2**10
"""
        self.__t_max_1D    = t_max
        self.__n_points_1D = n_points
        self.__time_1D     = numpy.linspace(0.,self.__t_max_1D,self.__n_points_1D) * 2.*math.pi
        self.__dt_1D       = self.__time_1D[1]-self.__time_1D[0]
        return
    
    def set_ftir(self,freq,ftir,w_01_1D,w_min,w_max,thresh=0.990):
        """
Set the benchmark FTIR data and acceptance threshold. 

Usage:

freq         - <x> frequencies in [cm-1]
ftir         - <y> experimental FTIR spectrum [dimensionless]
w_min, w_max - spectral window to compare with simulation
thresh       - R² threshold for simulation fit to experimental 
               FTIR spectrum. Default is 0.990
w_01_1D      - experimental center frequency (now it is useless)
"""
        ### set threshold
        assert 1.0>thresh>0.4, "Threshold has to be between 0.4 and 1.0"
        self.__ftir_threshold = thresh
        ### frequency range
        freq_1d        = numpy.fft.fftshift( numpy.fft.fftfreq(self.__n_points_1D,
                                             d=self.__dt_1D*self.ToCmRec) ) + self.w_cent
        self.__k_1D    = numpy.where(numpy.logical_and(freq_1d<w_max, freq_1d>w_min))[0]
        self.__freq_1D = freq_1d.copy()[self.__k_1D]
        self.__w_01_1D = w_01_1D
        ### interpolate the data
        f = scipy.interpolate.interp1d(freq,ftir)
        self.__ftir = f(self.__freq_1D)
        return
        
    def get_parameters(self):
        """get the parameters of the fitting"""
        #print " No fitting performed - no parameters."
        return self.param
    
    def get_fit(self):
        """return the fitted data"""
        if self.param is None: 
           #print " No fitting performed - no parameters."
           return None
        else:
           return self.func(**self.args)

    def get_peaks(self):
        """return all separate peaks"""
        peaks = []
        for i in xrange(self.n):
            if self.__func == 'r':
               w_01   = self.param[8*i+0]
               anh    = self.param[8*i+1]
               delta_1= self.param[8*i+2]
               delta_2= self.param[8*i+3]
               tau_1  = self.param[8*i+4]
               tau_2  = self.param[8*i+5]
               T1     = self.param[8*i+6]
               T2     = self.param[8*i+7]
               peak   = self._r1_no_exch(w_01,anh,delta_1,delta_2,tau_1,tau_2,T1,T2)

            peaks.append(peak)
        return numpy.array(peaks,dtype=numpy.float64)

    def get_r2(self):
        """return R^2 coefficient of fitting"""
        data_av = numpy.average(self.Z)
        sse = numpy.sum((self.func(**self.args)-self.Z)**2)
        sst = numpy.sum((self.Z-data_av)**2)
        return 1.0 - sse/sst
    
    def get_ftir(self,normalize=False):
        """
Return FTIR spectrum from simulated parameters.
        
Options:

normalize (default: False)
"""
        args = self.args.copy()
        args.update({'normalize': normalize, 'w_01':self.__w_01_1D})
        freq, ftir = self._spectrum_1D(**args)
        return freq, ftir
    
    def get_ftir_exp(self):
        """return experimental FTIR (after interpolation)"""
        _args = self.args.copy()
        _args.update({'normalize':True,'w_01':self.__w_01_1D})
        return self.__freq_1D, self.__ftir, self._spectrum_1D(**_args)[1][self.__k_1D]
    
    ### METHODS
    
    def fit(self,opts,method='slsqp',bounds=[],ieqcons=[],
            epsilon=1e-06,acc=1e-6,max_iter=1000,disp=2,
            ftircons=False,full_output=True):
        """perform fitting using [options] list"""
        self.__set_param(opts)

        if method=='slsqp':
           if ftircons: ieqcons.append( self._ftir_constr )
           
           result  = scipy.optimize.fmin_slsqp(self._residsq,
                                          self.init_list,iter=max_iter,
                                          acc=acc,disp=disp,bounds=bounds,
                                          args=(opts,),full_output=full_output,
                                          epsilon=epsilon,ieqcons=ieqcons)
           param ,a,b,c,s = result
           self.status = c
        #
        self.param = param
        self.__update_args(param,opts)
    
    
    # protected
    def _resid(self,params,opts):
        """residual function for optimization"""
        self.__update_args(params,opts)
        return self.Z - self.func(**self.args)

    def _residsq(self,params,opts):
        """square residual function for optimization"""
        self.__update_args(params,opts)
        return numpy.sum((self.Z - self.func(**self.args))**2.)
    
    ### Response Functions
    def __response_1D(self, t, w_01, mu_01, Delta, tau, T1, T2):
        """First-order response function"""
        w_off = (w_01-self.w_cent) * self.FromAngToCmRec
        R = 0.
        Tor = 1000000.0
        def g(t):
            G=0.
            for i in range(2):
                G+= tau[i]*Delta[i] * \
          ( numpy.exp(-t/tau[i])*tau[i] - tau[i] + t ) 
            G += t/T2
            return G
        r = mu_01 * mu_01 * numpy.exp(-g(t) -t/(Tor*3.) -t/(T1*2.) -1.j*w_off*t)
        #r*= pop_ratio[n]
        R+= r       
        return R 
    
    def __response_3D(self,t3,t1,
                   Tw, w_01, mu_01, mu_12, anh, Delta, tau, T1, T2):
        """3-rd order response function for 3-level system (in time domain)"""
        w_off = -(w_01-self.w_cent) * self.FromAngToCmRec
        anh  *= self.FromAngToCmRec
        # di-Kubo line-shape function
        def g(t):
            G=0.
            for i in range(2):
                G+= tau[i]*Delta[i] * \
          ( numpy.exp(-t/tau[i])*tau[i] - tau[i] + t ) 
            G += t/T2 + t/(2.*T1) + t/(3.*10000000.0)
            return G
        #
        M = numpy.exp(-Tw/T1) * (1.0 + 0.8 * numpy.exp(-Tw/10000000.0))
        RR = 0.0; NR = 0.0
        
        mu_01_2= mu_01*mu_01
        mu_12_2= mu_12*mu_12
        
        gt1     = g(t1)
        gTw     = g(Tw)
        gt3     = g(t3)
        gt1Tw   = g(t1+Tw)
        gTwt3   = g(Tw+t3)
        gt1Twt3 = g(t1+Tw+t3)
        #print gt1Twt3
        #try: print gt1Twt3[:,:,0]
        #except: pass
        #M =1
        # --- Rephasing
        # R1 and R2
        RR+= mu_01_2*mu_01_2*numpy.exp(-1.j*w_off*(-t1+t3))*M*2.0*\
             numpy.exp(-gt1+gTw-gt3-gt1Tw-gTwt3+gt1Twt3)
        # R3
        RR-= mu_01_2*mu_12_2*numpy.exp(-1.j*(w_off*(-t1+t3)-anh*t3))*M*\
             numpy.exp(-gt1+gTw-gt3-gt1Tw-gTwt3+gt1Twt3)
        # --- Nonrephasing
        # R4 and R5
        NR+= mu_01_2*mu_01_2*numpy.exp(-1.j*w_off*(t1+t3))*M*2.0*\
             numpy.exp(-gt1-gTw-gt3+gt1Tw+gTwt3-gt1Twt3)
        # R6
        NR-= mu_01_2*mu_12_2*numpy.exp(-1.j*(w_off*(t1+t3)-anh*t3))*M*\
             numpy.exp(-gt1-gTw-gt3+gt1Tw+gTwt3-gt1Twt3)
        #
        #try: print real(RR[:,:,0])
        #except: pass
        #print real(RR)
        return RR, NR
    
    def _ftir_constr(self, *args):
        """constraint for FTIR spectrum"""
        data_av = numpy.average(self.__ftir)
        _args = self.args.copy()
        _args.update({'normalize':True,'w_01':self.__w_01_1D})
        sim_dat = self._spectrum_1D(**_args)[1][self.__k_1D]
        sse = numpy.sum((sim_dat-self.__ftir)**2)
        sst = numpy.sum((self.__ftir-data_av)**2)
        r2 = 1.0 - sse/sst
        return r2 - self.__ftir_threshold
    
    def _spectrum_1D(self, w_01, anh, delta_1, delta_2, tau_1, tau_2, T1, T2,
                     mu_01=1., normalize=False):
        #             t_max=None, n_points=None):
        """1D FTIR spectrum. Returns tuple (x [cm-1], y [intensity])"""
        #if t_max is None:
        #    t_max    = self.__t_max_1D
        #    n_points = self.__n_points_1D

        #time_1D = linspace(0.,t_max,n_points) * 2.*math.pi
        #dt = time_1D[1]-time_1D[0]
        spectrum = self.__response_1D( self.__time_1D, w_01, mu_01, 
                                      numpy.array([delta_1,delta_2]), 
                                      numpy.array([tau_1,tau_2]), T1, T2)
        ftir = numpy.fft.fftshift(numpy.fft.fft(spectrum))
        freq = numpy.fft.fftshift(numpy.fft.fftfreq(self.__n_points_1D,d=self.__dt_1D*self.ToCmRec)) + self.w_cent
        data_f = real(ftir)[::-1]
        data_f-= data_f.min()
        if normalize: data_f/=data_f.max()
        return freq, data_f
    
    def _r1_no_exch(self,w_01,anh,delta_1,delta_2,tau_1,tau_2,T1,T2):
        """3-rd order response without exchange and coupling"""
        
        ### signal in time-domain
        rr, nr = self.sim_grid.eval(self.__response_3D, 
                                    # assumed parameters
                                    Tw=self.Tw, mu_01=1., mu_12=math.sqrt(2.),
                                    # optimizing parameters
                                    w_01=w_01, 
                                    anh=anh, 
                                    Delta=numpy.array([delta_1,delta_2]), 
                                    tau=numpy.array([tau_1,tau_2]), 
                                    T1=T1, T2=T2, )
        #print rr
        ### rephasing and non-rephasing spectras (2D FFT)
        data_rr_f = numpy.fft.fftshift( numpy.fft.fft2(rr,s=(self.__n_points,self.__n_points)) )
        data_nr_f = numpy.fft.fftshift( numpy.fft.fft2(nr,s=(self.__n_points,self.__n_points)) )
        data_rr_f = data_rr_f[:,::-1]
        data_rr_f = numpy.roll(data_rr_f,1,axis=1)
        
        ### total signal
        data_f = numpy.real(data_rr_f + data_nr_f)
        
        data_f = data_f[self.__kx,:]
        data_f = data_f[:,self.__ky]
        data_f/= data_f.max()
        data_f = data_f.transpose()
        
        return data_f
    
    def _r1_no_exch_m(self,w_01,anh,delta_1,delta_2,tau_1,tau_2,T1,T2):
        """3-rd order response without exchange and coupling"""
        
        ### signal in time-domain
        #self.__response_3D_f = vectorize(self.__response_3D)
        rr, nr = self.sim_grid.eval(self.__response_3D, 
                                    # assumed parameters
                                    mu_01=1., mu_12=math.sqrt(2.),
                                    # optimizing parameters
                                    w_01=w_01, 
                                    anh=anh, 
                                    Delta=numpy.array([delta_1,delta_2]), 
                                    tau=numpy.array([tau_1,tau_2]), 
                                    T1=T1, T2=T2, )
        #print rr[:,:,0]
        ### rephasing and non-rephasing spectras (2D FFT)
        #data_rr_f = numpy.zeros((self.__n_points,self.__n_points,self.__nTw),numpy.complex64)
        #data_nr_f = data_rr_f.copy()
        
        #for i in range(self.__nTw):
        #    data_rr_f[:,:,i] = numpy.fft.fftshift( numpy.fft.fft2(rr[:,:,i],s=(self.__n_points,self.__n_points)) )
        #    data_nr_f[:,:,i] = numpy.fft.fftshift( numpy.fft.fft2(nr[:,:,i],s=(self.__n_points,self.__n_points)) )
        
        data_rr_f = numpy.fft.fftshift( numpy.fft.fft2(rr,axes=(0,1),s=(self.__n_points,self.__n_points)), axes=(0,1))
        data_nr_f = numpy.fft.fftshift( numpy.fft.fft2(nr,axes=(0,1),s=(self.__n_points,self.__n_points)), axes=(0,1))
        data_rr_f = data_rr_f[:,::-1,:]
        #for i in range(3):
        #      data_rr_f[:,:,i] = data_rr_f[:,::-1,i]
        #    data_rr_f[:,:,i] = numpy.roll(data_rr_f[:,:,i],1,axis=1)
        data_rr_f = numpy.roll(data_rr_f,1,axis=1)
        
        ### total signal
        data_f = numpy.real(data_rr_f + data_nr_f)
        
        data_f = data_f[self.__kx,:,:]
        data_f = data_f[:,self.__ky,:]
        data_f/= data_f[:,:,0].max()
        #for i in range(1):
        #    data_f[:,:,i]/=data_f[:,:,i].max()
        #max_vals = numpy.amax(amax(data_f,1),0)
        #data_f/= max_vals
        data_f = numpy.transpose(data_f,(1,0,2))
        #for i in range(self.__nTw):
        #    print data_f[:,:,i].max()
        #print self.sim_grid.zcoorv
        return data_f
        
    # private
    def __set_param(self,param):
        args, init_list = self.__make_args(param)
        self.opts = param
        self.args = args
        self.init_list = init_list
        ### change the unit of frequencies
        

    def __make_args(self,init):
        args = dict(init)
        init_list=[]
        for opt in init:
            init_list.append(opt[1])
        return dict(init), init_list

    def __update_args(self,params,opts):
        for i in range(len(opts)):
            self.args[opts[i][0]] = params[i]

    def __repr__(self):
        """print the data report"""
        log = "\n"
        log+= " ====================================================\n"
        log+= "                  SIGNAL DESCRIPTION\n"
        log+= " ====================================================\n"
        d_name = {'r'  :'3-rd order Response without exchange'    ,
                  'l'  :'3-rd order Response with exchange'       ,
                  None :'Not assigned'          }
        log+= '\n'
        log+= " - Profile: %s\n" % d_name[self.__func]
        # peak generall overview
        if self.n is not None:
           log+= " - Peak No: %i\n" % self.n
        else:
           log+= " - Peak No: %s\n" % 'Not assigned'
        log+= '\n'
        # peak parameters and analysis
        if self.param is not None:
           log+= "   PARAMETERS\n"
           log+= "\n"
           p = self.param
           # ... function types ...
           if self.__func == 'r':
              l1 = " %6s"%('Peak'.rjust(6))
              l2 = " %6s"%(('%s'%letters.greek.omega+'_01').rjust(6))
              l3 = " %6s"%(letters.greek.Delta.rjust(6))
              l4 = " %6s"%('Peak'.rjust(6))
              l5 = " %6s"%('Peak'.rjust(6))
              l6 = " %6s"%('Peak'.rjust(6))
              l7 = " %6s"%('Peak'.rjust(6))
              l8 = " %6s"%('Peak'.rjust(6))
              l9 = " %6s"%('Peak'.rjust(6))
           for line in [l1,l2,l3,l4,l5,l6,l7,l8,l9]:
               log+= line + '\n'
           log+= "\n"
           log+= " ----------------------------------------------------\n"
           log+= (" R² = %10.6f" % self.get_r2()).center(52)
           log+= "\n"
        else:
           log+= " ----------------------------------------------------\n"
           log+= " No fitting was performed".center(52)
           log+= "\n"
        log+= " ----------------------------------------------------\n"            
        return str(log)
    
        
class Peak:
    """
Represent a peak in the spectrum of signals

Usage:

a = Peak(x=<coordsx>,y=<coordsy>)
a.set_peak(<n>=1,<func_name>='g')
a.fit(<parameters>,[method='slsqp',bounds=[],**kwargs])
a.get_r2()
a.get_parameters()
a.get_fwhm()
print peak

Notes:

 1) <x> and <y> are lists or numpy arrays
 2) <n> is number of subpeaks subject to extract
 3) 'g' - Gaussian distributions (default)
 4) <parameters> is the list of lists or tuples
    of the form:

      [['par1',val1],
       ['par2',val2],
           . . .
       ['parn',valn]]

    Parameter (par) have to be strings 
    whereas initial values (valn) floats.

    <parameters> depend on a type of function
    selected. Here I list {parn} names:
    a) Gaussian 1 peak:
       A, sigma, x_o
    b) Gaussian 2 peaks:
       A_1, sigma_1, x_o1,
       A_2, sigma_2, x_o2
       and so on
"""
    def __init__(self,x,y):
        self.x = x
        self.y = y
        self.func = None
        self.init = None
        self.out  = None
        self.status = -1
        self.param = None
        self.__func = None
        self.n     = None

    # public
    def set_peak(self,n=1,func_name='g'):
        """set the type of peak"""
        self.__func = func_name
        
        if func_name=='g':
           if   n==1: self.func = self._gauss1
           elif n==2: self.func = self._gauss2
           elif n==3: self.func = self._gauss3
           elif n==4: self.func = self._gauss4
           
        elif func_name=='l':
           if   n==1: self.func = self._lorentz1
           elif n==2: self.func = self._lorentz2
           elif n==3: self.func = self._lorentz3
           elif n==4: self.func = self._lorentz4
           
        elif func_name=='lg1':
           if   n==1: self.func = self._lg11
           elif n==2: self.func = self._lg12
           elif n==3: self.func = self._lg13
           elif n==4: self.func = self._lg14

        elif func_name=='lg2':
           if   n==1: self.func = self._lg21
           elif n==2: self.func = self._lg22
           elif n==3: self.func = self._lg23
           elif n==4: self.func = self._lg24

        elif func_name=='v':
           if   n==1: self.func = self._voigt1
           elif n==2: self.func = self._voigt2
           
        self.n = n
    
    def get_parameters(self):
        #print " No fitting performed - no parameters."
        return self.param
    
    def get_fit(self):
        if self.param is None: 
           #print " No fitting performed - no parameters."
           return None
        else:
           return self.func(**self.args)

    def get_fwhm(self):
        """calculate FWHM for peaks"""
        fwhm = []
        for i in xrange(self.n):
            if (self.__func == 'g' or self.__func == 'l'):
               fwhm.append(self.param[3*i+1])
            elif self.__func == 'lg1':
               fwhm.append(self.param[4*i+1])
            elif self.__func == 'lg2':
               mL    = self.param[5*i+4]
               mG    = 1. - mL
               sigmaL= self.param[5*i+1]
               sigmaG= self.param[5*i+2]
               value = mL * sigmaL + mG * sigmaG
               fwhm.append(value)
            elif self.__func == 'v':
               sigmaL= self.param[4*i+1]
               sigmaG= self.param[4*i+2]
               value = 0.5346 * sigmaL + math.sqrt(0.2166 * sigmaL**2. + sigmaG**2.)
               fwhm.append(value)
        return numpy.array(fwhm,numpy.float64)
               
    def get_peaks(self):
        peaks = []
        for i in xrange(self.n):
            if self.__func == 'g':
               x_0   = self.param[3*i+0]
               sigma = self.param[3*i+1]
               A     = self.param[3*i+2]
               peak  = self._gauss1  (x_0,sigma,A)
            elif self.__func == 'l':
               x_0   = self.param[3*i+0]
               sigma = self.param[3*i+1]
               A     = self.param[3*i+2]
               peak  = self._lorentz1(x_0,sigma,A)
            elif self.__func == 'lg1':
               x_0   = self.param[4*i+0]
               sigma = self.param[4*i+1]
               A     = self.param[4*i+2]
               m     = self.param[4*i+3]
               peak  = self._lg11(x_0,sigma,A,m)
            elif self.__func == 'lg2':
               x_0   = self.param[5*i+0]
               sigmaL= self.param[5*i+1]
               sigmaG= self.param[5*i+2]
               A     = self.param[5*i+3]
               m     = self.param[5*i+4]
               peak  = self._lg21(x_0,sigmaL,sigmaG,A,m)
            elif self.__func == 'v':
               x_0   = self.param[4*i+0]
               sigmaL= self.param[4*i+1]
               sigmaG= self.param[4*i+2]
               A     = self.param[4*i+3]
               peak  = self._voigt1(x_0,sigmaL,sigmaG,A)
            peaks.append(peak)
        return numpy.array(peaks,dtype=numpy.float64)

    def get_r2(self):
        """return R^2 coefficient of fitting"""
        data_av = numpy.average(self.y)
        sse = sum((self.func(**self.args)-self.y)**2)
        sst = sum((self.y-data_av)**2)
        return 1 - sse/sst
        
    def fit(self,opts,method='leastsq',disp=1,bounds=[],
            epsilon=1e-08,pgtol=1e-012,factr=100.0,m=8000,
            approx_grad=True,fprime=None,maxfun=2000000000,acc=1e-13):
        """perform fitting using [options] list"""
        self.__set_param(opts)
        if method=='leastsq':
           param, flag = scipy.optimize.leastsq(self._resid,
                                          self.init_list,
                                          args=(opts,),
                                          ftol=1e-12,
                                          xtol=1e-12,maxfev=2000000000,)
           self.status= flag
        elif method=='l-bfgs-b':
           if bounds  ==[]: bounds=None
           param, f, d = scipy.optimize.fmin_l_bfgs_b(self._residsq,
                                          self.init_list,factr=factr,
                                          args=(opts,),pgtol=pgtol,
                                          bounds=bounds,fprime=fprime,
                                          disp=disp,m=m,approx_grad=approx_grad,
                                          epsilon=epsilon,maxfun=maxfun,)
        elif method=='slsqp':
           result  = scipy.optimize.fmin_slsqp(self._residsq,
                                          self.init_list,iter=maxfun,
                                          acc=acc,disp=2,bounds=bounds,
                                          args=(opts,),full_output=True)
           param ,a,b,c,s = result
           self.status = c
        #
        self.param = param
        self.__update_args(param,opts)

    # protected 
    def _resid(self,params,opts):
        """residual function for optimization"""
        self.__update_args(params,opts)
        return self.y - self.func(**self.args)

    def _residsq(self,params,opts):
        """square residual function for optimization"""
        self.__update_args(params,opts)
        return numpy.sum((self.y - self.func(**self.args))**2.)
    
    ### Normal distribution

    def _normal(self,xo_1,sigma_1,A_1):
        """single Gaussian distribution"""
        return (A_1/(sigma_1*math.sqrt(2*math.pi)))\
               * numpy.exp(-(self.x-xo_1)**2/(2*sigma_1**2))
               
    ### pure Gaussian profiles
                           
    def _gauss1(self,xo_1,sigma_1,A_1):
        """single Gaussian distribution"""
        return A_1/sigma_1 * math.sqrt(4.*math.log(2.)/math.pi) * \
               numpy.exp(-4.*math.log(2.)/(sigma_1**2.)*(self.x-xo_1)**2.)

    def _gauss2(self,xo_1,sigma_1,A_1,
                     xo_2,sigma_2,A_2):
        """bimodal gaussian distribution"""
        g1 = A_1/sigma_1 * math.sqrt(4.*math.log(2.)/math.pi) * \
               numpy.exp(-4.*math.log(2.)/(sigma_1**2.)*(self.x-xo_1)**2.)
        g2 = A_2/sigma_2 * math.sqrt(4.*math.log(2.)/math.pi) * \
               numpy.exp(-4.*math.log(2.)/(sigma_2**2.)*(self.x-xo_2)**2.)
        return A_1*g1 + A_2*g2

    def _gauss3(self,xo_1,sigma_1,A_1,
                     xo_2,sigma_2,A_2,
                     xo_3,sigma_3,A_3):
        """trimodal gaussian distribution"""
        g1 = A_1/sigma_1 * math.sqrt(4.*math.log(2.)/math.pi) * \
               numpy.exp(-4.*math.log(2.)/(sigma_1**2.)*(self.x-xo_1)**2.)
        g2 = A_2/sigma_2 * math.sqrt(4.*math.log(2.)/math.pi) * \
               numpy.exp(-4.*math.log(2.)/(sigma_2**2.)*(self.x-xo_2)**2.)
        g3 = A_3/sigma_3 * math.sqrt(4.*math.log(2.)/math.pi) * \
               numpy.exp(-4.*math.log(2.)/(sigma_3**2.)*(self.x-xo_3)**2.)
        return A_1*g1 + A_2*g2 + A_3*g3

    def _gauss4(self,xo_1,sigma_1,A_1,
                     xo_2,sigma_2,A_2,
                     xo_3,sigma_3,A_3,
                     xo_4,sigma_4,A_4):
        """trimodal gaussian distribution"""
        g1 = A_1/sigma_1 * math.sqrt(4.*math.log(2.)/math.pi) * \
               numpy.exp(-4.*math.log(2.)/(sigma_1**2.)*(self.x-xo_1)**2.)
        g2 = A_2/sigma_2 * math.sqrt(4.*math.log(2.)/math.pi) * \
               numpy.exp(-4.*math.log(2.)/(sigma_2**2.)*(self.x-xo_2)**2.)
        g3 = A_3/sigma_3 * math.sqrt(4.*math.log(2.)/math.pi) * \
               numpy.exp(-4.*math.log(2.)/(sigma_3**2.)*(self.x-xo_3)**2.)
        g4 = A_4/sigma_4 * math.sqrt(4.*math.log(2.)/math.pi) * \
               numpy.exp(-4.*math.log(2.)/(sigma_4**2.)*(self.x-xo_4)**2.)

        return A_1*g1 + A_2*g2 + A_3*g3 + A_4*g4

    ### pure Lorentzian profiles
    
    def _lorentz1(self,xo_1,sigma_1,A_1):
        """single Lorenzian distribution"""
        return 2.*A_1/math.pi * sigma_1/(4.*(self.x-xo_1)**2.+sigma_1**2.)

    def _lorentz2(self,xo_1,sigma_1,A_1,
                       xo_2,sigma_2,A_2):
        """bimodal Lorenzian distribution"""
        l1 = 2.*A_1/math.pi * sigma_1/(4.*(self.x-xo_1)**2.+sigma_1**2.)
        l2 = 2.*A_2/math.pi * sigma_2/(4.*(self.x-xo_2)**2.+sigma_2**2.)
        return l1 + l2

    def _lorentz3(self,xo_1,sigma_1,A_1,
                       xo_2,sigma_2,A_2,
                       xo_3,sigma_3,A_3):
        """trimodal Lorenzian distribution"""
        l1 = 2.*A_1/math.pi * sigma_1/(4.*(self.x-xo_1)**2.+sigma_1**2.)
        l2 = 2.*A_2/math.pi * sigma_2/(4.*(self.x-xo_2)**2.+sigma_2**2.)
        l3 = 2.*A_3/math.pi * sigma_3/(4.*(self.x-xo_3)**2.+sigma_3**2.)
        return l1 + l2 + l3
    
    def _lorentz4(self,xo_1,sigma_1,A_1,
                       xo_2,sigma_2,A_2,
                       xo_3,sigma_3,A_3,
                       xo_4,sigma_4,A_4):
        """4-modal Lorenzian distribution"""
        l1 = 2.*A_1/math.pi * sigma_1/(4.*(self.x-xo_1)**2.+sigma_1**2.)
        l2 = 2.*A_2/math.pi * sigma_2/(4.*(self.x-xo_2)**2.+sigma_2**2.)
        l3 = 2.*A_3/math.pi * sigma_3/(4.*(self.x-xo_3)**2.+sigma_3**2.)
        l4 = 2.*A_4/math.pi * sigma_4/(4.*(self.x-xo_4)**2.+sigma_4**2.)
        return l1 + l2 + l3 + l4
    
    ### pseudo-Voigt-1 profiles
    
    def _lg11(self,xo_1,sigma_1,A_1,m_1):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigma_1/(4.*(self.x-xo_1)**2. + sigma_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigma_1**2.)*(self.x-xo_1)**2.)/sigma_1
        return A_1 * lg1
    
    def _lg12(self,xo_1,sigma_1,A_1,m_1,
                   xo_2,sigma_2,A_2,m_2):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigma_1/(4.*(self.x-xo_1)**2. + sigma_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigma_1**2.)*(self.x-xo_1)**2.)/sigma_1
              
        lg2 = m_2*2./math.pi * sigma_2/(4.*(self.x-xo_2)**2. + sigma_2**2.) + \
              (1.-m_2)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigma_2**2.)*(self.x-xo_2)**2.)/sigma_2
        return A_1 * lg1 + A_2 * lg2

    def _lg13(self,xo_1,sigma_1,A_1,m_1,
                   xo_2,sigma_2,A_2,m_2,
                   xo_3,sigma_3,A_3,m_3,):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigma_1/(4.*(self.x-xo_1)**2. + sigma_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigma_1**2.)*(self.x-xo_1)**2.)/sigma_1
              
        lg2 = m_2*2./math.pi * sigma_2/(4.*(self.x-xo_2)**2. + sigma_2**2.) + \
              (1.-m_2)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigma_2**2.)*(self.x-xo_2)**2.)/sigma_2

        lg3 = m_3*2./math.pi * sigma_3/(4.*(self.x-xo_3)**2. + sigma_3**2.) + \
              (1.-m_3)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigma_3**2.)*(self.x-xo_3)**2.)/sigma_3
        return A_1 * lg1 + A_2 * lg2 + A_3 * lg3
    
    def _lg14(self,xo_1,sigma_1,A_1,m_1,
                   xo_2,sigma_2,A_2,m_2,
                   xo_3,sigma_3,A_3,m_3,
                   xo_4,sigma_4,A_4,m_4,):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigma_1/(4.*(self.x-xo_1)**2. + sigma_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigma_1**2.)*(self.x-xo_1)**2.)/sigma_1
              
        lg2 = m_2*2./math.pi * sigma_2/(4.*(self.x-xo_2)**2. + sigma_2**2.) + \
              (1.-m_2)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigma_2**2.)*(self.x-xo_2)**2.)/sigma_2

        lg3 = m_3*2./math.pi * sigma_3/(4.*(self.x-xo_3)**2. + sigma_3**2.) + \
              (1.-m_3)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigma_3**2.)*(self.x-xo_3)**2.)/sigma_3
              
        lg4 = m_4*2./math.pi * sigma_4/(4.*(self.x-xo_4)**2. + sigma_4**2.) + \
              (1.-m_4)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigma_4**2.)*(self.x-xo_4)**2.)/sigma_4
        return A_1 * lg1 + A_2 * lg2 + A_3 * lg3 + A_4 * lg4

    ### pseudo-Voigt-2 profiles
    
    def _lg21(self,xo_1,sigmaL_1,sigmaG_1,A_1,m_1):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigmaL_1/(4.*(self.x-xo_1)**2. + sigmaL_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigmaG_1**2.)*(self.x-xo_1)**2.)/sigmaG_1
        return A_1 * lg1

    def _lg22(self,xo_1,sigmaL_1,sigmaG_1,A_1,m_1,
                   xo_2,sigmaL_2,sigmaG_2,A_2,m_2):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigmaL_1/(4.*(self.x-xo_1)**2. + sigmaL_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigmaG_1**2.)*(self.x-xo_1)**2.)/sigmaG_1
              
        lg2 = m_2*2./math.pi * sigmaL_2/(4.*(self.x-xo_2)**2. + sigmaL_2**2.) + \
              (1.-m_2)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigmaG_2**2.)*(self.x-xo_2)**2.)/sigmaG_2
        return A_1 * lg1 + A_2 * lg2

    def _lg23(self,xo_1,sigmaL_1,sigmaG_1,A_1,m_1,
                   xo_2,sigmaL_2,sigmaG_2,A_2,m_2,
                   xo_3,sigmaL_3,sigmaG_3,A_3,m_3):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigmaL_1/(4.*(self.x-xo_1)**2. + sigmaL_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigmaG_1**2.)*(self.x-xo_1)**2.)/sigmaG_1
              
        lg2 = m_2*2./math.pi * sigmaL_2/(4.*(self.x-xo_2)**2. + sigmaL_2**2.) + \
              (1.-m_2)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigmaG_2**2.)*(self.x-xo_2)**2.)/sigmaG_2

        lg3 = m_3*2./math.pi * sigmaL_3/(4.*(self.x-xo_3)**2. + sigmaL_3**2.) + \
              (1.-m_3)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigmaG_3**2.)*(self.x-xo_3)**2.)/sigmaG_3
        return A_1 * lg1 + A_2 * lg2 + A_3 * lg3
    
    def _lg24(self,xo_1,sigmaL_1,sigmaG_1,A_1,m_1,
                   xo_2,sigmaL_2,sigmaG_2,A_2,m_2,
                   xo_3,sigmaL_3,sigmaG_3,A_3,m_3,
                   xo_4,sigmaL_4,sigmaG_4,A_4,m_4):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigmaL_1/(4.*(self.x-xo_1)**2. + sigmaL_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigmaG_1**2.)*(self.x-xo_1)**2.)/sigmaG_1
              
        lg2 = m_2*2./math.pi * sigmaL_2/(4.*(self.x-xo_2)**2. + sigmaL_2**2.) + \
              (1.-m_2)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigmaG_2**2.)*(self.x-xo_2)**2.)/sigmaG_2

        lg3 = m_3*2./math.pi * sigmaL_3/(4.*(self.x-xo_3)**2. + sigmaL_3**2.) + \
              (1.-m_3)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigmaG_3**2.)*(self.x-xo_3)**2.)/sigmaG_3
              
        lg4 = m_4*2./math.pi * sigmaL_4/(4.*(self.x-xo_4)**2. + sigmaL_4**2.) + \
              (1.-m_4)*math.sqrt(4.*math.log(2.)/math.pi)*numpy.exp( (-4.*math.log(2.)/sigmaG_4**2.)*(self.x-xo_4)**2.)/sigmaG_4
        return A_1 * lg1 + A_2 * lg2 + A_3 * lg3 + A_4 * lg4

    ### pure Voigt profiles
    
    def __v(self,t,x,xc,wl,wg,a):
        A = math.exp(-t**2.)
        B = (math.sqrt(math.log(2.))*wl/wg)**2.
        C = (math.sqrt(4.*math.log(2.)) * (x-xc)/wg - t)**2.
        return a*2.*math.log(2.)/math.pi**(3./2.) * wl/wg**2 * A/(B+C)

    def _voigt1(self,xo_1,sigmaL_1,sigmaG_1,A_1):
        y = numpy.zeros(len(self.x),dtype=numpy.float64)
        for i in xrange(len(self.x)):
            y[i] = scipy.integrate.quad(self.__v,-inf,inf,full_output=0,args=(self.x[i],xo_1,sigmaL_1,sigmaG_1,A_1))[0]
        return y
    
    def _voigt2(self,xo_1,sigmaL_1,sigmaG_1,A_1,
                     xo_2,sigmaL_2,sigmaG_2,A_2):
        y = numpy.zeros(len(self.x),dtype=numpy.float64)
        for i in xrange(len(self.x)):
            val = scipy.integrate.quad(self.__v,-inf,inf,full_output=0,args=(self.x[i],xo_1,sigmaL_1,sigmaG_1,A_1))[0]
            val+= scipy.integrate.quad(self.__v,-inf,inf,full_output=0,args=(self.x[i],xo_2,sigmaL_2,sigmaG_2,A_2))[0]
            y[i] = val
        return y

                                                  
    # private
    def __set_param(self,param):
        args, init_list = self.__make_args(param)
        self.opts = param
        self.args = args
        self.init_list = init_list

    def __make_args(self,init):
        args = dict(init)
        init_list=[]
        for opt in init:
            init_list.append(opt[1])
        return dict(init), init_list

    def __update_args(self,params,opts):
        for i in range(len(opts)):
            self.args[opts[i][0]] = params[i]

    def __repr__(self):
        """print the data report"""
        log = "\n"
        log+= " ====================================================\n"
        log+= "                  SIGNAL DESCRIPTION\n"
        log+= " ====================================================\n"
        d_name = {'g'  :'Pure Gaussian'         ,
                  'l'  :'Pure Lorentzian'       ,
                  'lg1':'Pseudo-Voigt 1 profile',
                  'lg2':'Pseudo-Voigt 2 profile',
                  'v'  :'Exact Voigt profile'   ,
                  None :'Not assigned'          }
        log+= '\n'
        log+= " - Profile: %s\n" % d_name[self.__func]
        # peak generall overview
        if self.n is not None:
           log+= " - Peak No: %i\n" % self.n
        else:
           log+= " - Peak No: %s\n" % 'Not assigned'
        log+= '\n'
        # peak parameters and analysis
        if self.param is not None:
           log+= "   PARAMETERS\n"
           log+= "\n"
           p = self.param
           f = self.get_fwhm()
           if (self.__func == 'g' or self.__func == 'l'):
               log+= " %6s %8s %8s %8s %8s\n"%('Peak'.rjust(6),
                                               'Freq'.rjust(8),
                                               'sigm'.rjust(8),
                                               'Area'.rjust(8),
                                               'FWHM'.rjust(8))
               for i in range(self.n):
                   log+= " %6i %8.2f %8.2f %8.4f %8.2f\n" % (i+1,p[i*3+0],
                                                             p[i*3+1],p[i*3+2],
                                                             f[i])
           elif self.__func == 'lg1':
               log+= " %6s %8s %8s %8s %8s %8s\n"%('Peak'.rjust(6),
                                                   'Freq'.rjust(8),
                                                   'sigm'.rjust(8),
                                                   'Area'.rjust(8),
                                                   'mixL'.rjust(8),
                                                   'FWHM'.rjust(8))
               for i in range(self.n):
                   log+= " %6i %8.2f %8.2f %8.4f %8.4f %8.2f\n" % (i+1,p[i*4+0],
                                                                   p[i*4+1],p[i*4+2],
                                                                   p[i*4+3],f[i])

           elif self.__func == 'lg2':
               log+= " %6s %8s %8s %8s %8s %8s %8s\n"%('Peak'.rjust(6),
                                                       'Freq'.rjust(8),
                                                       'sigL'.rjust(8),
                                                       'sigG'.rjust(8),
                                                       'Area'.rjust(8),
                                                       'mixL'.rjust(8),
                                                       'FWHM'.rjust(8))
               for i in range(self.n):
                   log+= " %6i %8.2f %8.2f %8.2f %8.4f %8.4f %8.2f\n" % (i+1,p[i*5+0],
                                                                         p[i*5+1],p[i*5+2],
                                                                         p[i*5+3],p[i*5+4],
                                                                         f[i])

           elif self.__func == 'v':
               log+= " %6s %8s %8s %8s %8s %8s\n"%('Peak'.rjust(6),
                                                   'Freq'.rjust(8),
                                                   'sigL'.rjust(8),
                                                   'sigG'.rjust(8),
                                                   'Area'.rjust(8),
                                                   'FWHM'.rjust(8))
               for i in range(self.n):
                   log+= " %6i %8.2f %8.2f %8.2f %8.4f %8.2f\n" % (i+1,p[i*4+0],
                                                                   p[i*4+1],p[i*4+2],
                                                                   p[i*4+3],f[i])
           log+= "\n"
           log+= " ----------------------------------------------------\n"
           log+= (" R² = %10.6f" % self.get_r2()).center(52)
           log+= "\n"
        else:
           log+= " ----------------------------------------------------\n"
           log+= " No fitting was performed".center(52)
           log+= "\n"
        log+= " ----------------------------------------------------\n"
        return str(log)
            
def interchange(T,ind):
    """\
interchange rows according to order list

Usage: ingerchange(array,order_list)

Returns: permuted array

Example: 

Assume we have initial matrix T. We want obtain T':

T =  [[ 1  6  8]         T'=  [[ 3 -4 -4]
      [ 2  5  7]               [ 7 -1  9]
      [ 3 -4 -4]      ind      [ 1  6  8]
      [ 4  1  0]     ---->     [ 4  1  0]
      [ 5 -7 -8]               [ 2  5  7]
      [ 6  6  2]               [ 5 -7 -8]
      [ 7 -1  9]]              [ 6  6  2]]
      
This is accomplished by run:
T_prime = interchange(T,ind=[3, 7, 1, 4, 2, 5, 6])
"""
    ind = numpy.array(ind)-1
    B = T.copy()
    for i,index in enumerate(ind):
        B[i] = T[index]
    return B

def choose(a,ids):
    """\
choose the array barray

Usage:

choose(barray,index_list)
it returns the barray similar to array,
but without elements of array placed in
indices given by index_list (indices are 
Python-like!)"""
    inranges = len(ids)-1
    t = []
    # first range
    t.append(a[:ids[0]])
    # inranges
    for r in xrange(inranges):
        t.append( a[ ids[r]+1: ids[r+1] ] )
    # last range
    t.append(a[ids[-1]+1:])
    #
    return numpy.concatenate(t)

def get_pmloca(natoms,mapi,sao,vecin,nae,
               maxit=1000,conv=1.0E-06,lprint=False,
               freeze=None):
    """\
===============================================================
Performs the Pipek-Mezey molecular orbital localization.
Reference: 
J. PIPEK AND P. G. MEZEY  J. CHEM. PHYS. 90, 4916 (1989)

Usage:
get_pmloca(natoms,mapi,sao,vecin,nae,
           [maxit=1000,conv=1.0E-06,lprint=False,
            freeze=None])

returns:
transormation matrix: 2-d ndarray of shape(nmos,nmos)
transformed orbitals: 2-d ndarray of shape(nmos,nbasis)

NMOS and NBASIS are taken from the dimensions of VECIN.
If you want to exclude some MOs provide MOs indices in
freeze list (Python convenction,N-1). The program will
get a slice ov VECIN and return full transformation U
with frozen orbitals too (not written yet).
---------------------------------------------------------------
arguments:
natoms - number of atoms in molecule
mapi   - list of atoms in basis set order (LIST1 in PyQuanteM)
sao    - array of AO overlap integrals of size (nbasis,nbasis)
vecin  - input MO coefficients
nae    - number of alpha electrons

optional:
maxit  - maximum number of iterations
conv   - convergence for electron population
lprint - whether print no of iteration or not after finish
===============================================================
"""
    # freeze the orbitals requested
    if freeze is not None: 
       vecin = choose(vecin,freeze)
    # mapi in FORTRAN default indexes, nmos and no of elements
    # in triangular matrix (n2) for PM localizator matrix elements
    mapi = numpy.array(mapi,int) + 1
    nmos = len(vecin)
    n2   = (nmos+1)*nmos/2
    #
    tran = qm.pmloca.pmloca(natoms=natoms,mapi=mapi,sao=sao,vecin=vecin,
                         maxit=maxit,cvgloc=conv,n2=n2,nae=nae,
                         lprint=lprint)
    #
    tran = numpy.transpose(tran)
    vecout = numpy.dot(tran,vecin)
    #
    return tran, vecout

def reorder(P,sim,axis=0):
    """Reorders the tensor according to <axis> (default is 0). 
<sim> is the list of pairs from 'order' function. 
In normal numbers (starting from 1...)"""
    P_new = numpy.zeros(P.shape,dtype=numpy.float64)
    if   axis==0:
         for i,j in sim:
             P_new[i-1] = P[j-1]
    elif axis==1:
         for i,j in sim:
             P_new[:,i-1] = P[:,j-1]
    elif axis==2:
         for i,j in sim:
             P_new[:,:,i-1] = P[:,:,j-1]
    return P_new

def order(R,P,start=0,lprint=1):
    """order list"""
    new_P = P.copy()
    sim   = []
    rad =  []
    for i in range(len(R)-start):
        J = 0+start
        r = 1.0E+100
        rads = []
        for j in range(len(P)-start):
            r_ = numpy.sum(( R[i+start]-P[j+start])**2)
            r__= numpy.sum((-R[i+start]-P[j+start])**2)
            if r__<r_: r_=r__
            rads.append(r_)
            if r_<r:
               r=r_
               J = j
        sim.append((i+1,J+1))
        new_P[i+start] = P[J+start]
        rad.append(rads)
    for i in xrange(len(R)-start):
        s = numpy.sum(numpy.sign(new_P[i])/numpy.sign(R[i]))
        if lprint: print "%10d %f" %(i+1,s)
        r_ = sum(( R[i+start]-new_P[i+start])**2)
        r__= sum((-R[i+start]-new_P[i+start])**2)
       
        #if s < -154: 
        #   print "TUTAJ s < -154"
        #   #new_P[i]*=-1.
        if r__<r_:
          if lprint: print "    HERE r__ < r_ (sign reversal)"
          new_P[i]*=-1.
    return new_P, sim#, array(rad,dtype=float)

class GROUPS:
      """ 
 Grouping algorithm from numerical project:
 ------------------------------------------
      Assignemt for the purpose of:        
 ==========================================
 ADVANCED PROGRAMMING AND NUMERICAL METHODS
 ==========================================
 Teacher    : dr inż. Paweł Szarek         
 Author     : inż.    Bartosz Błasiak      
 Affiliation:                              
     - - - - - - - - - - - - - - - - -     
     Wrocław University of Technology      
    - - - - - - - - - - - - - - - - -      
                 CUBEFILER (c) 2012        
"""
      def __init__(self,A):
          self.groups = self.__make_GROUPS(A)
          
      def __make_GROUPS(self,A):
          """group items according to boolean relation matrix A_ij"""
          def add(a,b):
              c = []
              for i in b: c.append(i)
              for i in a:
                  if i not in b: c.append(i)
              return c

          def add_line(i,group,A):
              for j in range(len(A)):
                  if A[i][j] and j not in group: group.append(j)
              if i not in group: group.append(i)
              return group

          def make_group(i,A):
              g_old = add_line(i,[],A)
              g_new1 = []
              g_new2 = []
              while len(g_new1)!=len(g_old):
                     for j in g_old:
                         g_new1 = add_line(j,g_old,A)
                         g_new1 = add(g_new1,g_old)
                     for j in g_new1:
                         g_new2 = add_line(j,g_new1,A)
                         g_new2 = add(g_new2,g_new1)
                     g_old = g_new1
                     g_new = g_new2
              return g_new2

          GROUPS = []
          G_0 = make_group(0,A)
          GROUPS.append(G_0)
          for i in range(len(A)-1):
              t=i+1
              Q=[]
              for G in GROUPS:
                  Q+=G
              if t not in Q:
                 GROUPS.append(make_group(t,A))
          return GROUPS
         

def get_tcf(file,nmax,norgns,ndels,
                 nskip=0,lprint=False,
                 save=False,outfile='tcf.out'):
    """\
Computes time auto-correlation function.

Usage:
get_tcf(file,nmax,norgns,ndels,
            [nkip=0,lprint=False,
             save=False,outfile='tcf.out'])

returns:
2-d ndarray of shape(ndels,2)

arguments:
file   - input file with observable in two-column
         FORTRAN format (I6,13.5D) 
nmax   - the number of time steps in the file
norgns - number of time origins for tcf
         norgns < nmax
ndels  - maximum delay time in time steps
         for which tcf is to be computed
nskip  - number of time steps to be omitted
         from trajectory in input file
lprint - whether print the results or not
save   - whether save the tcf file or not
outfile- if save: provide the name of output
         tcf file 
"""
    # compute tcf from data input file
    r = qm.gentcf.gentcf(file,nmax,nskip,norgns,lprint,ndels)
    r = r[:ndels]
    # write the tcf file on disk
    if save:
       out = open(outfile,'w')
       for i in range(ndels):
           print >> out, "%10i %13.5E" % ((i+1),r[i])
       out.close()
    # return results:
    tcf = numpy.zeros((ndels,2),dtype=numpy.float64)
    tcf[:,0] = numpy.linspace(1,ndels,ndels)
    tcf[:,1] = numpy.array(r)

    return tcf

def DistanceRelationMatrix(xyz,threshold=1):
    """calculate boolean relation matrix for structure xyz (array).
    Threshold = 1 Bohr and coordinates of xyz have to be Bohr too!
    You can then search for groups using: GROUPS(A_ij).groups"""
    
    K = len(xyz)
    A_ij = numpy.zeros((K,K),dtype=bool)
    for i in range(K):
        for j in range(i):
            if numpy.sqrt(numpy.sum((xyz[i]-xyz[j])**2))<=threshold:
               A_ij[i,j] = True
               A_ij[j,i] = True
               
    return A_ij

### SVDSuperimposer from BIOPYTHON PACKAGE
# Copyright (C) 2002, Thomas Hamelryck (thamelry@vub.ac.be)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# [!] zmiana w oryginalnym kodzie:
#        w run() zmieniono na:     
#           av1=sum(coords,axis=0)/self.n  
#           av2=sum(reference_coords,axis=0)/self.n 
#        z powodu złej wersji metody sum()
class SVDSuperimposer(object):
    """
    SVDSuperimposer finds the best rotation and translation to put
    two point sets on top of each other (minimizing the RMSD). This is 
    eg. useful to superimpose crystal structures.  

    SVD stands for Singular Value Decomposition, which is used to calculate
    the superposition.

    Reference:

    Matrix computations, 2nd ed. Golub, G. & Van Loan, CF., The Johns 
    Hopkins University Press, Baltimore, 1989
    """
    def __init__(self):
        self._clear()

    # Private methods

    def _clear(self):
        self.reference_coords=None
        self.coords=None
        self.transformed_coords=None
        self.rot=None
        self.tran=None
        self.rms=None
        self.init_rms=None

    def _rms(self, coords1, coords2):
        "Return rms deviations between coords1 and coords2."
        diff=coords1-coords2
        l=coords1.shape[0]
        return numpy.sqrt(numpy.sum(numpy.sum(diff*diff))/l)

    # Public methods
    
    def set(self, reference_coords, coords):
        """
        Set the coordinates to be superimposed.
        coords will be put on top of reference_coords.

        o reference_coords: an NxDIM array
        o coords: an NxDIM array

        DIM is the dimension of the points, N is the number
        of points to be superimposed.
        """
        # clear everything from previous runs
        self._clear()
        # store cordinates
        self.reference_coords=reference_coords
        self.coords=coords
        n=reference_coords.shape
        m=coords.shape
        if n!=m or not(n[1]==m[1]==3):
            raise Exception("Coordinate number/dimension mismatch.")
        self.n=n[0]

    def run(self):
        "Superimpose the coordinate sets."
        if self.coords is None or self.reference_coords is None:
            raise Exception("No coordinates set.")
        coords=self.coords
        reference_coords=self.reference_coords
        # center on centroid
        av1=numpy.sum(coords,axis=0)/self.n  
        av2=numpy.sum(reference_coords,axis=0)/self.n    
        coords=coords-av1
        reference_coords=reference_coords-av2
        # correlation matrix
        a=numpy.dot(numpy.transpose(coords), reference_coords)
        u, d, vt=numpy.linalg.svd(a)
        self.rot=numpy.transpose(numpy.dot(numpy.transpose(vt), numpy.transpose(u)))
        # check if we have found a reflection
        if numpy.linalg.det(self.rot)<0:
            vt[2]=-vt[2]
            self.rot=numpy.transpose(numpy.dot(numpy.transpose(vt), numpy.transpose(u)))
        self.tran=av2-numpy.dot(av1, self.rot)

    def get_transformed(self):
        "Get the transformed coordinate set."
        if self.coords is None or self.reference_coords is None:
            raise Exception("No coordinates set.")
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        if self.transformed_coords is None:
            self.transformed_coords=numpy.dot(self.coords, self.rot)+self.tran
        return self.transformed_coords

    def get_rotran(self):
        "Right multiplying rotation matrix and translation."
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        return self.rot, self.tran

    def get_init_rms(self):
        "Root mean square deviation of untransformed coordinates."
        if self.coords is None:
            raise Exception("No coordinates set yet.")
        if self.init_rms is None:
            self.init_rms=self._rms(self.coords, self.reference_coords)
        return self.init_rms

    def get_rms(self):
        "Root mean square deviation of superimposed coordinates."
        if self.rms is None:
            transformed_coords=self.get_transformed()
            self.rms=self._rms(transformed_coords, self.reference_coords)
        return self.rms


def MakeMol(atno,coord,name='dummy',**opts):
    """Make Molecule object from a a list of atomic numbers and coordinates.
The default units are Bohrs.
Usage: mol = MakeMol(atno,coord,name='dummy',**opts)
Notes: atno is a list or array of atomic numbers, 
       coord is an array of size (natoms,3) in Bohrs,
       unless units='Angstroms' was also specified."""
    coords = []
    for i in xrange(len(atno)):
        atom  = (atno[i], tuple(coord[i]))
        coords.append(atom)
    Mol = PyQuante.Molecule(name,coords,**opts)
    return Mol

def Read_xyz_file(file,ar=False,mol=False,mult=1,charge=0,name='dummy',
                  units='Angstrom',method='RHF',basis='6-311++G**'):
    """\
reads xyz or fchk file and returns coords and atoms. 
Coordinates are returned in AU!

Usage:
get_mol(file,[ar=False,mol=False,mult=1,charge=0,name='dummy',
              units='Angstrom',method='RHF',basis='6-311++G**'])

Returns: 
+ list of coords plus ...
+ numpy.array of xyz (if ar=True) 
+ or Molecule object

Description:
file   - xyz file
units  - input units (output is always in Bohr!)
mol    - wheather return Molecule object
name   - name of Molecule object
mult   - multiplicity
charge - charge
ar     - return also array with only coordinates
"""
    ### handle *.xyz files
    if file[-4:].lower() == '.xyz':
       data = open(file).readlines()
       n_atoms = int(data[0])
       data.pop(0);data.pop(0)
       coord = []
       for i in range(n_atoms):
           coord.append(data[i].split()[:4])
           coord[i][1:] = map(numpy.float64,coord[i][1:])
           if units.lower()=='angstrom':
              for j in range(3):
                  coord[i][j+1]*= uUNITS.AngstromToBohr
            
       if ar:
           data = [map(numpy.float64,[x for x in coord[y][1:]]) \
                                  for y in range( len(coord))]
           data = numpy.array(data,dtype=numpy.float64)
       if mol:
           Coords = []
           for i in range(n_atoms):
               atom  = (uAtom(coord[i][0]).atno, (coord[i][1], 
                                                  coord[i][2],
                                                  coord[i][3]) )
               Coords.append(atom)
           Mol = PyQuante.Molecule(name,Coords,units='Bohr',
                                   multiplicity=mult,charge=charge,
                                   basis=basis,method=method)
    
       if   mol : return Mol                 
       elif ar  : return coord, data
       else     : return coord

    ### handle g09.fchk files
    elif file[-4:].lower() == 'fchk':
       print "WARNING: set charge and multiplicity correctly because they are not read from FCHK!"
       file = open(file)
       line = file.readline()
       g = lambda n,m: n/m+bool(n%m)
       
       # search for atomic numbers
       querry = "Atomic numbers"
       while True:
           if querry in line: break
           line = file.readline()       
       n_atoms = int(line.split()[-1])
       line = file.readline()
       
       atnos = []
       for i in xrange(g(n_atoms,6)):
           atnos += line.split()
           line = file.readline()
           
       atnos = numpy.array(atnos,dtype=int)
       
       # search for atomic coordinates       
       querry = "Current cartesian coordinates"
       while True:
           if querry in line: break
           line = file.readline()
       N = int(line.split()[-1])
       line = file.readline()
       
       coord = []
       for i in xrange(g(N,5)):
           coord += line.split()
           line = file.readline()
           
       coord = numpy.array(coord,dtype=numpy.float64).reshape(n_atoms,3)

       # create Molecule object
       if mol:
           Coords = []
           for i in range(n_atoms):
               atom  = (atnos[i], (coord[i][0],
                                   coord[i][1],
                                   coord[i][2]) )
               Coords.append(atom)
           Mol = PyQuante.Molecule(name,Coords,units='Bohr',
                                   multiplicity=mult,charge=charge,
                                   basis=basis,method=method)
    
           return Mol                        
       else: return None
       
def Vr_dma(dma,Rb,is_full=False):
    """calculates electrostatic potential in point Rb 
from dma distribution."""

    dmac = dma.copy()
    if not is_full:
       dmac.MAKE_FULL()
       dmac.MakeTraceless()
    Ra,qa,Da,Qa,Oa = dmac.DMA_FULL
    V=0
    for i in range(len(Ra)):
        R=Rb-Ra[i] 
        Rab=numpy.sqrt(numpy.sum(R**2,axis=0))
        V+=qa[i]/Rab
        V+=numpy.tensordot(Da[i],R,(0,0))/Rab**3 
        V+=numpy.tensordot(R,numpy.tensordot(Qa[i],R,(1,0)),(0,0))/Rab**5 
        V+=numpy.tensordot(R,numpy.tensordot(R,numpy.tensordot(Oa[i],R,(0,0)),(0,0)),(0,0))/Rab**7 
    return V

def Energy_density(dma,Rb,is_full=False):
    """calculates electrostatic potential in point Rb 
from dma distribution.
Energy_density(dma,R,full=False/True)
if full - the dma object is turned into traceless object
"""
 
    dmac = dma.copy()
    if not is_full:
       dmac.MAKE_FULL()
       dmac.MakeTraceless()
    Ra,qa,Da,Qa,Oa = dmac.DMA_FULL
    e2=0
    for i in range(len(Ra)):
        R=Rb-Ra[i]
        Rab=numpy.sqrt(numpy.sum(R**2,axis=0))
        e2+=qa[i]**2/Rab**4
        #
        e2+=3*(numpy.dot(Da[i],R))**2/Rab**8
        #
        e2+=(numpy.dot(Da[i],Da[i]))**2 / Rab**6
        #
        t  =numpy.tensordot(R,numpy.tensordot(Qa[i],R,(0,0)),(0,0))
        g  =numpy.tensordot(Qa[i],R,(0,0))
        e2+=5*t**2/Rab**12 +  4*numpy.dot(g,g)/Rab**10
        #
        t=numpy.tensordot(R,numpy.tensordot(R,numpy.tensordot(Oa[i],R,(0,0)),(0,0)),(0,0))
        e2+= 7*t**2 / Rab**16
        #
        t  =numpy.tensordot(R,numpy.tensordot(Oa[i],R,(0,0)),(0,0))
        g  =numpy.dot(t,t)
        e2+= 9*g / Rab**14
        
    return e2

def ElectricField(dma,Rb,is_full=False):
    """calculates electrostatic field in point Rb
from dma distribution. Usage:
ElectricField(dma,R,full=False/True)
if not is_full - the dma object is turned into traceless object
"""

    dmac = dma.copy()
    if not is_full:
       dmac.MAKE_FULL()
       dmac.MakeTraceless()
    Ra,qa,Da,Qa,Oa = dmac.DMA_FULL
    field=numpy.zeros(3,dtype=numpy.float64)
    for i in range(len(Ra)):
        R=Rb-Ra[i] 
        Rab=numpy.sqrt(numpy.sum(R**2,axis=0))
        field+= qa[i] * R / Rab**3
        
        field+= 3 * R * numpy.dot(R,Da[i]) / Rab**5
        field-= Da[i]/Rab**3
        
        t  =numpy.tensordot(R,numpy.tensordot(Qa[i],R,(0,0)),(0,0))
        field+= 5* t * R / Rab**7
        field-= 2* numpy.tensordot(Qa[i],R,(0,0)) / Rab**5
        
        c=numpy.tensordot(R,numpy.tensordot(Oa[i],R,(0,0)),(0,0))
        g=numpy.tensordot(R,c,(0,0))
        field+= 7 * g * R / Rab**9
        field-= 3 * c / Rab**7
        
    return field
    
def ParseDMA(file,type='coulomb',hexadecapoles=True):
    """\
============================================================================
Parse DMA from GAUSSIAN, GAMESS or COULOMB.py file. It returns a DMA object.
Usage:
ParseDMA(file, type='coulomb', hexadecapoles=False)
----------------------------------------------------------------------------
<type>s:
1) coulomb  or c (or just nothing - it is default)
2) gaussian or gau
3) gamess   or gms
Notes:
Gaussian file has to have pop=ChelpG printout in log
Gamess reads Stone's DMA analysis
============================================================================
"""
    if   type.lower() == 'slv':
         pass
         
    elif type.lower() == 'gamess' or type.lower() == 'gms':
         data = open(file)
         # ----------------------------------
         #querry1 = "PROPERTIES FOR THE B3LYP    DFT FUNCTIONAL (RHF  TYPE) DENSITY MATRIX"
         #querry2 = "MP2 PROPERTIES...FOR THE FIRST ORDER WAVEFUNCTION"
         line = data.readline()
         #while 1:
         #      if ((querry1 in line) or (querry2 in line)): break
         #      line = data.readline()
         # ----------------------------------
         querry = " NET CHARGES AT POINTS"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(3): line = data.readline()
         ZerothMoments = []
         Structure = []
         while line.split()!=[]:
               ZerothMoments.append( numpy.float64( line.split()[2] ) )
               Structure.append( numpy.array( line.split()[-3:],dtype=numpy.float64)  )
                    
               line = data.readline()
         # ----------------------------------
         querry = " FIRST MOMENTS AT POINTS"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(3): line = data.readline()
         FirstMoments = []
         while line.split()!=[]:
               FirstMoments.append( map( numpy.float64, line.split()[1:]))
               line = data.readline()
         # ----------------------------------
         querry = " SECOND MOMENTS AT POINTS"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(3): line = data.readline()
         SecondMoments = []
         while line.split()!=[]:
               SecondMoments.append( map( numpy.float64, line.split()[1:] ))
               line = data.readline()
         # ----------------------------------
         querry = " THIRD MOMENTS AT POINTS"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(4): line = data.readline()
         ThirdMoments = []
         while 'CPU' not in line.split():
               A = map( numpy.float64, line.split()[1:] )
               line = data.readline()
               B = map( numpy.float64, line.split() )
               ThirdMoments.append( A+B )
               line = data.readline()

         # add up the nuclear and electronic charges
         npoints = len(FirstMoments)
         natoms  = len(ZerothMoments) - npoints
         #if
         for i in range(natoms):
             ZerothMoments[i+natoms] += ZerothMoments[i]
         for i in range(natoms):
             ZerothMoments.pop(0)
             Structure.pop(0)

         data.close()
         return dma.DMA( q  =numpy.array(ZerothMoments)   ,
                         m  =numpy.array(FirstMoments )   ,
                         T  =numpy.array(SecondMoments)   , 
                         O  =numpy.array(ThirdMoments )   ,
                         pos=numpy.array(Structure)       ,
                         hexadecapoles=False              )

    # -----------------------------------------------------------------------------
    elif type.lower() == 'coulomb' or type.lower() == 'c':
         # determine if the file contains hexadecapole moments
         has_hexadecapoles = False
         data = open(file)
         if " Distributed fourth-order property" in data.read(): has_hexadecapoles = True
         data.close()

         # now read the file once again line by line
         data = open(file)
         line = data.readline()

         # start reading the sections
         querry = " Distributed zeroth-order property"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(3): line = data.readline()
         ZerothMoments = []
         Structure = []
         Origin = []
         while line.split()!=[]:
               ZerothMoments.append( numpy.float64( line.split()[0] ) )
               line = data.readline()
         nfrag = len(ZerothMoments)
         # ----------------------------------
         querry = " Distributed first-order property"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(4): line = data.readline()
         FirstMoments = []
         while line.split()!=[]:
               FirstMoments.append( map( numpy.float64, line.split()[:]))
               line = data.readline()
         # ----------------------------------
         querry = " Distributed second-order property"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(4): line = data.readline()
         SecondMoments = []
         while line.split()!=[]:
               SecondMoments.append( map( numpy.float64, line.split()[:] ))
               line = data.readline()
         # ----------------------------------
         querry = " Distributed third-order property"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(5): line = data.readline()
         ThirdMoments = []
         #while '-----' not in line or :
         for i in range(nfrag):
               A = map( numpy.float64, line.split()[:] )
               line = data.readline()
               B = map( numpy.float64, line.split()[:] )
               ThirdMoments.append( A+B )
               line = data.readline()
         # ----------------------------------
         FourthMoments = None
         if has_hexadecapoles and hexadecapoles:
            querry = " Distributed fourth-order property" 
            while 1:
                if querry in line: break
                line = data.readline()
            FourthMoments = []
            for i in range(6): line = data.readline() 
            for i in range(nfrag):
                A = map( numpy.float64, line.split()[:] )
                line = data.readline()
                B = map( numpy.float64, line.split()[:] )
                line = data.readline()
                C = map( numpy.float64, line.split()[:] )
                line = data.readline()
                FourthMoments.append( A+B+C )
            FourthMoments = numpy.array(FourthMoments)
         # ----------------------------------
         querry = " form@"
         while 1:
             if querry in line: break
             line = data.readline()
         typ = line.split()[0]
         if   typ == 'primitive': is_traceless = False
         elif typ == 'traceless': is_traceless = True
         # ------------------------------------
         querry = " Structure"
         struct = True
         while 1:
               if querry in line: break
               #print line
               line = data.readline()
               if not line: 
                  struct = False
                  break

         atoms = []
         if struct:
            for i in range(3): line = data.readline()
            while line.split()!=[]:
                  coord = line.split()
                  atoms.append(uAtom(coord[0]))
                  Structure.append( map( numpy.float64, coord[1:] ) )
                  line = data.readline()

         Structure = numpy.array(Structure,dtype=numpy.float64)                  
         
         querry = " Origins"
         origins = True
         while 1:
               if querry in line: break
               line = data.readline()
               if not line: 
                  origins = False
                  break              

         if origins:
            for i in range(3): line = data.readline()
            while line.split()!=[]:
                  coord = line.split()
                  Origin.append( map( numpy.float64, coord[1:] ) )
                  line = data.readline()  
            
            Origin = numpy.array(Origin   ,dtype=numpy.float64)

         else:
            Origin    = Structure.copy()
         
         data.close() 
         return dma.DMA( q=numpy.array(ZerothMoments)   ,
                         m=numpy.array(FirstMoments )   , 
                         T=numpy.array(SecondMoments)   ,
                         O=numpy.array(ThirdMoments )   , 
                         H=FourthMoments                ,
                         atoms=atoms              ,
                         pos=Structure            ,
                         origin=Origin,
                         traceless=is_traceless )
    # -----------------------------------------------------------------------------
    elif type.lower() == 'gaussian' or type.lower() == 'gau':
         data = open(file)
         # seek for atomic positions!
         querry = "Electrostatic Properties Using The SCF Density"
         line = data.readline()
         while 1:
               if querry in line: break
               line = data.readline()
               if line =='': raise Exception('No CHELPG population found!')
         for i in range(4): line = data.readline()
         Structure = []
         while ('Atomic Center' in line or 'Ghost Center' in line):
               Structure.append( numpy.array(line.split()[-3:],dtype=numpy.float64)  )
               line = data.readline()

         # seek for charges!
         querry = " Charge="
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(2): line = data.readline()
         ZerothMoments = []
         for i in range(len(Structure)):
             ZerothMoments.append( numpy.float64( line.split()[-1] ) )    
             line = data.readline()
        
         Result = dma.DMA(nfrag=len(Structure), hexadecapoles=False)
         Result.pos = numpy.array(Structure) * units.UNITS.AngstromToBohr
         Result.DMA[0] = numpy.array(ZerothMoments)
         
         data.close()
         return Result

def ParseDMAFromGamessEfpFile(f):
    """parse DMA and their centers from GAMESS *.efp file"""
    STR = []
    chg = []
    dip = []
    qad = []
    oct = []
    d = open(f)
    l = d.readline()
    # STRUCTURE
    while not ("COORDINATES (BOHR)" in l): l = d.readline()
    l = d.readline()
    while not l.startswith(' STOP'):
      STR.append(l.split()[1:4])
      l = d.readline()
    STR = numpy.array(STR,numpy.float64)
    # MONOPOLES
    l = d.readline();l = d.readline()
    while not l.startswith(' STOP'):
       a,b,c = l.split()
       chg.append(numpy.float64(b)+numpy.float64(c))
       l = d.readline()
    chg = numpy.array(chg,numpy.float64)
    # DIPOLES
    l = d.readline();l = d.readline()
    while not l.startswith(' STOP'):
        dip.append(l.split()[1:])
        l = d.readline()
    dip = numpy.array(dip,numpy.float64)
    # QUADRUPOLES
    l = d.readline();l = d.readline()
    while not l.startswith(' STOP'):
         q = l.split()[1:-1]
         l = d.readline()
         q+= l.split()
         qad.append(q)
         l = d.readline()
    qad = numpy.array(qad,numpy.float64)
    # OCTUPOLES
    l = d.readline();l = d.readline()
    while not l.startswith(' STOP'):
         q = l.split()[1:-1]
         l = d.readline()
         q+= l.split()[ :-1]
         l = d.readline()
         q+= l.split()
         oct.append(q)
         l = d.readline()
    oct = numpy.array(oct,numpy.float64)
    return STR, chg, dip, qad, oct

def ParseDistributedPolarizabilitiesFromGamessEfpFile(f):
    """parse distributed polarizabilities and their centers from GAMESS *.efp file"""
    STR = []
    A = []
    d = open(f)
    l = d.readline()
    while not ("POLARIZABLE POINTS" in l): l = d.readline()
    l = d.readline()
    N=0
    while not ("STOP" in l):
      s = l.split()[1:]
      STR.append(s)
      #
      l = d.readline()
      a1 = l.split()
      l = d.readline()
      a2 = l.split()
      l = d.readline()
      a3 = l.split()
      A.append(a1[0])    # XX
      A.append(a1[3])    # XY
      A.append(a2[0])    # XZ
      A.append(a2[2])    # YX
      A.append(a1[1])    # YY
      A.append(a2[1])    # YZ
      A.append(a2[3])    # ZX
      A.append(a3[0])    # ZY
      A.append(a1[2])    # ZZ
      #  
      l = d.readline()
      N+=1

    STR=numpy.array(STR,numpy.float64).reshape(N,3)
    A = numpy.array(A,numpy.float64).reshape(N,3,3)
    return STR, A

def ParseDistributedPolarizabilitiesWrtImFreqFromGamessEfpFile(f):
    """parse distributed polarizabilities wrt imaginary frequency from GAMESS *.efp file"""
    STR = []
    A = []
    d = open(f)
    l = d.readline()
    while not ("DYNAMIC POLARIZABLE POINTS" in l): l = d.readline()
    l = d.readline()
    N=0
    while not ("STOP" in l):
      s = l.split()[2:5]
      STR.append(s)
      #
      l = d.readline()
      a1 = l.split()
      l = d.readline()
      a2 = l.split()
      l = d.readline()
      a3 = l.split()
      A.append(a1[0])    # XX
      A.append(a1[3])    # XY
      A.append(a2[0])    # XZ
      A.append(a2[2])    # YX
      A.append(a1[1])    # YY
      A.append(a2[1])    # YZ
      A.append(a2[3])    # ZX
      A.append(a3[0])    # ZY
      A.append(a1[2])    # ZZ
      #  
      l = d.readline()
      N+=1

    assert not N%12, 'Error in reading DPOLI from EFP file!'
    nmos = N/12
    STR=numpy.array(STR,numpy.float64).reshape(12,nmos,3)
    A = numpy.array(A,numpy.float64).reshape(12,nmos,3,3)
    return STR, A

def ParsePolDerFromFchk(f):
    """parse polarizability derivatives wrt nuclear coordinates from file"""
    g = lambda n: n/5+bool(n%5)
    A = []
    d = open(f)
    l = d.readline()
    while not ("Polarizability Derivatives" in l): l = d.readline()
    N = int(l.split()[-1])
    l = d.readline()
    for i in range(g(N)):
      A+= l.split()
      l = d.readline()

    A = numpy.array(A,numpy.float64).reshape(N,3,3)
    return A

def ParseEnergyFromFchk(file,type='SCF'):
    """parse total energies from g09 fchk given type='SCF' or 'MP2'"""
    data = open(file)
    line = data.readline()
    ### look for energy
    if   type.lower()=='scf': querry = "SCF Energy"
    elif type.lower()=='mp2': querry = "MP2 Energy"
    elif type.lower()=='mp3': querry = "MP3 Energy"
    elif type.lower()=='mp4': querry = "MP4 Energy"
    elif type.lower()=='mp4d': querry = "MP4D Energy"
    elif type.lower()=='mp4dq': querry = "MP4DQ Energy" 
    elif type.lower()=='mp4sdq': querry = "MP4SDQ Energy"
    elif type.lower()=='cc': querry = "Cluster Energy"
    else: raise Exception(" Type <%s> is invalid" % type)
    while querry not in line:
          line = data.readline()
    E = numpy.float64(line.split()[-1])
    data.close()
    return E

def ParseElectronsFromFchk(file):
    """parse number of alpha and beta electrons"""
    data = open(file)
    line = data.readline()
    ### look for energy
    querry1 = 'Number of alpha electrons'
    querry2 = 'Number of beta electrons'
    while querry1 not in line:
          line = data.readline()
    a = int(line.split()[-1])
    line = data.readline()
    b = int(line.split()[-1])
    data.close()
    return a, b


def ParseDipoleMomentFromFchk(file):
    """parse total dipole moment from g09 fchk. 
Note that level of theory is dependent of the keyword used in g09 input file!"""
    data = open(file)
    line = data.readline()
    ### look for dipole moment
    querry = "Dipole Moment"
    while querry not in line:
          line = data.readline()
    line = data.readline()
    D = numpy.float64(line.split())
    data.close()
    return D

def ParseChargesFromFchk(file, type='Mulliken'):
    """Parse the charges from the g09 FCHK file. 
Please specify type if other than 'Mulliken' charges
are to be parsed (by name from FCHK so 'ESP Charges', 'NPA Charges' 
and so on are valid"""
    data = open(file)
    line = data.readline()
    if   type.lower().startswith('mul'): querry = 'Mulliken Charges'
    elif type.lower().startswith('esp'): querry = 'ESP Charges'
    elif type.lower().startswith('npa'): querry = 'NPA Charges'
    elif type.lower().startswith('oni'): querry = 'ONIOM Charges'
    while 1:
        if querry in line: break
        line = data.readline()
    M = int(line.split()[-1])
    line = data.readline()
    G = []
    g = lambda n: n/5+bool(n%5)
    for i in range(g(M)):
        G+= line.split()
        line = data.readline()
    data.close()
    
    G = numpy.array(G,numpy.float64)
    return G

def ParseVecFromFchk(file):
    """parse Ci\mu coeeficients from g09 fchk"""
    data = open(file)
    line = data.readline()
    
    ### look for basis set size
    querry = "Number of basis functions"
    while querry not in line:
          line = data.readline()
    M = int(line.split()[-1])
    
    ### look for MO size
    querry = "Alpha Orbital Energies"
    while querry not in line:
          line = data.readline()
    N = int(line.split()[-1])
    
    ### look for vectors
    querry = "Alpha MO coefficients"
    while querry not in line:
          line = data.readline()
    line = data.readline()
    C = []
    g = lambda n: n/5+bool(n%5)
    for i in range(g(N*M)):
        C+= line.split()
        line = data.readline()
    C = numpy.array(C,dtype=numpy.float64).reshape(N,M)
    data.close()
    return C

def ParseAlphaOrbitalEnergiesFromFchk(file, spin='alpha'):
    """parse alpha or beta orbital energies coeeficients from g09 fchk"""
    if spin!='alpha': raise NotImplementedError, "UTILITIES: Not Implemented beta spins!"

    data = open(file)
    line = data.readline()
    
    ### look for basis set size
    querry = "Number of basis functions"
    while querry not in line:
          line = data.readline()
    M = int(line.split()[-1])
    
    ### look for Alpha MO size
    querry = "Alpha Orbital Energies"
    while querry not in line:
          line = data.readline()
    N = int(line.split()[-1])
    
    line = data.readline()
    C = []
    g = lambda n: n/5+bool(n%5)
    for i in range(g(N)):
        C+= line.split()
        line = data.readline()
    C = numpy.array(C,dtype=numpy.float64)
    data.close()
    return C

def ParseFockFromGamessLog(file,interpol=False):
    """parses Fock matrix from GAMESS log file"""
    data = open(file)
    line = data.readline()
    querry = 'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS ='
    while 1:
       if querry in line: break
       line = data.readline()
    nbasis = int(line.split()[-1])
    data.close()
    data = open(file)
    if interpol:
          querry = 'DIIS INTERPOLATED ALPHA FOCK'
    else: querry = 'TOTAL FOCK OPERATOR'
    for i in range(lind(file,querry)[-1]):
        data.readline()    
    while 1:
       if querry in line: break
       line = data.readline()

    g = lambda n: n/5+bool(n%5)
    fock = []
    fock = numpy.zeros((nbasis,nbasis),dtype=numpy.float64)
    for i in xrange(g(nbasis)):
        line = data.readline()
        line = data.readline()
        nxses= numpy.array(line.split(),int)-1
        line = data.readline()

        for j in xrange(nbasis-i*5):
            line = data.readline()
            line = re_templates.re_dbl_fort_c.sub(r'\1E\2', line)
            ny = int(line.split()[0])-1
            values = line.split()[4:]
            for k in xrange(len(values)):
                nx = nxses[k]
                v  = values[k]
                fock[nx,ny] = v
                fock[ny,nx] = v
    data.close()
    fock = numpy.array(fock,dtype=numpy.float64)
    return fock

def ParseDmatFromFchk(file, basis_size=None, type='SCF'):
    """parses density matrix from Gaussian fchk file"""
        
    data = open(file)
    line = data.readline()
    if   type.lower()=='scf': querry = "Total SCF Density"
    elif type.lower()=='mp2': querry = "Total MP2 Density"
    elif type.lower()=='cc' : querry = "Total CC Density"
    else: raise Exception(" Type <%s> is invalid" % type)

    if basis_size is None:
       ### look for basis set size          
       querry1 = "Number of basis functions"
       while querry1 not in line:
             line = data.readline()
       basis_size = int(line.split()[-1])

    while 1:
        if querry in line: break
        line = data.readline()
    N = int(line.split()[-1])

    line = data.readline()
    dmat = []
    for i in range(int(numpy.ceil(N/5.))): 
        dmat+=[x for x in line.split()] 
        line = data.readline()
    #dmat = numpy.array(dmat,dtype=numpy.float64)
        
    # construct explicit 2D density matrix
    P = numpy.zeros((basis_size,basis_size),dtype=numpy.float64)
    #I = 0
    for i in range(basis_size):
        for j in range(i+1):
            P[i,j] = numpy.float64(dmat.pop(0))#dmat[I]
            P[j,i] = P[i,j] #dmat[I]
            #I += 1
    data.close()
    return numpy.array(P)

def ParseFCFromFchk(file):
    """parses cartesian force constants from Gaussian fchk file"""
        
    data = open(file)
    line = data.readline()
    querry = "Number of atoms"
    while 1:
        if querry in line: break
        line = data.readline()
    N = int(line.split()[-1])
    querry = "Cartesian Force Constants"
    while 1:
        if querry in line: break
        line = data.readline()
    M = int(line.split()[-1])
    line = data.readline()
    FC = []
    g = lambda n: n/5+bool(n%5)
    for i in range(g(M)):
        FC+= line.split()
        line = data.readline()
    data.close()
    
    FC = numpy.array(FC,numpy.float64)
    H = numpy.zeros((N*3,N*3),dtype=numpy.float64)
    I = 0
    for i in xrange(N*3):
        for j in xrange(i+1):
            H[i,j] = FC[I]
            H[j,i] = H[i,j]
            I+=1
    
    return H

def ParseGradFromFchk(file):
    """parses cartesian gradients from Gaussian fchk file"""
        
    data = open(file)
    line = data.readline()
    querry = "Number of atoms"
    while 1:
        if querry in line: break
        line = data.readline()
    N = int(line.split()[-1])
    querry = "Cartesian Gradient"
    while 1:
        if querry in line: break
        line = data.readline()
    M = int(line.split()[-1])
    line = data.readline()
    G = []
    g = lambda n: n/5+bool(n%5)
    for i in range(g(M)):
        G+= line.split()
        line = data.readline()
    data.close()
    
    G = numpy.array(G,numpy.float64)
    return G

def ParseDipoleDerivFromFchk(file):
    """parses dipole derivatives from Gaussian fchk file"""
        
    data = open(file)
    line = data.readline()
    querry = "Dipole Derivatives"
    while 1:
        if querry in line: break
        line = data.readline()
    N = int(line.split()[-1])
    line = data.readline()
    fd = []
    g = lambda n: n/5+bool(n%5)
    for i in range(g(N)):
        fd+= line.split()
        line = data.readline()
    data.close()
    M = N/9
    fd = numpy.array(fd,numpy.float64).reshape(M*3,3)
    return fd

def Parse_EDS_InteractionEnergies(file, method='HF'):
    """parses EDS interaction energies from GAMESS log file"""
    
    data = open(file).read()
    #line = data.readline()
    querry = r'.*VARIATIONAL-PERTURBATIONAL DECOMPOSITION.*'
    if method.lower()=='hf':
       E = ['DE\\(HL\\)'   , 'E\\(EL,10\\)', 'E\\(EL,M,1\\)',
            'E\\(EL,P,1\\)', 'E\\(EX,HL\\)', 'DE\\(DEL,HF\\)', 
            'DE\\(HF\\)']
    elif method.lower()=='mp2':
       E = ['DE\\(HL\\)'   , 'E\\(EL,10\\)', 'E\\(EL,M,1\\)',
            'E\\(EL,P,1\\)', 'E\\(EX,HL\\)', 'DE\\(DEL,HF\\)', 
            'DE\\(HF\\)'   ,
            'E\\(MP,2\\)'  , 'E\\(EL,R,12\\)', 'E\\(EL,M,2\\)', 
            'E\\(EL,P,2\\)', 'E\\(DS,20\\)'  , 'DE\\(EX-DEL,2\\)',
            'DE\\(MP2\\)'  ]
    else: raise ValueError, 'Incorrect method %s specified! Quitting...' % method.upper()
         
    for term in E:
        querry+= '\s*%s\s+(%s).*\n' % (term,re_templates.re_real_e)
    querry = re.compile(querry,re.DOTALL)
    match = re.search(querry,data)
    energies = numpy.array(match.groups(),dtype=numpy.float64)    
    
    return energies

def ParseEFPInteractionEnergies(file):
    """parses EFP interaction energies from GAMESS log file"""
    
    dat  = open(file)
    data = dat.read()
    dat.close()
    
    querry = r'.*FRAGMENT-FRAGMENT INTERACTION ENERGIES.*'
    E = ['ELECTROSTATIC ENERGY','REPULSION ENERGY','POLARIZATION ENERGY',
         'DISPERSION ENERGY','CHARGE TRANSFER ENRGY','FINAL EFP ENERGY',]
         
    for term in E:
        querry+= '\s*%s\s+=\s+(%s).*\n' % (term,re_templates.re_real)
    querry = re.compile(querry,re.DOTALL)
    match = re.search(querry,data)
    energies = numpy.array(match.groups(), dtype=numpy.float64)
    
    return energies

def CalcStep(step,reduced_mass):
    """return step in Angstroms when given step in normal coordinate unit [Bohr*me-1/2] 
    and reduced mass [AMU]"""
    return step * units.UNITS.BohrToAngstrom / numpy.sqrt(units.UNITS.AmuToElectronMass*reduced_mass)

def CalculateCAMM(basis='6-311++G**'): 
    """calculates CAMMs from density matrix from GAUSSIAN09
       using COULOMB.py routines. Usage:
       1) gather all the .fchk and and .log files with pop=chelpg
       2) type ./diff 
       3) the files .camm are creating!!! """
       
    import os, glob, sys
       
    pliki_fchk  = glob.glob('./*_.fchk')
    pliki_fchk.sort()
    pliki_log   = glob.glob('./*_.log')
    pliki_log .sort()    
    print "\n Kolejność plików. Sprawdź czy się zgadzają!\n"  
    for i in range(len(pliki_log)):
        print pliki_log[i], pliki_fchk[i]
    print
       
    for i,file_log in enumerate(pliki_log):
        #dma, fragment = ParseDMA( file_log, 'gaussian' )
        dma = ParseDMA( file_log, 'gaussian' )
        fragment = numpy.array(dma.pos)
        
        frag_file = open('slv.frags','r')
        frag_names = []
        line = frag_file.readline()
        while 1:
            if not line: break
            frag_names.append( line.split()[-1])
            line = frag_file.readline()  

        ### create Molecule object
        structure = []
        for j in range(len(fragment)):
            structure.append( (units.UNITS.atomic_numbers[frag_names[j]],
                                                            fragment[j]) ) 
        molecule = PyQuante.Molecule('mol',
                                      structure,
                                      multiplicity=1,
                                      charge=0,
                                      units='Bohr')
                            
        basis_size = len(PyQuante.Ints.getbasis(molecule,basis))
        print " - basis size= ",basis_size
        dmat = ParseDmatFromFchk(pliki_fchk[i],basis_size)
       
        ### calculate CAMMs                    
        CAMM = multip.MULTIP(molecule=molecule,
                        basis=basis,
                        #basis='sto-3g',
                        method='b3lyp',
                        matrix=dmat,
                        transition=False)
        CAMM.camms()
        CAMM.mmms()
        CAMM.__printMMMs__()
        #CAMM.__printCAMMs__()
       
        result = dma.DMA(nfrag=len(structure))
        result.DMA[0][:] = CAMM.Mon
        #
        result.DMA[1][:] = CAMM.Dip
        #
        result.DMA[2][:,0] = numpy.array(CAMM.Quad)[:,0,0]
        result.DMA[2][:,1] = numpy.array(CAMM.Quad)[:,1,1]
        result.DMA[2][:,2] = numpy.array(CAMM.Quad)[:,2,2]
        result.DMA[2][:,3] = numpy.array(CAMM.Quad)[:,0,1]
        result.DMA[2][:,4] = numpy.array(CAMM.Quad)[:,0,2]
        result.DMA[2][:,5] = numpy.array(CAMM.Quad)[:,1,2]
        #
        result.DMA[3][:,0] = numpy.array(CAMM.Oct)[:,0,0,0]
        result.DMA[3][:,1] = numpy.array(CAMM.Oct)[:,1,1,1]
        result.DMA[3][:,2] = numpy.array(CAMM.Oct)[:,2,2,2]
        result.DMA[3][:,3] = numpy.array(CAMM.Oct)[:,0,0,1]
        result.DMA[3][:,4] = numpy.array(CAMM.Oct)[:,0,0,2]
        result.DMA[3][:,5] = numpy.array(CAMM.Oct)[:,0,1,1]
        result.DMA[3][:,6] = numpy.array(CAMM.Oct)[:,1,1,2]
        result.DMA[3][:,7] = numpy.array(CAMM.Oct)[:,0,2,2]
        result.DMA[3][:,8] = numpy.array(CAMM.Oct)[:,1,2,2]
        result.DMA[3][:,9] = numpy.array(CAMM.Oct)[:,0,1,2]
        #
        #print result
        out = open(file_log[:-4]+'.camm','w')
        out.write(str(result))
        out.close()
        print " Writing file:  :", file_log[:-4]+'.camm'
    print

def RotationMatrix(initial=None,final=None):
    """returns rotation matrix and rms from SVD superposition of two structures.
    The initial structure is rotated into final one. The transformation is defined as follows:
    final = numpy.dot(initial, rot) + transl
    Returns: rot, rms"""
    sup = SVDSuperimposer()
    sup.set(final,initial)
    sup.run()
    rms = sup.get_rms()
    rot, transl = sup.get_rotran()
    return rot, rms

def Periodic(mendeleiev):
    '''Returns the mendeleiev table as a python list of tuples. Each cell
    contains either None or a tuple (symbol, atomic number), or a list of pairs
    for the cells * and **. Requires: "import re". Source: Gribouillis at
    www.daniweb.com - 2008 '''

    # L is a consecutive list of tuples ('Symbol', atomic number)
    L = [ (e,i+1) for (i,e) in enumerate( re.compile ("[A-Z][a-z]*").findall('''
    HHeLiBeBCNOFNeNaMgAlSiPSClArKCaScTiVCrMnFeCoNiCuZnGaGeAsSeBrKr
    RbSrYZrNbMoTcRuRhPdAgCdInSnSbTeIXeCsBaLaCePrNdPmSmEuGdTbDyHoEr
    TmYbLuHfTaWReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaUNpPuAmCmBkCfEsFm
    MdNoLrRfDbSgBhHsMtDsRgUubUutUuqUupUuhUusUuo'''))]

    # The following fills the void with nones and returns the list of lists
    mendeleiev = 0

    if mendeleiev:
        for i,j in ( (88,103), (56,71) ):
            L[i] = L[i:j]
            L[i+1:] = L[j:]
        for i,j in ( (12,10 ), (4,10), (1,16) ):
            L[i:i]=[None]*j 

        return [ L[18*i:18*(i+1)] for i in range (7) ]

    # Return a plain list of tuples
    else:
        return L

def Atomn(s,ptable):
    '''Returns the atomic number based on atomic symbol string
    ptable is a list of consecutive (symbol, atomic number) tuples.'''
    for n,a in enumerate(ptable):
        if a[0].lower().find(s.strip().lower()) !=-1 :
            return float(n+1)

def ParseUnitedAtoms(string):
    """parses united atom indices given in a format:
       i,j,k,...[-m,n,o,...-...] where first indice 
       denotes UA atom and the others before '-' are
       atoms to be contracted. Returns list of lists
       of UA atoms"""
       
    string = string.split('-')
    for i in range(len(string)):
        string[i] = map(int, string[i].split(','))
    return string

def DMAMadrixMultiply(matrix,dma_list):
    """Returns the output of matrix multiplication
    between DMA objects and a matrix"""
    
    # number of cartesian coordinates
    N = len(dma_list)
    positions = dma_list[0].pos
    # number of centers
    K = len(dma_list[0][0])
    charges = numpy.zeros((N ,K  ),dtype=numpy.float64)
    dipoles = numpy.zeros((3 ,N,K),dtype=numpy.float64)
    qdrples = numpy.zeros((6 ,N,K),dtype=numpy.float64)
    octples = numpy.zeros((10,N,K),dtype=numpy.float64)
    for i,dmai in enumerate(dma_list):
        charges[i,:]   = dmai.DMA[0]
        dipoles[:,i,:] = numpy.transpose(dmai.DMA[1])
        qdrples[:,i,:] = numpy.transpose(dmai.DMA[2])
        octples[:,i,:] = numpy.transpose(dmai.DMA[3])
 
    ### TRANSFORMATION!    
    charges = numpy.dot(matrix,charges)
    dipoles = numpy.dot(matrix,dipoles)
    qdrples = numpy.dot(matrix,qdrples)
    octples = numpy.dot(matrix,octples)
    
    result = []
    for i in range(len(matrix)):
        dmai = dma.DMA(nfrag=K)
        dmai.DMA[0] = charges[i]
        dmai.DMA[1] = numpy.transpose(dipoles[i])
        dmai.DMA[2] = numpy.transpose(qdrples[i])
        dmai.DMA[3] = numpy.transpose(octples[i])
        dmai.pos = positions
        result.append( dmai )
        
    return result

### COMPILATION CODES 
code = """
converter=UNITS.HartreePerHbarToCmRec
#
dma1=DMA1.copy()
dma2=DMA2.copy()
# make FULL format of DMA distribution
dma1.MAKE_FULL() # hexadecapole integrals not implemented yet
dma2.MAKE_FULL()
# transform FULL format to fraceless forms for quadrupoles and octupoles
dma1.MakeTraceless()
dma2.MakeTraceless()
#
Ra,qa,Da,Qa,Oa = dma1.DMA_FULL
Rb,qb,Db,Qb,Ob = dma2.DMA_FULL
#
qq = 0
qD = 0 ; Dq = 0
qQ = 0 ; Qq = 0
DQ = 0 ; QD = 0
QQ = 0 
qO = 0 ; Oq = 0
DD = 0 
DO = 0 ; OD = 0
QO = 0 ; OQ = 0
OO = 0 ;
qH = 0 ; Hq = 0
#
Tensordot = numpy.tensordot
Dot = numpy.dot
for i in xrange(len(Ra)):
         for j in xrange(len(Rb)):
            R    = Rb[j]-Ra[i]
            Rab=numpy.sqrt(numpy.sum(R**2,axis=0))
            if (Rab < threshold and Rab !=0):
             qq  +=   qa[i]*qb[j]/Rab                                                               # qa - qb  | R1
             #if not hash:
             qD  +=  -qa[i]*Tensordot(Db[j],R,(0,0))/Rab**3                                         # qa - Db  | R2
             Dq  +=  +qb[j]*Tensordot(Da[i],R,(0,0))/Rab**3                                         # qb - Da  | R2
             DD  +=-3*Tensordot(Da[i],R,(0,0))*Tensordot(Db[j],R,(0,0))/Rab**5                      # Da - Db  | R3
             DD  +=   Tensordot(Da[i],Db[j],(0,0))/Rab**3                                           # Da - Db  | R3
             qQ  +=   qa[i]*Tensordot(R,Tensordot(Qb[j],R,(0,0)),(0,0))/Rab**5                      # qa - Qb  | R3
             Qq  +=   qb[j]*tensordot(R,Tensordot(Qa[i],R,(0,0)),(0,0))/Rab**5                      # qb - Qa  | R3
             DQ  +=-2*Tensordot(Da[i],Tensordot(Qb[j],R,(0,0)),(0,0))/Rab**5                        # Da - Qb  | R4
             QD  += 2*Tensordot(Db[j],Tensordot(Qa[i],R,(0,0)),(0,0))/Rab**5                        # Db - Qa  | R4
             DQ  += 5*Tensordot(Da[i],R,(0,0))*Tensordot(R,Tensordot(Qb[j],R,(0,0)),(0,0))/Rab**7   # Da - Qb  | R4
             QD  +=-5*Tensordot(Db[j],R,(0,0))*Tensordot(R,Tensordot(Qa[i],R,(0,0)),(0,0))/Rab**7   # Db - Qa  | R4
             qO  +=  -qa[i]*Tensordot(R,Tensordot(R,Tensordot(Ob[j],R,(0,0)),(0,0)),(0,0))/Rab**7   # qa - Ob  | R4
             Oq  +=   qb[j]*Tensordot(R,Tensordot(R,Tensordot(Oa[i],R,(0,0)),(0,0)),(0,0))/Rab**7   # qb - Oa  | R4
             QQ  += (35.)/(3.)* (Tensordot(R,Tensordot(Qa[i],R,(0,0)),(0,0)) *
                                 Tensordot(R,Tensordot(Qb[j],R,(0,0)),(0,0))  ) / Rab**9            # Qa - Qb  | R5
             OD  +=-7*(Tensordot(Db[j],R,(0,0)) *
                       Tensordot(R,Tensordot(R,Tensordot(Oa[i],R,(0,0)),(0,0)),(0,0)) ) / Rab**9    # Db - Oa  | R5
             DO  +=-7*(Tensordot(Da[i],R,(0,0)) *
                       Tensordot(R,Tensordot(R,Tensordot(Ob[j],R,(0,0)),(0,0)),(0,0)) ) / Rab**9    # Da - Ob  | R5
             QQ  +=-(20.)/(3.) * Tensordot(Tensordot(R,Qa[i],(0,0)),
                                           Tensordot(R,Qb[j],(0,0)),(0,0)) / Rab**7                 # Qa - Qb  | R5
             QQ  +=(2.)/(3.)  * Tensordot(Qa[i],Qb[j])  / Rab**5                                    # Qa - Qb  | R5
             OD  +=3 * Tensordot(R,Tensordot(R,Tensordot(Oa[i],Db[j],(0,0)),(0,0)),(0,0)) / Rab**7  # Db - Oa  | R5
             DO  +=3 * Tensordot(R,Tensordot(R,Tensordot(Ob[j],Da[i],(0,0)),(0,0)),(0,0)) / Rab**7  # Da - Ob  | R5
             ### The remaining terms with hexadecapoles are not implemented yet
             #Eint+= qb[j] * Tensordot(R,Tensordot(R,Tensordot(R,Tensordot(R,Ha[i],
             #                (0,0)),(0,0)),(0,0)),(0,0))   / Rab**9                                 # Ha - qb  | R5
             #Eint+= qa[i] * Tensordot(R,Tensordot(R,Tensordot(R,Tensordot(R,Hb[j],
             #                (0,0)),(0,0)),(0,0)),(0,0))   / Rab**9                                 # Hb - qj  | R5
             ### these are implemented already !
             #OQ  += 2* Tensordot(Tensordot(Oa[i],Qb[j],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qb - Oa  | R6
             #QO  +=-2* Tensordot(Tensordot(Ob[j],Qa[i],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qa - Ob  | R6
             #OQ  +=-14*Tensordot(Tensordot(R,Tensordot(Oa[i],R,(1,0)),(0,0)) ,                      # Qb - Oa  | R6
             #                    Tensordot(R,Qb[j],(0,0)) ,(0,0)) / Rab**9                          
             #QO  += 14*Tensordot(Tensordot(R,Tensordot(Ob[j],R,(1,0)),(0,0)) ,                      # Qa - Ob  | R6
             #                    Tensordot(R,Qa[i],(0,0)) ,(0,0)) / Rab**9
             #OQ  +=( 21*Tensordot(Tensordot(R,Tensordot(Oa[i],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
             #         * Tensordot(R,Tensordot(Qb[j],R,(0,0)),(0,0))) / Rab**11
             #QO  +=(-21*Tensordot(Tensordot(R,Tensordot(Ob[j],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
             #         * Tensordot(R,Tensordot(Qa[i],R,(0,0)),(0,0))) / Rab**11   
             #OO  +=(2.)/(5.)*Tensordot(Oa[i],Ob[j],((0,1,2),(0,1,2))) / Rab**7                      # Ob - Oa  | R7
             #OO  +=(-42./5.)*Tensordot(Tensordot(R,Oa[i],(0,0)),
             #                          Tensordot(R,Ob[j],(0,0)),
             #                          ((0,1),(0,1))) / Rab**9                                      # Ob - Oa  | R7
             #OO  +=(189.)/(5.)*Tensordot(
             #                            Tensordot(Tensordot(R,Oa[i],(0,0)),R,(0,0)),
             #                            Tensordot(Tensordot(R,Ob[j],(0,0)),R,(0,0)),
             #                            (0,0)) /Rab**11
             #OO  +=-(231./5.)*(Tensordot(Tensordot(Tensordot(R,Oa[i],(0,0)),R,(0,0)),R,(0,0)) *
             #                  Tensordot(Tensordot(Tensordot(R,Ob[j],(0,0)),R,(0,0)),R,(0,0)) ) /\
             #                  Rab**13
             
             Eint = qq + qD + Dq + qQ + Qq + qO + Oq + DD + DQ + QD + DO + OD + QQ + QO + OQ + OO
             
             ### save the partitioning for current usage
             Emtp.qq = qq;Emtp.qD = qD;Emtp.qQ = qQ;Emtp.qO = qO;Emtp.QO = QO;
             Emtp.DD = DD;Emtp.DQ = DQ;Emtp.DO = DO;Emtp.QQ = QQ;Emtp.OO = OO;
#
Emtp.A = (         qq )
Emtp.B = (Emtp.A + qD + Dq )
Emtp.C = (Emtp.B + DD + qQ + Qq )
Emtp.D = (Emtp.C + qO + Oq + DQ + QD)
Emtp.E = (Emtp.D + DO + OD + QQ + qH + Hq)
Emtp.F = (Emtp.A + Dq + qD + DD )
Emtp.G = (Emtp.F + Qq + qQ + QD + DQ + QQ)
Emtp.H = (Emtp.G + Oq + qO + OD + DO + OQ + QO + OO)
#
Emtp.A *= converter
Emtp.B *= converter
Emtp.C *= converter
Emtp.D *= converter
Emtp.E *= converter
"""
command = compile(code,'<string>','exec')             

def Emtp(DMA1,DMA2):
    """calculates E(EL)MTP from two DMA distributions.
dma1 and dma2 are the objects of the class DMA. Calculations are
in atomic units and a respective interaction energy is in 
a.u. as well. Uses FORTRAN subroutine CLEMTP"""
    
    converter=units.UNITS.HartreePerHbarToCmRec
    #
    dma1=DMA1.copy()
    dma2=DMA2.copy()
    # make FULL format of DMA distribution
    dma1.MAKE_FULL() # hexadecapole integrals not implemented yet
    dma2.MAKE_FULL()
    # transform FULL format to fraceless forms for quadrupoles and octupoles
    dma1.MakeTraceless()
    dma2.MakeTraceless()
    #
    Ra,qa,Da,Qa,Oa = dma1.DMA_FULL
    Rb,qb,Db,Qb,Ob = dma2.DMA_FULL
    #
    Eint,A,B,C,D,E,CC,CD,CQ,CO,DD,DQ=qm.clemtp.clemtp(Ra,qa,Da,Qa,Oa,Rb,qb,Db,Qb,Ob)
    DO = 0; QQ = 0; QO = 0; OO = 0
    #
    Emtp.A = A
    Emtp.B = B
    Emtp.C = C
    Emtp.D = D
    Emtp.E = E
    Emtp.F = 0
    Emtp.G = 0
    Emtp.H = 0
    #
    Emtp.A *= converter
    Emtp.B *= converter
    Emtp.C *= converter
    Emtp.D *= converter
    Emtp.E *= converter
    #
    log = "\n" 
    log+= " --------------------------------:--------------------------\n"
    log+= " INTERACTION ENERGY TERMS [cm-1] : PARTITIONING TERMS [cm-1]\n"
    log+= " --------------------------------:--------------------------\n"
    log+= "%6s %20.2f      :\n" % ("Total".rjust(6),Eint*converter)
    log+= " "+"-"*32+":"+"-"*26+"\n"
    log+= "%7s %19.2f      :  1        %10.2f\n" % ("q-q".rjust(6), CC *converter,Emtp.A)
    log+= "%7s %19.2f      :  1+2      %10.2f\n" % ("q-D".rjust(6), CD *converter,Emtp.B)
    log+= "%7s %19.2f      :  1+2+3    %10.2f\n" % ("q-Q".rjust(6), CQ *converter,Emtp.C)
    log+= "%7s %19.2f      :  1+2+3+4  %10.2f\n" % ("q-O".rjust(6), CO *converter,Emtp.D)
    log+= "%7s %19.2f      :  1+2+3+4+5%10.2f\n" % ("D-D".rjust(6), DD *converter,Emtp.E)
    log+= "%7s %19.2f      :\n"                  % ("D-Q".rjust(6), DQ *converter)
    log+= "%7s %19.2f      :\n"                  % ("D-O".rjust(6), DO *converter)
    log+= "%7s %19.2f      :\n"                  % ("Q-Q".rjust(6), QQ *converter)
    log+= "%7s %19.2f      :\n"                  % ("Q-O".rjust(6), QO *converter)
    log+= "%7s %19.2f      :\n"                  % ("O-O".rjust(6), OO *converter)
    log+= " "+"-"*32+":"+"-"*26+"\n"
    log+= "\n"
    Emtp.log = log

    return Emtp.A, Emtp.B, Emtp.C, Emtp.D, Emtp.E

def Emtp_charges(DMA1,DMA2):
    """calculates E(EL)MTP from two DMA CHARGE distributions.
dma1 and dma2 are the objects of the class DMA. Calculations are
in atomic units and a respective interaction energy is in 
a.u. as well."""
    
    converter=units.UNITS.HartreePerHbarToCmRec
    #
    dma1=DMA1.copy()
    dma2=DMA2.copy()
    #
    Ra,Rb = dma1.get_origin(),dma2.get_origin()
    qa,qb = dma1[0],dma2[0]
    #
    R = Ra[:,:,numpy.newaxis] - Rb.T[numpy.newaxis,:,:]
    R = R*R; R = numpy.sum(R, axis=1)
    R = numpy.sqrt(R)
    qq = numpy.sum(numpy.outer(qa,qb)/R,axis=None)
    #
    qq *= converter
    Emtp_charges.A = qq
    Emtp_charges.B = qq
    Emtp_charges.C = qq
    Emtp_charges.D = qq
    Emtp_charges.E = qq
    Emtp_charges.F = 0
    Emtp_charges.G = 0
    Emtp_charges.H = 0
    #
    #Q = qq*converter
    #log = "\n" 
    #log+= " --------------------------------:--------------------------\n"
    #log+= " INTERACTION ENERGY TERMS [cm-1] : PARTITIONING TERMS [cm-1]\n"
    #log+= " --------------------------------:--------------------------\n"
    #log+= "%6s %20.2f      :\n" % ("Total".rjust(6),Q)
    #log+= " "+"-"*32+":"+"-"*26+"\n"
    #log+= "%7s %19.2f      :  1        %10.2f\n" % ("q-q".rjust(6), Q,Q)
    #log+= "%7s %19.2f      :  1+2      %10.2f\n" % ("q-D".rjust(6), 0,Q)
    #log+= "%7s %19.2f      :  1+2+3    %10.2f\n" % ("q-Q".rjust(6), 0,Q)
    #log+= "%7s %19.2f      :  1+2+3+4  %10.2f\n" % ("q-O".rjust(6), 0,Q)
    #log+= "%7s %19.2f      :  1+2+3+4+5%10.2f\n" % ("D-D".rjust(6), 0,Q)
    #log+= "%7s %19.2f      :\n"                  % ("D-Q".rjust(6), 0)
    #log+= "%7s %19.2f      :\n"                  % ("D-O".rjust(6), 0)
    #log+= "%7s %19.2f      :\n"                  % ("Q-Q".rjust(6), 0)
    #log+= "%7s %19.2f      :\n"                  % ("Q-O".rjust(6), 0)
    #log+= "%7s %19.2f      :\n"                  % ("O-O".rjust(6), 0)
    #log+= " "+"-"*32+":"+"-"*26+"\n"
    log = "\n"
    Emtp_charges.log = log

    return qq,qq,qq,qq,qq
    
def get_elmtp(DMA1,DMA2, return_all=False):
    """calculates E(EL)MTP from two DMA distributions.
dma1 and dma2 are the objects of the class DMA. Calculations are
in atomic units and a respective interaction energy is in 
a.u. as well. """

    converter=units.UNITS.HartreePerHbarToCmRec
    #
    dma1=DMA1.copy()
    dma2=DMA2.copy()
    # make FULL format of DMA distribution
    dma1.MAKE_FULL() # hexadecapole integrals not implemented yet
    dma2.MAKE_FULL()
    # transform FULL format to fraceless forms for quadrupoles, octupoles and hexadecapoles
    dma1.MakeTraceless()
    dma2.MakeTraceless()
    #
    if dma1.has_hexadecapoles and dma2.has_hexadecapoles:
       hexadecapoles = True
       Ra,qa,Da,Qa,Oa,Ha = dma1.DMA_FULL
       Rb,qb,Db,Qb,Ob,Hb = dma2.DMA_FULL
    else:
       hexadecapoles = False
       Ra,qa,Da,Qa,Oa = dma1.DMA_FULL 
       Rb,qb,Db,Qb,Ob = dma2.DMA_FULL
    #
    qq = 0
    qD = 0 ; Dq = 0
    qQ = 0 ; Qq = 0
    DQ = 0 ; QD = 0
    QQ = 0 
    qO = 0 ; Oq = 0
    DD = 0 
    DO = 0 ; OD = 0
    QO = 0 ; OQ = 0
    OO = 0 ;
    qH = 0 ; Hq = 0
    #
    Tensordot = numpy.tensordot
    for i in xrange(len(Ra)):
         for j in xrange(len(Rb)):
             R    = Rb[j]-Ra[i]
             Rab=numpy.sqrt(numpy.sum(R**2,axis=0))
             qq  +=   qa[i]*qb[j]/Rab                                                               # qa - qb  | R1
             qD  +=  -qa[i]*Tensordot(Db[j],R,(0,0))/Rab**3                                         # qa - Db  | R2
             Dq  +=  +qb[j]*Tensordot(Da[i],R,(0,0))/Rab**3                                         # qb - Da  | R2
             DD  +=-3*Tensordot(Da[i],R,(0,0))*Tensordot(Db[j],R,(0,0))/Rab**5                      # Da - Db  | R3
             DD  +=   Tensordot(Da[i],Db[j],(0,0))/Rab**3                                           # Da - Db  | R3
             qQ  +=   qa[i]*Tensordot(R,Tensordot(Qb[j],R,(0,0)),(0,0))/Rab**5                      # qa - Qb  | R3
             Qq  +=   qb[j]*Tensordot(R,Tensordot(Qa[i],R,(0,0)),(0,0))/Rab**5                      # qb - Qa  | R3
             DQ  +=-2*Tensordot(Da[i],Tensordot(Qb[j],R,(0,0)),(0,0))/Rab**5                        # Da - Qb  | R4
             QD  += 2*Tensordot(Db[j],Tensordot(Qa[i],R,(0,0)),(0,0))/Rab**5                        # Db - Qa  | R4
             DQ  += 5*Tensordot(Da[i],R,(0,0))*Tensordot(R,Tensordot(Qb[j],R,(0,0)),(0,0))/Rab**7   # Da - Qb  | R4
             QD  +=-5*Tensordot(Db[j],R,(0,0))*Tensordot(R,Tensordot(Qa[i],R,(0,0)),(0,0))/Rab**7   # Db - Qa  | R4
             qO  +=  -qa[i]*Tensordot(R,Tensordot(R,Tensordot(Ob[j],R,(0,0)),(0,0)),(0,0))/Rab**7   # qa - Ob  | R4
             Oq  +=   qb[j]*Tensordot(R,Tensordot(R,Tensordot(Oa[i],R,(0,0)),(0,0)),(0,0))/Rab**7   # qb - Oa  | R4
             QQ  += (35.)/(3.)* (Tensordot(R,Tensordot(Qa[i],R,(0,0)),(0,0)) *
                                 Tensordot(R,Tensordot(Qb[j],R,(0,0)),(0,0))  ) / Rab**9            # Qa - Qb  | R5
             OD  +=-7*(Tensordot(Db[j],R,(0,0)) *
                       Tensordot(R,Tensordot(R,Tensordot(Oa[i],R,(0,0)),(0,0)),(0,0)) ) / Rab**9    # Db - Oa  | R5
             DO  +=-7*(Tensordot(Da[i],R,(0,0)) *
                       Tensordot(R,Tensordot(R,Tensordot(Ob[j],R,(0,0)),(0,0)),(0,0)) ) / Rab**9    # Da - Ob  | R5
             QQ  +=-(20.)/(3.) * Tensordot(Tensordot(R,Qa[i],(0,0)),
                                           Tensordot(R,Qb[j],(0,0)),(0,0)) / Rab**7                 # Qa - Qb  | R5
             QQ  +=(2.)/(3.)  * Tensordot(Qa[i],Qb[j])  / Rab**5                                    # Qa - Qb  | R5
             OD  +=3 * Tensordot(R,Tensordot(R,Tensordot(Oa[i],Db[j],(0,0)),(0,0)),(0,0)) / Rab**7  # Db - Oa  | R5
             DO  +=3 * Tensordot(R,Tensordot(R,Tensordot(Ob[j],Da[i],(0,0)),(0,0)),(0,0)) / Rab**7  # Da - Ob  | R5
             if hexadecapoles:
                Hq  += qb[j] * Tensordot(R,Tensordot(R,Tensordot(R,Tensordot(R,Ha[i],                                  
                                (0,0)),(0,0)),(0,0)),(0,0))   / Rab**9                              # Ha - qb  | R5
                qH  += qa[i] * Tensordot(R,Tensordot(R,Tensordot(R,Tensordot(R,Hb[j],
                                (0,0)),(0,0)),(0,0)),(0,0))   / Rab**9                              # Hb - qa  | R5
             ### these are implemented already !
             OQ  += 2* Tensordot(Tensordot(Oa[i],Qb[j],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qb - Oa  | R6
             QO  +=-2* Tensordot(Tensordot(Ob[j],Qa[i],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qa - Ob  | R6
             OQ  +=-14*Tensordot(Tensordot(R,Tensordot(Oa[i],R,(1,0)),(0,0)) ,                      # Qb - Oa  | R6
                                 Tensordot(R,Qb[j],(0,0)) ,(0,0)) / Rab**9                          
             QO  += 14*Tensordot(Tensordot(R,Tensordot(Ob[j],R,(1,0)),(0,0)) ,                      # Qa - Ob  | R6
                                 Tensordot(R,Qa[i],(0,0)) ,(0,0)) / Rab**9
             OQ  +=( 21*Tensordot(Tensordot(R,Tensordot(Oa[i],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
                      * Tensordot(R,Tensordot(Qb[j],R,(0,0)),(0,0))) / Rab**11
             QO  +=(-21*Tensordot(Tensordot(R,Tensordot(Ob[j],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
                      * Tensordot(R,Tensordot(Qa[i],R,(0,0)),(0,0))) / Rab**11   
             OO  +=(2.)/(5.)*Tensordot(Oa[i],Ob[j],((0,1,2),(0,1,2))) / Rab**7                      # Ob - Oa  | R7
             OO  +=(-42./5.)*Tensordot(Tensordot(R,Oa[i],(0,0)),
                                       Tensordot(R,Ob[j],(0,0)),
                                       ((0,1),(0,1))) / Rab**9                                      # Ob - Oa  | R7
             OO  +=(189.)/(5.)*Tensordot(
                                         Tensordot(Tensordot(R,Oa[i],(0,0)),R,(0,0)),
                                         Tensordot(Tensordot(R,Ob[j],(0,0)),R,(0,0)),
                                         (0,0)) /Rab**11                                            # Oa - Ob  | R7
             OO  +=-(231./5.)*(Tensordot(Tensordot(Tensordot(R,Oa[i],(0,0)),R,(0,0)),R,(0,0)) *
                               Tensordot(Tensordot(Tensordot(R,Ob[j],(0,0)),R,(0,0)),R,(0,0)) ) /\
                               Rab**13                                                              # Oa - Ob  | R7
             
             Eint = qq + qD + Dq + qQ + Qq + qO + Oq + DD + DQ + QD + DO + OD + QQ + QO + OQ + OO + qH + Hq
             
             ### save the partitioning for current usage
             get_elmtp.qq = qq;get_elmtp.qD = qD;get_elmtp.qQ = qQ;get_elmtp.qO = qO;get_elmtp.QO = QO;get_elmtp.qH = qH
             get_elmtp.DD = DD;get_elmtp.DQ = DQ;get_elmtp.DO = DO;get_elmtp.QQ = QQ;get_elmtp.OO = OO;get_elmtp.Hq = Hq
    #
    get_elmtp.A = (         qq )
    get_elmtp.B = (get_elmtp.A + qD + Dq )
    get_elmtp.C = (get_elmtp.B + DD + qQ + Qq )
    get_elmtp.D = (get_elmtp.C + qO + Oq + DQ + QD)
    get_elmtp.E = (get_elmtp.D + DO + OD + QQ + qH + Hq)
    get_elmtp.F = (get_elmtp.A + Dq + qD + DD )
    get_elmtp.G = (get_elmtp.F + Qq + qQ + QD + DQ + QQ)
    get_elmtp.H = (get_elmtp.G + Oq + qO + OD + DO + OQ + QO + OO)
    #
    get_elmtp.A *= converter
    get_elmtp.B *= converter
    get_elmtp.C *= converter
    get_elmtp.D *= converter
    get_elmtp.E *= converter
    get_elmtp.F *= converter
    get_elmtp.G *= converter
    get_elmtp.H *= converter
    #
    log = "\n" 
    log+= " --------------------------------:--------------------------\n"
    log+= " INTERACTION ENERGY TERMS [cm-1] : PARTITIONING TERMS [cm-1]\n"
    log+= " --------------------------------:--------------------------\n"
    log+= "%6s %20.2f      :\n" % ("Total".rjust(6),Eint*converter)
    log+= " "+"-"*32+":"+"-"*26+"\n"
    log+= "%7s %19.2f      :  1        %10.2f\n" % ("q-q".rjust(6), qq    *converter,get_elmtp.A)
    log+= "%7s %19.2f      :  1+2      %10.2f\n" % ("q-D".rjust(6),(qD+Dq)*converter,get_elmtp.B)
    log+= "%7s %19.2f      :  1+2+3    %10.2f\n" % ("q-Q".rjust(6),(qQ+Qq)*converter,get_elmtp.C)
    log+= "%7s %19.2f      :  1+2+3+4  %10.2f\n" % ("q-O".rjust(6),(qO+Oq)*converter,get_elmtp.D)
    log+= "%7s %19.2f      :  1+2+3+4+5%10.2f\n" % ("D-D".rjust(6), DD    *converter,get_elmtp.E)
    log+= "%7s %19.2f      :\n"                  % ("D-Q".rjust(6),(DQ+QD)*converter)
    log+= "%7s %19.2f      :\n"                  % ("D-O".rjust(6),(DO+OD)*converter)
    log+= "%7s %19.2f      :\n"                  % ("Q-Q".rjust(6), QQ    *converter)
    log+= "%7s %19.2f      :\n"                  % ("Q-O".rjust(6),(QO+OQ)*converter)
    log+= "%7s %19.2f      :\n"                  % ("O-O".rjust(6), OO    *converter)
    log+= "%7s %19.2f      :\n"                  % ("q-H".rjust(6),(qH+Hq)*converter)
    log+= " "+"-"*32+":"+"-"*26+"\n"
    log+= "\n"
    get_elmtp.log = log
    
    if return_all: return get_elmtp.A, get_elmtp.B, get_elmtp.C, get_elmtp.D, get_elmtp.E, get_elmtp.F, get_elmtp.G, get_elmtp.H
    else:          return get_elmtp.A, get_elmtp.B, get_elmtp.C, get_elmtp.D, get_elmtp.E

def Emtpc(DMA1,DMA2,threshold=1000,hash=True):
    """Compiled version of EMtp. calculates E(EL)MTP from two DMA distributions.
dma1 and dma2 are the objects of the class DMA. Calculations are
in atomic units and a respective interaction energy is in 
a.u. as well. """
    exec(command)
    #
    log = "\n" 
    log+= " --------------------------------:--------------------------\n"
    log+= " INTERACTION ENERGY TERMS [cm-1] : PARTITIONING TERMS [cm-1]\n"
    log+= " --------------------------------:--------------------------\n"
    log+= "%6s %20.2f      :\n" % ("Total".rjust(6),Eint*converter)
    log+= " "+"-"*32+":"+"-"*26+"\n"
    log+= "%7s %19.2f      :  1        %10.2f\n" % ("q-q".rjust(6), qq    *converter,Emtp.A)
    log+= "%7s %19.2f      :  1+2      %10.2f\n" % ("q-D".rjust(6),(qD+Dq)*converter,Emtp.B)
    log+= "%7s %19.2f      :  1+2+3    %10.2f\n" % ("q-Q".rjust(6),(qQ+Qq)*converter,Emtp.C)
    log+= "%7s %19.2f      :  1+2+3+4  %10.2f\n" % ("q-O".rjust(6),(qO+Oq)*converter,Emtp.D)
    log+= "%7s %19.2f      :  1+2+3+4+5%10.2f\n" % ("D-D".rjust(6), DD    *converter,Emtp.E)
    log+= "%7s %19.2f      :\n"                  % ("D-Q".rjust(6),(DQ+QD)*converter)
    log+= "%7s %19.2f      :\n"                  % ("D-O".rjust(6),(DO+OD)*converter)
    log+= "%7s %19.2f      :\n"                  % ("Q-Q".rjust(6), QQ    *converter)
    log+= "%7s %19.2f      :\n"                  % ("Q-O".rjust(6),(QO+OQ)*converter)
    log+= "%7s %19.2f      :\n"                  % ("O-O".rjust(6), OO    *converter)
    log+= " "+"-"*32+":"+"-"*26+"\n"
    log+= "\n"
    Emtp.log = log
    
    return Emtp.A, Emtp.B, Emtp.C, Emtp.D, Emtp.E

def FrequencyShiftPol(solvent,solpol,point):
    """calculates polarization contribution to frequency shift from simple
    1-center model. Solvent object has to be a traceless property!
    Returns frequency shift in [cm-1]"""
    
    field = ElectricField(solvent,point,is_full=True)
    print field
    shift = numpy.tensordot(field,numpy.tensordot(solpol,field,(0,0)),(0,0))
    shift*= -1./2.
    
    return shift * units.UNITS.HartreePerHbarToCmRec

def FrequencyShift(solute=0,solvent=0,solute_structure=0):
    """calculates frequency shift of solute (MCHO instance)
    as well as solvent (tuple of DMA instances)."""
    # all calculations are performed in atomic units.
    # The final cm-1 unit is obtained at the very end
    # of the computations.
    new = solute.copy()
    #new.pos = array(solute_structure)
    new.set_structure(pos=solute_structure,equal=True)
    A,B,C,D,E = Emtp(new,solvent.copy())
    #A,B,C,D,E = Emtp_charges(new,solvent.copy())
    result = numpy.array([A,B,C,D,E])
    # switch to cm-1
    # result  *= units.UNITS.HartreePerHbarToCmRec
    return result

class Allign:
    """
 Rotates and translates the molecule to allign it with respect to the new user-specified global coordinate frame. 
 
 Usage:

 a = Allign(xyz, atid=[], vec=[], axes=(0,1,2), dma=None)
 dma, xyz = a.get_transformed()

 where:
  - xyz       - ndarray of shape (natoms,3)
  - atid      - tuple/list of shape (3): atomic indices (not-Python convention). They specify the coordinate system.
                First atom : new origin 
                Second atom: z-axis
                Third atom : probably xz-plane
                If atid(Third atom)<1: compute this axis orthogonal to the two previous ones. 
  - vec       - ndarray of shape (3,3): alternative specification of coordinate system. 
                Rows of vec determine the axes of new coordinate system so they should be orthogonal to each other.
  - axes      - tuple/list of shape (3): type of allignment. Specifies the pivoting of the axes. 
  - dma       - libbbg.dma.DMA object. If specified, this object will be also alligned to the new orientation.

                                                                             Last Revision: 10 Jan 2015
"""
    def __init__(self,xyz=numpy.zeros(3),atid=[],vec=[],axes=(0,1,2),dma=None):
        self.xyz=xyz
        self.atid=atid
        self.vec=vec
        self.axes=axes
        self.__dma_alligned = None
        self.__allign() ### ---> init,final
        self.rot,self.rms = RotationMatrix(initial=self.initial,final=self.final)
        if abs(self.rms)>0.0001: print " Warning! Not orthogonal set! (rms=%f)"%self.rms

        self.xyz=numpy.dot(self.xyz,self.rot)  # rotate
        self.trans = self.xyz.copy()[self.atid[0]-1]
        self.xyz-=self.trans                   # translate
        if dma is not None: self.__dma_alligned = self.allignDMA(dma,self.rot,-self.trans); print " DMA is alligned!\n"

    def allignDMA(self,dma,rot,trans):
        dma_copy=dma.copy()
        dma_copy.MAKE_FULL()
        dma_copy.Rotate(rot)
        dma_copy.translate(trans)
        return dma_copy
        
    def get_transformed(self):
        return self.__dma_alligned, self.xyz
    
    def __allign(self):
        axes=numpy.identity(3,dtype=numpy.float64)[[self.axes]]
        if self.vec:
           init = numpy.zeros((3,3),dtype=numpy.float64)
           for c in [0,1,2]: 
               init[self.axes[c]] = self.vec[c]
               init[self.axes[c]]/=numpy.sqrt(numpy.sum(init[self.axes[c]]**2))
        
           self.initial=init
        elif self.atid:
           if self.atid[2]<1:
              P1 = self.xyz[self.atid[0]-1]
              P2 = self.xyz[self.atid[1]-1]
              X,Y,Z = P2-P1
              y=z=1.0; x = - (Y+Z)/X
              P3 = numpy.array([x,y,z]); P3/= numpy.linalg.norm(P3)
           else:
             P1 = self.xyz[self.atid[0]-1]
             P2 = self.xyz[self.atid[1]-1]
             P3 = self.xyz[self.atid[2]-1]
           C = P2 - P1
           B = numpy.cross(C,P3 - P1)
           A = numpy.cross(B,C)

           self.initial=numpy.array([A,B,C])
           for c in [0,1,2]: 
               self.initial[c]/=numpy.sqrt(sum(self.initial[c]**2))
        self.final=axes          

class ModifyStruct(object):
    """structure modifier"""
    def __init__(self,xyz):
        self.xyz = xyz
        self.ring = numpy.zeros((1,3),dtype=numpy.float64)
        self.n_atoms = len(xyz)
        self.rings = []
        
    def write(self,name,units='angs'):
        ring = self.ring.copy()
        if units=='angs': ring*= units.UNITS.BohrToAngstrom
        out = open(name,'w')
        out.write('%d\n\n' % len(ring))
        for i in range(len(ring)):
            out.write(" X %13.6f %13.6f %13.6f\n"%tuple(ring[i]))
        out.write('\n')
        return
        
    def makeRing(self,p1,p2,p3,n,r,scale=0):
        """kreuje obwolutek wokół atomu p1 składający się z n punktów
        oddalonych od atomu p1 o odległość r. p2 - punkt określający oś 'z'
        p3 - punkt określający oś horyzontalną (chyba 'x'). Numery atomów
        są normalne (zaczynąją się od 1)."""
        new, center, rot = self.__makeAxes(p1-1,p2-1,p3-1,scale)
        obw = numpy.zeros((n,3),dtype=numpy.float64)
        for i in range(n):
            obw[i,0] = r  * numpy.cos(2*math.pi*i/n)
            obw[i,1] = r  * numpy.sin(2*math.pi*i/n)
        obw = numpy.dot(obw,rot) + numpy.array([center])
        self.ring = numpy.concatenate((self.ring,obw),axis=0)
        self.rings.append(obw)
        return 
    
    def makeMidBonds(self,all=True,bonds=None):
        """adds dummy atom in the middle between atoms. Bonds is a list of (i,j) elements
where i and j is the atom ID (starting from 1)."""
        midBonds = []
        if all:
           for i in range(self.n_atoms):
               for j in range(i):
                   point = 0.5 * (self.xyz[i]+self.xyz[j])
                   midBonds.append(point)
        else:
            for i in bonds:
                point = 0.5 * (self.xyz[i[0]-1]+self.xyz[i[1]-1])
                midBonds.append(point)

        midBonds = numpy.array( midBonds, dtype=numpy.float64)
        #
        self.ring = numpy.concatenate((self.ring, midBonds),axis=0)
        return
    
    def add(self,xyz):
        """add points to ring"""
        pass
    
    def shrink(self,threshold=0.5):
        """find groups of points lying to close to one another
        and delete them"""
        n_points = len(self.ring[1:])
        A_ij = DistanceRelationMatrix(self.ring[1:],threshold=threshold)
        g = GROUPS(A_ij).groups
        n_groups = len(g)
        
        if n_points==n_groups:
           print "\n No groups found for thershold = %.5f a.u.\n" % threshold
        else:
           n_del = n_points - n_groups
           ring = [numpy.zeros(3,dtype=numpy.float64)]
           for group in g:
               average_point = numpy.zeros(3,dtype=numpy.float64)
               for i in group:
                   average_point+=self.ring[i+1]
               average_point /= len(group)
               ring.append(average_point)
               #ring.append(self.ring[group[0]+1])
           self.ring = numpy.array( ring, dtype=numpy.float64)
           print "\n %i points deleted for thershold = %.5f a.u.\n" % (n_del,threshold)
           
        return
    
    def reset(self):
        """resets previous changes"""
        self.ring = numpy.zeros((1,3),dtype=numpy.float64)
        self.rings = []
        return
    
    def __makeAxes(self,p1,p2,p3,scale=0):
        """create the local axes with the center at p1, z-axis goint through p2 and (probably x)-axis throug p3"""
        # case for determine P3 point automatically
        if p3<0:
           P1,P2 = self.xyz[(p1,p2),]
           D = (P2-P1)/numpy.linalg.norm(P2-P1)
           X,Y,Z = D[:]
           if   X > 0.00001:
              y=1.0; z=0.0; x = -Y/X
           elif Y > 0.00001:
              z=1.0; x=0.0; y = -Z/Y
           else:
              x=1.0; y=0.0; z = -X/Z
           P3 = numpy.array([x,y,z]); P3/= numpy.linalg.norm(P3)
        # case, where P3 is provided
        else:
            P1,P2,P3 = self.xyz[(p1,p2,p3),]

        # determine the local axes
        c = P2 - P1
        c/= numpy.linalg.norm(c)
        b = numpy.cross(c,P3-P1)
        b/= numpy.linalg.norm(b)
        a = numpy.cross(b,c)
        a/= numpy.linalg.norm(a)
        
        old = numpy.identity(3,dtype=numpy.float64)
        new = numpy.array([a,b,c],dtype=numpy.float64)
        
        rot, rms = RotationMatrix(initial=old,final=new)

        return new, P1 + (P2 - P1) * scale , rot

    def __repr__(self):
        """print the status"""
        log = '\n'
        log+= ' Number of points: %10i\n' % (len(self.ring)-1)
        #log+= '\n'
        return str(log)
    
class status(object):
    """defines the status of the object"""
    def __init__(self,some_object):
        self.object = some_object
        # is object rotated or in its initial orientation?
        self.__rotated = False
        # is object translated or in its initial position?
        self.__translated = False
        # is object a copy of another object?
        self.__copied = False
    
    def get_object(self):
        return self.object
    
    def get_status(self): 
        return self.__rotated, self.__translated, self.__copied
    
    def set_status(rotated=None,copied=None,translated=None):
        if rotated     is not None: self.__rotated     = rotated
        if copied      is not None: self.__copied      = copied
        if translated  is not None: self.__translated  = translated
        
class ROTATE:
    """Rotates the object from gasphase orientation to target orientation in    
    solvent. The form or rotation depends of the type of object:                
    - list: is treated as eigenvector matrix of dimension 3NxM where N is number
      of atoms and M is number of modes. The shape of that matrix is the same   
      as from FREQ class;                                                       
    - ndarray: is treated as a list of DMA objects (e.g. first DMA derivatives);
    - DMA: is one DMA object obviously.                                         
    Final is the solute target structure and initial is gas phase structure.    
                                                                                
    The class is a container for storing rotated eigenvectors and DMA objects   
    in one place."""
    
    def __init__(self,initial=0,final=0,object=None):
        self.__initial=numpy.array(initial)
        self.__final=numpy.array(final)
        self.__objects = []
        if object is not None: self.__objects.append(status(object))
        ### superimpose structures
        self.__superimpose()

    def __superimpose(self):
        """compute rotation matrix and translation vector"""
        sup = SVDSuperimposer()
        sup.set(self.__final,self.__initial)
        sup.run()
        rms = sup.get_rms()
        rot, transl = sup.get_rotran()
        transformed_gas_phase_str = sup.get_transformed()
            
        self.__rot = rot
        self.__transformed_gas_phase_str = transformed_gas_phase_str 
        self.__transl = transl
        return
    
    def get(self):
        """return rotated objects in a list"""
        A = []
        for i in self.__objects:
            A.append(i.object)
        return A
    
    def rotate(self):         
        """rotate the object if is not rotated"""
        for object in self.__objects:
            if not object.get_status()[0]:
              
              if object.object.__class__.__name__ == 'list':
              ### rotate the DMA list
                for dmai in object.object:
                    dmai.pos  =numpy.array(self.__initial)
                    dmai.origin  = numpy.array(self.__initial)
                    dmai.MAKE_FULL()
                    dmai.Rotate(self.__rot)
              elif object.object.__class__.__name__ == 'ndarray':
              ### rotate the eigenvectors
                  N,M = object.object.shape; N/=3
                  object.object = object.object.reshape(N,3,M)
                  object.object = numpy.tensordot(object.object,self.__rot,(1,0))   # dimension: nstat,nmodes,3
                  object.object = numpy.transpose(object.object,(0,2,1))            # dimension: nstat,3,nmodes
                  object.object = object.object.reshape(N*3,M)                # dimension: nstat*3,nmodes
              elif object.object.__class__.__name__ == 'DMA':
              ### rotate the DMA object 
                  object.object.MAKE_FULL()
                  object.object.Rotate(self.__rot)
                  
              ### set status to rotated
              object.set_status(rotated=True)
if import_matplotlib:    
   class Grid2D:
       """represents 2D-grid of points"""
       def __init__(self,
                    xmin=0, xmax=1, dx=0.5,
                    ymin=0, ymax=1, dy=0.5,):
           # coordinates in each space direction
           nx = numpy.int64((xmax-xmin)/dx + 1)
           ny = numpy.int64((ymax-ymin)/dy + 1)
           
           x,y = numpy.mgrid[0:nx,0:ny]
           
           # store for convenience
           self.dx = dx; self.dy = dy
           self.nx = nx; self.ny = ny
           self.shape = (self.nx,self.ny)
           
           # make 3D versions of the coordinate arrays
           # (needed for vectorized  function evaluators)
           self.xcoorv = numpy.float64(x)*dx + xmin
           self.ycoorv = numpy.float64(y)*dy + ymin
   
       def __init________(self,
                    xmin=0, xmax=1, dx=0.5,
                    ymin=0, ymax=1, dy=0.5):
           # coordinates in each space direction
           self.xcoor = scitools.numpyutils.seq(xmin, xmax, dx)
           self.ycoor = scitools.numpyutils.seq(ymin, ymax, dy)
           
           # store for convenience
           self.dx = dx;  self.dy = dy
           self.nx = self.xcoor.size;  self.ny = self.ycoor.size
           self.shape = (self.nx,self.ny)
           # make 2D versions of the coordinate arrays
           # (needed for vectorized  function evaluators)
           ###self.xcoorv = self.xcoor[:, numpy.newaxis]
           ###self.ycoorv = self.ycoor[numpy.newaxis, :]
           self.ycoorv, self.xcoorv = numpy.meshgrid(self.ycoor,self.xcoor)
       
       def eval(self,f,**kwargs):
           """Evaluate vectorized function f at each grid point"""
           return f(self.xcoorv,self.ycoorv,**kwargs)
   
       def newaxis(self,axis,val):
           """change values of axis to new values"""
           # x-axis
           if   axis==0:
                self.xcoorv.fill(1.)
                self.xcoorv*= val[:,numpy.newaxis]
           # y-axis
           elif axis==1:
                self.ycoorv.fill(1.)
                self.ycoorv*= val[numpy.newaxis,:]
           else: raise IndexError
           return
    
class Grid3D:
    """represents 3D-grid of points"""
    def __init__(self,
                 xmin=0, xmax=1, dx=0.5,
                 ymin=0, ymax=1, dy=0.5,
                 zmin=0, zmax=1, dz=0.5):
        # coordinates in each space direction
        nx = numpy.int64((xmax-xmin)/dx + 1)
        ny = numpy.int64((ymax-ymin)/dy + 1)
        nz = numpy.int64((zmax-zmin)/dz + 1)
        
        x,y,z = numpy.mgrid[0:nx,0:ny,0:nz]
        
        # store for convenience
        self.dx = dx; self.dy = dy; self.dz = dz
        self.nx = nx; self.ny = ny; self.nz = nz
        self.shape = (self.nx,self.ny,self.nz)
        
        # make 3D versions of the coordinate arrays
        # (needed for vectorized  function evaluators)
        self.xcoorv = numpy.float64(x)*dx + xmin
        self.ycoorv = numpy.float64(y)*dy + ymin
        self.zcoorv = numpy.float64(z)*dz + zmin
            
    def eval(self,f,**kwargs):
        """Evaluate vectorized function f at each grid point"""
        return f(self.xcoorv,self.ycoorv,self.zcoorv,**kwargs)
    
    def newaxis(self,axis,val):
        """change values of axis to new values"""
        # x-axis
        if   axis==0:
             self.xcoorv.fill(1.)
             self.xcoorv*= val[:,numpy.newaxis,numpy.newaxis]
        # y-axis
        elif axis==1:
             self.ycoorv.fill(1.)
             self.ycoorv*= val[numpy.newaxis,:,numpy.newaxis]
        # z-axis
        elif axis==2:
             self.zcoorv.fill(1.)
             self.zcoorv*= val[numpy.newaxis,numpy.newaxis,:]
        else: raise IndexError
        return

def func_to_method(func, class_, method_name=None): 
    """inserts a method to a given class!:"""
    setattr(class_, method_name or func.__name__, func)
    
class Integrator:
    def __init__(self):
        self.setup()
    
    def setup(self):
        self.weights  = None
        self.points   = None
        
    def eval(self,f):
        sum = 0.0
        for i in xrange(len(self.points)): sum += self.weights[i]*f(self.points[i])
        return sum
    
class Trapezoidal(Integrator):
    def setup(self):
        self.weights = (1,1)
        self.points  = (-1,1)
        
class Simpson(Integrator):
    def setup(self):
        self.weights = (1./3.,4./3.,1./3.)
        self.points  = (-1,0,1)
        
class GaussLegendre2(Integrator):
    def setup(self):
        p = 1./math.sqrt(3)
        self.weights = (-p,p)
        self.points  = (1,1)
        
class TransFunc:
    def __init__(self,f,h,a):
        self.f=f; self.h =h; self.a=a
    def coor_mapping(self,xi):
        return self.a + (self.j - 0.5)*self.h + 0.5*self.h*xi
    def __call__(self,xi):
        x = self.coor_mapping(xi)
        return self.f(x)
        
def integrate(integrator, a, b, f, n):
    sum = 0.0
    h = (b-a)/float64(n)
    g = TransFunc(f,h,a)
    for j in xrange(1,n+1):
        g.j = j
        sum += integrator.eval(g)
    return 0.5*sum*h

def Newton( n, k ):
    """Newton symbol n nad k huhaha"""
    Wynik=1
    for i in range( 1, k+1 ):
        Wynik = Wynik * ( n - i + 1 ) / i
    return Wynik



def PRINT(vec):
    """ print helper 4 """
    log=''
    for i in range(len(vec)):
        log+= " %3d : %12.6f " % ( i+1 , vec[i] )
        if (i+1)%6==0: log+= "\n"
        else: continue
    print log
    print



def PRINTV(M,list1,list2,list3):
    """ print helper 3 """
    d = 6
    L = len(numpy.transpose(M))

    if L % d == 0:
       n = L / d
    else:
       n = L / d + 1

    if list1 == '' and list2 == '':
       for b in range(n):
           try:
               m = M[:,(b*d):((b+1)*d)]
           except IndexError:
               m = M[:,(b*d):-1]

           for u in range(len(m)):
             for i in range(len(numpy.transpose(m))):
               v = "%.6f" % m[u][i]
               print "%10s" % v.rjust(10),
             print
           print
    elif list2 == '':
        for b in range(n):
           try:
               l1 = list1[(b*d):((b+1)*d)]
               m  =   M[:,(b*d):((b+1)*d)]
           except IndexError:
               l1 = list1[(b*d):-1]
               m =    M[:,(b*d):-1]

           for i in range(len(l1)):
               t1 = "%s" % l1[i]
               print "%15s" % t1.rjust(15),
           print
           for i in range(len(l1)):
               kk = '-'*13
               print "%s" % kk.rjust(15),
           print

           for u in range(len(m)):
             for i in range(len(numpy.transpose(m))):
               v = "%.3f" % m[u][i]
               print "%15s" % v.rjust(15),
             t3 = "%s" % list3[u]
             print ': %4s' % t3.rjust(4)
           print
    else:

       for b in range(n):
           try:
               l1 = list1[(b*d):((b+1)*d)]
               l2 = list2[(b*d):((b+1)*d)]
               m  =   M[:,(b*d):((b+1)*d)]
           except IndexError:
               l1 = list1[(b*d):-1]
               l2 = list2[(b*d):-1]
               m =    M[:,(b*d):-1]

           for i in range(len(l1)):
               t1 = "%4d" % round(l1[i],0) # oryginalna wersja: "%s" % l1[i]
               print "%15s" % t1.rjust(15),
           print
           for i in range(len(l1)):
               t2 = "%s" % l2[i]
               print "%15s" % t2.rjust(15),
           print
           for i in range(len(l1)):
               kk = '-'*13
               print "%s" % kk.rjust(15),
           print
           
           for u in range(len(m)):
             for i in range(len(numpy.transpose(m))):
               v = "%.3f" % m[u][i] # oryginalna wersja: %.5E
               print "%15s" % v.rjust(15),
             t3 = "%4d" % round(list3[u],0)   # oryginalna wersja: "%s" % list3[u]
             print ': %4s' % t3.rjust(4)
           print
           
def PUPA(a):
    """ print helper 1 """
    log=""
    for i in range(len(a)):
        for j in range(len(a[0])):
            log+= "%5.2f " % a[i][j]
            if j==len(a[0])-1: log+='\n'
    return log

def PRINTL(M,list1="",list2=""):
    """ print helper 2 """
    d = 5
    L = len(numpy.transpose(M))

    if L % d == 0:
       n = L / d
    else:
       n = L / d + 1

    if list1 == '' and list2 == '':
       for b in range(n):
           try:
               m = M[:,(b*d):((b+1)*d)]
           except IndexError:
               m = M[:,(b*d):-1]

           for u in range(len(m)):
             for i in range(len(numpy.transpose(m))):
               v = "%12.6f" % m[u][i]
               print "%14s" % v.rjust(14),
             print
           print

    elif list2 == '':

       for b in range(n):
           try:
               l1 = list1[(b*d):((b+1)*d)]
               m  = M[:,(b*d):((b+1)*d)]
           except IndexError:
               l1 = list1[(b*d):-1]
               m = M[:,(b*d):-1]

           for i in range(len(l1)):
               t = "%s" % l1[i]
               print "%10s" % t.rjust(10),
           print
           for i in range(len(l1)):
               kk = '-'*8
               print "%s" % kk.rjust(10),
           print

           for u in range(len(m)):
             for i in range(len(numpy.transpose(m))):
               v = "%.6e" % m[u][i]
               print "%10s" % v.rjust(10),
             print
           print

    else:

       for b in range(n):
           try:
               l1 = list1[(b*d):((b+1)*d)]
               l2 = list2[(b*d):((b+1)*d)]
               m  =   M[:,(b*d):((b+1)*d)]
           except IndexError:
               l1 = list1[(b*d):-1]
               l2 = list2[(b*d):-1]
               m =    M[:,(b*d):-1]

           for i in range(len(l1)):
               t1 = "%s" % l1[i]
               print "%15s" % t1.rjust(15),
           print
           for i in range(len(l1)):
               t2 = "%s" % l2[i]
               print "%15s" % t2.rjust(15),
           print
           for i in range(len(l1)):
               kk = '-'*13
               print "%s" % kk.rjust(15),
           print

           for u in range(len(m)):
             for i in range(len(numpy.transpose(m))):
               v = "%.6e" % m[u][i]
               print "%15s" % v.rjust(15),
             print
           print


def Histogram(data=[],npoints=100,out="histogram.dat"):
    """make a nice histogram"""
    
    a = min(data)
    b = max(data)
    X = numpy.linspace(a,b,npoints)
    spacing = numpy.abs(b-a)/(npoints-1)
    
    histogram = []
    for i in X:
        ni = 0
        for j,Y in enumerate(data):
            if (Y > i and Y <= (i+spacing)): ni+=1
        histogram.append(ni)
        
    result = open(out,'w')
    
    return 

def MakeSoluteAndSolventFiles(file,typ,solute_ids,charges=0):
    """create solute and solvent files"""

    # initiate target files
    solute = open( 'solute_%s'%typ , 'w')
    solvent= open( 'solvent_%s'%typ , 'w')

    # find the relevant data
    data = open(file)
    line = data.readline()
    querry = "CHELPG"
    while 1:
       if querry in line: break
       line = data.readline()
    text = []
    querry = " ---------------"
    while 1:
       if querry in line: break 
       text.append(line)
       line = data.readline()

    end = querry*2 + ' SLV' + querry*2
    text.append(end)
    text_solvent = text[:]
    text_solute = text
    # positions
    N1 = 12
    I = 0
    Imax = max(solute_ids)
    dupa = '       Atomic Center'
    while 1:
        if I>max(solute_ids):  break
        if I in solute_ids: 
           text_solvent.pop(N1)
        I+=1

    I = 0
    N1+= Imax+1
    while 1:
        if not text_solute[N1].startswith(dupa):  break
        if I not in solute_ids:
           text_solute.pop(N1)
        I+=1

    # charges
    I = 0
    querry = 'Charge='
    while 1:
       if querry in text_solvent[I]: break
       I+=1
    N1=I+2
    Sol_I = I+2
    I = 0
    while 1:
        if I>max(solute_ids):  break
        if I in solute_ids: 
           text_solvent.pop(N1)
        I+=1

    I = 0
    while 1:
       if querry in text_solute[I]: break
       I+=1
    N1=I+2; I=0
    N1+=Imax+1
    while 1:
        if text_solute[N1].startswith(' ---'): break
        if I not in solute_ids:
           text_solute.pop(N1)
        I+=1
 
    # change charges values
    if charges:
       I = Sol_I
       while 1:
          if text_solvent[I].startswith(' ---'): break
          line = text_solvent[I].split()
          charge = charges[line[1]]
          text_solvent[I] = "     %s   %s   %13.6f\n" % (line[0],line[1],charge)
          I+=1
          
    # write!
    for i in text_solvent: solvent.write(i)
    for i in text_solute: solute.write(i)
    print "\n The files: \n      < solute_%s > and < solvent_%s > \n have been saved in the actuall directory! \n" %(typ,typ)
    
    return #text_solvent, text_solute

#def bua(file_fchk,basis,bonds,vec,vec_ref,method,mult,charge):
#    """helper temporary function from solvshift.diff"""
#    molecule = Read_xyz_file(file_fchk,mol=True,
#                             mult=mult,charge=charge,
#                             name='happy dummy molecule')
#    
#    bfs        = PyQuante.Ints.getbasis(molecule,basis)
#    basis_size = len(bfs)
#    #print " - basis size= ", basis_size
#    print " - parsing %s density matrix" % method
#    dmat = ParseDmatFromFchk(file_fchk,basis_size,method)
#    def check_sim(l):
#        """check the sim list"""
#        for x,y in l:
#            i=0;j=0
#            for a,b in l:
#                if a==x: i+=1
#                if b==y: j+=1
#            if (i>1 or j>1): 
#                print " --- !ERROR! --- "
#                break
#
#    ### parse vectors and make Pipek-Mezey transformation
#    if vec is not None:
#       natoms= len(molecule.atoms)
#       SAO   = PyQuante.Ints.getS(bfs)
#       print " - ovelrap AO matrix evaluation..."
#       nae = vec
#       vec = ParseVecFromFchk(file_fchk)[:nae,:]
#       
#       print " - Pipek-Mezey localization of %i orbitals..." %nae
#       tran, vec = get_pmloca(natoms,mapi=bfs.LIST1,sao=SAO,
#                                        vecin=vec,nae=nae,
#                                        maxit=100000,conv=1.0E-10,
#                                        lprint=False,
#                                        freeze=None)
#       vec, sim = order(vec_ref,vec,start=0)
#       print sim
#       check_sim(sim)
#    ### calculate CAMMs
#    print " - multipole integrals in AO basis evaluation..."
#    camm = coulomb.multip.MULTIP(molecule=molecule,
#                                 basis=basis,
#                                 method='b3lyp',
#                                 matrix=dmat,
#                                 transition=False,
#                                 bonds=bonds,vec=vec)
#    print " - calculation of %s"%camm.operation
#    camm.camms()
#    #camm.mmms()
#    #camm.__printMMMs__()
#    #CAMM.__printCAMMs__()
#    
#    dma = camm.get()[0]
#    dma.write(file_fchk[:-5]+'.camm')
#    print " --- Writing file:  :", file_fchk[:-5]+'.camm'    
#    return
#
#def gen_camm(file=None, basis='6-311++G**',bonds=[],ncpus=4,
#             vec=None,vec_ref=None,natoms=7,method='SCF',mult=1,charge=0,
#             fchk_search='*_.fchk', fchk_dir=os.environ['PWD']): 
#    """
#Temporary utility:
#Calculates CAMMs/CBAMMs/LMTPs from density matrix from GAUSSIAN09.
#The DMA file with the appropriate distribution is saved.
#
#
#Usage: 
# 1) calculating CAMM: gen_camm(file=..., basis=..., method=..., fchk_search='*_.fchk', fchk_dir='./')
#
#    where 
#       o file is a FCHK file. If it is set to None (default) then it is assumed
#         you have several FCHK files in the directory specified in fchk_dir. 
#         In this case FCHK files will be searched according to fchk_search wildcard.
#         The wildcard can be set to any bash-type wildcard.
#       o basis is a string defining basis set name compatible with PyQuante
#       o method refers to density level from FCHK (SCF, MP2, CC etc)
#
# 2) calculating CBAMM: the same as above but add bonds (see in the code of this function)
# 3) calculating LMTP: the same as above but specify vec (LCAO-MO coefficients for canonical set of orbitals)
#    and natoms (number of atoms). 
#
#Wanring: Parallel run (for many fchk files) is not working at present 
#because of problem with administration of some shared libraries.
#
#Notes:
#Note that for the purpose of differentiation the option with file=None is handy
#since the FCHK files are first always sorted. Then the first FCHK is assumed to be the reference
#(not perturbed) file. The option vec_ref is not implemented yet but it is trivial.
#"""
#    if vec_ref is not None: raise NotImplementedError, 'vec_ref is not supported yet! The first FCHK file (after sorting) is assumed to be the reference!'
#    if file is None:
#       pliki_fchk  = glob.glob(fchk_dir+'/'+fchk_search)
#       pliki_fchk.sort()
#       print "\n Kolejność plików. Sprawdź czy się zgadzają!\n"  
#       for i in range(len(pliki_fchk)):
#           print pliki_fchk[i]
#       print
#    else:
#       pliki_fchk = [file]
#    
#    # compute reference vectors
#    if vec is not None:
#       ref_mol = Read_xyz_file(pliki_fchk[0],mol=True,
#                                         mult=mult,charge=charge,
#                                         name='happy dummy molecule')
#       bfs_ref    = PyQuante.Ints.getbasis(ref_mol,basis)
#       basis_size = len(bfs_ref)
#       sao_ref    = PyQuante.Ints.getS(bfs_ref)
#       print " - basis size= ", basis_size
#       nae = vec
#       print " - nae       = ", len(vec)
#       vec_ref = ParseVecFromFchk(pliki_fchk[0])[:nae,:]
#       t, vec_ref = get_pmloca(natoms,mapi=bfs_ref.LIST1,
#                               sao=sao_ref,
#                               vecin=vec_ref,nae=nae,
#                               maxit=100000,conv=1.0E-19,
#                               lprint=False,
#                               freeze=None)
#    # submit the jobs!
#    for file_fchk in pliki_fchk:
#        bua (file_fchk,basis,bonds,vec,vec_ref,method,mult,charge)
#    print
#    return
#
#
if __name__ == '__main__':
   from sys import argv
   a=ParseDMA(argv[1],argv[2])#[0]
   #b=ParseDMA(argv[1][:-4]+'log','gaussian')[0]
   #a.pos = array(b.pos)
   #a.origin = array(b.pos)
   print a
   print a.OverallMoments_old()
