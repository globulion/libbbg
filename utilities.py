# ------------------------------------------ #
#                UTILITIES                   #
# ------------------------------------------ #

__all__=['SVDSuperimposer','ParseDMA','RotationMatrix',
         'Periodic','Atomn','DMAMadrixMultiply',
         'PRINT','PRINTV','PRINTL','ParseDmatFromFchk',
         'Emtp','FrequencyShift','Grid2D','func_to_method',
         'Trapezoidal','Simpson','GaussLegendre2','integrate',
         'Read_xyz_file','Vr_dma','Allign','ElectricField',
         'FrequencyShiftPol','Newton','Parse_EDS_InteractionEnergies',
         'CalcStep','ModifyStruct','ParseUnitedAtoms',
         'MakeSoluteAndSolventFiles','GROUPS','DistanceRelationMatrix',
         'status','ROTATE','get_tcf','choose','get_pmloca',
         'ParseVecFromFchk','interchange','Peak']

import re, gentcf, orbloc, PyQuante, clemtp, scipy.optimize
from numpy import transpose, zeros, dot, \
                  float64, shape, array, \
                  sqrt, ceil, tensordot, \
                  cross, sum, where    , \
                  concatenate, average , \
                  exp
from numpy.linalg import svd, det, norm
from dma   import DMA
from units import *
from re_templates import *
import copy, os, math
#if bool(os.environ.get('__IMPORT_EASYVIZ__')):
from scitools.all import *

class Peak:
    """
Represent a peak in the spectrum of signals

Usage:

a = Peak(x=<coordsx>,y=<coordsy>)
a.set_peak(<n>=1,<func_name>='g')
a.fit([options])
a.get_r2()
a.get_param()

Notes:

 1) <x> and <y> are lists or numpy arrays
 2) <n> is number of subpeaks subject to extract
 3) 'g' - Gaussian distributions (default)
 4) [options] is the list of lists or tuples
    of the form:

      [['par1',val1],
       ['par2',val2],
           . . .
       ['parn',valn]]

    Parameter (par) have to be strings 
    whereas initial values (valn) floats.

    [options] depend on a type of function
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
        self.__func = 'g'

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
           
        self.n = n
    
    def get_parameters(self):
        print " No fitting performed - no parameters."
        return self.param
    
    def get_fit(self):
        if self.param is None: 
           print " No fitting performed - no parameters."
           return None
        else:
           return self.func(**self.args)

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
            peaks.append(peak)
        return array(peaks,dtype=float64)

    def get_r2(self):
        """return R^2 coefficient of fitting"""
        data_av = average(self.y)
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
        return sum((self.y - self.func(**self.args))**2.)

    ### pure Gaussian profiles
        
    def _gauss1(self,xo_1,sigma_1,A_1):
        """single Gaussian distribution"""
        return (A_1/(sigma_1*math.sqrt(2*math.pi)))\
               * exp(-(self.x-xo_1)**2/(2*sigma_1**2))

    def _gauss2(self,xo_1,sigma_1,A_1,
                     xo_2,sigma_2,A_2):
        """bimodal gaussian distribution"""
        g1 = (1./(sigma_1*math.sqrt(2*math.pi)))\
               * exp(-(self.x-xo_1)**2/(2*sigma_1**2))
        g2 = (1./(sigma_2*math.sqrt(2*math.pi)))\
               * exp(-(self.x-xo_2)**2/(2*sigma_2**2))
        return A_1*g1 + A_2*g2

    def _gauss3(self,xo_1,sigma_1,A_1,
                     xo_2,sigma_2,A_2,
                     xo_3,sigma_3,A_3):
        """trimodal gaussian distribution"""
        g1 = (1./(sigma_1*math.sqrt(2*math.pi)))\
               * exp(-(self.x-xo_1)**2/(2*sigma_1**2))
        g2 = (1./(sigma_2*math.sqrt(2*math.pi)))\
               * exp(-(self.x-xo_2)**2/(2*sigma_2**2))
        g3 = (1./(sigma_3*math.sqrt(2*math.pi)))\
               * exp(-(self.x-xo_3)**2/(2*sigma_3**2))
        return A_1*g1 + A_2*g2 + A_3*g3

    def _gauss4(self,xo_1,sigma_1,A_1,
                     xo_2,sigma_2,A_2,
                     xo_3,sigma_3,A_3,
                     xo_4,sigma_4,A_4):
        """trimodal gaussian distribution"""
        g1 = (1./(sigma_1*math.sqrt(2*math.pi)))\
               * exp(-(self.x-xo_1)**2/(2*sigma_1**2))
        g2 = (1./(sigma_2*math.sqrt(2*math.pi)))\
               * exp(-(self.x-xo_2)**2/(2*sigma_2**2))
        g3 = (1./(sigma_3*math.sqrt(2*math.pi)))\
               * exp(-(self.x-xo_3)**2/(2*sigma_3**2))
        g4 = (1./(sigma_4*math.sqrt(2*math.pi)))\
               * exp(-(self.x-xo_4)**2/(2*sigma_4**2))

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
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigma_1**2.)*(self.x-xo_1)**2.)/sigma_1
        return A_1 * lg1
    
    def _lg12(self,xo_1,sigma_1,A_1,m_1,
                   xo_2,sigma_2,A_2,m_2):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigma_1/(4.*(self.x-xo_1)**2. + sigma_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigma_1**2.)*(self.x-xo_1)**2.)/sigma_1
              
        lg2 = m_2*2./math.pi * sigma_2/(4.*(self.x-xo_2)**2. + sigma_2**2.) + \
              (1.-m_2)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigma_2**2.)*(self.x-xo_2)**2.)/sigma_2
        return A_1 * lg1 + A_2 * lg2

    def _lg13(self,xo_1,sigma_1,A_1,m_1,
                   xo_2,sigma_2,A_2,m_2,
                   xo_3,sigma_3,A_3,m_3,):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigma_1/(4.*(self.x-xo_1)**2. + sigma_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigma_1**2.)*(self.x-xo_1)**2.)/sigma_1
              
        lg2 = m_2*2./math.pi * sigma_2/(4.*(self.x-xo_2)**2. + sigma_2**2.) + \
              (1.-m_2)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigma_2**2.)*(self.x-xo_2)**2.)/sigma_2

        lg3 = m_3*2./math.pi * sigma_3/(4.*(self.x-xo_3)**2. + sigma_3**2.) + \
              (1.-m_3)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigma_3**2.)*(self.x-xo_3)**2.)/sigma_3
        return A_1 * lg1 + A_2 * lg2 + A_3 * lg3
    
    def _lg14(self,xo_1,sigma_1,A_1,m_1,
                   xo_2,sigma_2,A_2,m_2,
                   xo_3,sigma_3,A_3,m_3,
                   xo_4,sigma_4,A_4,m_4,):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigma_1/(4.*(self.x-xo_1)**2. + sigma_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigma_1**2.)*(self.x-xo_1)**2.)/sigma_1
              
        lg2 = m_2*2./math.pi * sigma_2/(4.*(self.x-xo_2)**2. + sigma_2**2.) + \
              (1.-m_2)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigma_2**2.)*(self.x-xo_2)**2.)/sigma_2

        lg3 = m_3*2./math.pi * sigma_3/(4.*(self.x-xo_3)**2. + sigma_3**2.) + \
              (1.-m_3)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigma_3**2.)*(self.x-xo_3)**2.)/sigma_3
              
        lg4 = m_4*2./math.pi * sigma_4/(4.*(self.x-xo_4)**2. + sigma_4**2.) + \
              (1.-m_4)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigma_4**2.)*(self.x-xo_4)**2.)/sigma_4
        return A_1 * lg1 + A_2 * lg2 + A_3 * lg3 + A_4 * lg4

    ### pseudo-Voigt-2 profiles
    
    def _lg21(self,xo_1,sigmaL_1,sigmaG_1,A_1,m_1):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigmaL_1/(4.*(self.x-xo_1)**2. + sigmaL_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigmaG_1**2.)*(self.x-xo_1)**2.)/sigmaG_1
        return A_1 * lg1

    def _lg22(self,xo_1,sigmaL_1,sigmaG_1,A_1,m_1,
                   xo_2,sigmaL_2,sigmaG_2,A_2,m_2):
        """single Lorenzian distribution"""
        lg1 = m_1*2./math.pi * sigmaL_1/(4.*(self.x-xo_1)**2. + sigmaL_1**2.) + \
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigmaG_1**2.)*(self.x-xo_1)**2.)/sigmaG_1
              
        lg2 = m_2*2./math.pi * sigmaL_2/(4.*(self.x-xo_2)**2. + sigmaL_2**2.) + \
              (1.-m_2)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigmaG_2**2.)*(self.x-xo_2)**2.)/sigmaG_2
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
              (1.-m_1)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigmaG_1**2.)*(self.x-xo_1)**2.)/sigmaG_1
              
        lg2 = m_2*2./math.pi * sigmaL_2/(4.*(self.x-xo_2)**2. + sigmaL_2**2.) + \
              (1.-m_2)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigmaG_2**2.)*(self.x-xo_2)**2.)/sigmaG_2

        lg3 = m_3*2./math.pi * sigmaL_3/(4.*(self.x-xo_3)**2. + sigmaL_3**2.) + \
              (1.-m_3)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigmaG_3**2.)*(self.x-xo_3)**2.)/sigmaG_3
              
        lg4 = m_4*2./math.pi * sigmaL_4/(4.*(self.x-xo_4)**2. + sigmaL_4**2.) + \
              (1.-m_4)*math.sqrt(4.*math.log(2.)/math.pi)*exp( (-4.*math.log(2.)/sigmaG_4**2.)*(self.x-xo_4)**2.)/sigmaG_4
        return A_1 * lg1 + A_2 * lg2 + A_3 * lg3 + A_4 * lg4
                                              
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
    ind = array(ind)-1
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
    return concatenate(t)

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
    mapi = array(mapi,int) + 1
    nmos = len(vecin)
    n2   = (nmos+1)*nmos/2
    #
    tran = orbloc.pmloca(natoms=natoms,mapi=mapi,sao=sao,vecin=vecin,
                         maxit=maxit,cvgloc=conv,n2=n2,nae=nae,
                         lprint=lprint)
    #
    tran = transpose(tran)
    vecout = dot(tran,vecin)
    #
    return tran, vecout

def order(R,P,start=0):
    new_P = P.copy()
    sim   = []
    rad =  []
    for i in range(len(R)-start):
        J = 0+start
        r = 1.0E+100
        rads = []
        for j in range(len(P)-start):
            r_ = sum(( R[i+start]-P[j+start])**2)
            r__= sum((-R[i+start]-P[j+start])**2)
            if r__<r_: r_=r__
            rads.append(r_)
            if r_<r:
               r=r_
               J = j
        sim.append((i+1,J+1))
        new_P[i+start] = P[J+start]
        rad.append(rads)
    return new_P, sim#, array(rad,dtype=float)

class GROUPS:
      """ 
grouping algorithm from numerical project: 
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
    r = gentcf.gentcf(file,nmax,nskip,norgns,lprint,ndels)
    r = r[:ndels]
    # write the tcf file on disk
    if save:
       out = open(outfile,'w')
       for i in range(ndels):
           print >> out, "%10i %13.5E" % ((i+1),r[i])
       out.close()
    # return results:
    tcf = zeros((ndels,2),dtype=float64)
    tcf[:,0] = linspace(1,ndels,ndels)
    tcf[:,1] = array(r)

    return tcf

def DistanceRelationMatrix(xyz,threshold=1):
    """calculate boolean relation matrix for structure xyz (array).
    Threshold = 1 Bohr and coordinates of xyz have to be Bohr too!
    You can then search for groups using: GROUPS(A_ij).groups"""
    
    K = len(xyz)
    A_ij = zeros((K,K),dtype=bool)
    for i in range(K):
        for j in range(i):
            if sqrt(sum((xyz[i]-xyz[j])**2))<=threshold:
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
        return sqrt(sum(sum(diff*diff))/l)

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
        av1=sum(coords,axis=0)/self.n  
        av2=sum(reference_coords,axis=0)/self.n    
        coords=coords-av1
        reference_coords=reference_coords-av2
        # correlation matrix
        a=dot(transpose(coords), reference_coords)
        u, d, vt=svd(a)
        self.rot=transpose(dot(transpose(vt), transpose(u)))
        # check if we have found a reflection
        if det(self.rot)<0:
            vt[2]=-vt[2]
            self.rot=transpose(dot(transpose(vt), transpose(u)))
        self.tran=av2-dot(av1, self.rot)

    def get_transformed(self):
        "Get the transformed coordinate set."
        if self.coords is None or self.reference_coords is None:
            raise Exception("No coordinates set.")
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        if self.transformed_coords is None:
            self.transformed_coords=dot(self.coords, self.rot)+self.tran
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


def Read_xyz_file(file,ar=False,mol=False,mult=1,charge=0,name='dummy',
                  units='Angstrom'):
    """\
reads xyz or fchk file and returns coords and atoms. 
Coordinates are returned in AU!

Usage:
get_mol(file,[ar=False,mol=False,mult=1,charge=0,name='dummy',
              units='Angstrom'])

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
           coord[i][1:] = map(float64,coord[i][1:])
           if units.lower()=='angstrom':
              for j in range(3):
                  coord[i][j+1]*= UNITS.AngstromToBohr
            
       if ar:
           data = [map(float64,[x for x in coord[y][1:]]) \
                                  for y in range( len(coord))]
           data = array(data,dtype=float64)
       if mol:
           Coords = []
           for i in range(n_atoms):
               atom  = (Atom(coord[i][0]).atno, (coord[i][1], 
                                                 coord[i][2],
                                                 coord[i][3]) )
               Coords.append(atom)
           Mol = PyQuante.Molecule(name,Coords,units=units,
                                   multiplicity=mult,charge=charge)
    
       if   mol : return Mol                 
       elif ar  : return coord, data
       else     : return coord

    ### handle g09.fchk files
    elif file[-4:].lower() == 'fchk':
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
           
       atnos = array(atnos,dtype=float64)
       
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
           
       coord = array(coord,dtype=float64).reshape(n_atoms,3)

       # create Molecule object
       if mol:
           Coords = []
           for i in range(n_atoms):
               atom  = (atnos[i], (coord[i][0],
                                   coord[i][1],
                                   coord[i][2]) )
               Coords.append(atom)
           Mol = PyQuante.Molecule(name,Coords,units='Bohr',
                                   multiplicity=mult,charge=charge)
    
           return Mol                        
       else: return None
       
def Vr_dma(dma,Rb,is_full=False):
    """calculates electrostatic potential in point Rb 
from dma distribution."""

    if not is_full:
       dma.MAKE_FULL()
       dma.MakeTraceless()
    Ra,qa,Da,Qa,Oa = dma.DMA_FULL
    V=0
    for i in range(len(Ra)):
        R=Rb-Ra[i] 
        Rab=sqrt(sum(R**2,axis=0))
        V+=qa[i]/Rab
        V+=tensordot(Da[i],R,(0,0))/Rab**3 
        V+=tensordot(R,tensordot(Qa[i],R,(1,0)),(0,0))/Rab**5 
        V+=tensordot(R,tensordot(R,tensordot(Oa[i],R,(0,0)),(0,0)),(0,0))/Rab**7 
    return V

def Energy_density(dma,Rb,is_full=False):
    """calculates electrostatic potential in point Rb 
from dma distribution.
Energy_density(dma,R,full=False/True)
if full - the dma object is turned into traceless object
"""

    if not is_full:
       dma.MAKE_FULL()
       dma.MakeTraceless()
    Ra,qa,Da,Qa,Oa = dma.DMA_FULL
    e2=0
    for i in range(len(Ra)):
        R=Rb-Ra[i]
        Rab=sqrt(sum(R**2,axis=0))
        e2+=qa[i]**2/Rab**4
        #
        e2+=3*(dot(Da[i],R))**2/Rab**8
        #
        e2+=(dot(Da[i],Da[i]))**2 / Rab**6
        #
        t  =tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0))
        g  =tensordot(Qa[i],R,(0,0))
        e2+=5*t**2/Rab**12 +  4*dot(g,g)/Rab**10
        #
        t=tensordot(R,tensordot(R,tensordot(Oa[i],R,(0,0)),(0,0)),(0,0))
        e2+= 7*t**2 / Rab**16
        #
        t  =tensordot(R,tensordot(Oa[i],R,(0,0)),(0,0))
        g  =dot(t,t)
        e2+= 9*g / Rab**14
        
    return e2

def ElectricField(dma,Rb,is_full=False):
    """calculates electrostatic field in point Rb
from dma distribution. Usage:
ElectricField(dma,R,full=False/True)
if not is_full - the dma object is turned into traceless object
"""

    if not is_full:
       dma.MAKE_FULL()
       dma.MakeTraceless()
    Ra,qa,Da,Qa,Oa = dma.DMA_FULL
    field=zeros(3,dtype=float64)
    for i in range(len(Ra)):
        R=Rb-Ra[i] 
        Rab=sqrt(sum(R**2,axis=0))
        field+= qa[i] * R / Rab**3
        
        field+= 3 * R * dot(R,Da[i]) / Rab**5
        field-= Da[i]/Rab**3
        
        t  =tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0))
        field+= 5* t * R / Rab**7
        field-= 2* tensordot(Da[i],R,(0,0)) / Rab**5
        
        c=tensordot(R,tensordot(Oa[i],R,(0,0)),(0,0))
        g=tensordot(R,c,(0,0))
        field+= 7 * g * R / Rab**9
        field-= 3 * c / Rab**7
        
    return field
    
def ParseDMA(file,type='coulomb'):
    """\
============================================================================
Parse DMA from GAUSSIAN, GAMESS or COULOMB.py file. It returns a DMA object.
Usage:
ParseDMA(type='coulomb')
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
               ZerothMoments.append( float64( line.split()[2] ) )
               Structure.append( array(line.split()[-3:],dtype=float64)  )
                    
               line = data.readline()
         # ----------------------------------
         querry = " FIRST MOMENTS AT POINTS"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(3): line = data.readline()
         FirstMoments = []
         while line.split()!=[]:
               FirstMoments.append( map(float64, line.split()[1:]))
               line = data.readline()
         # ----------------------------------
         querry = " SECOND MOMENTS AT POINTS"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(3): line = data.readline()
         SecondMoments = []
         while line.split()!=[]:
               SecondMoments.append( map( float64, line.split()[1:] ))
               line = data.readline()
         # ----------------------------------
         querry = " THIRD MOMENTS AT POINTS"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(4): line = data.readline()
         ThirdMoments = []
         while 'CPU' not in line.split():
               A = map( float64, line.split()[1:] )
               line = data.readline()
               B = map( float64, line.split() )
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

         return DMA( q=array(ZerothMoments)   ,
                     m=array(FirstMoments )   ,
                     T=array(SecondMoments)   ,
                     O=array(ThirdMoments )   ,
                     pos=array(Structure)    )#,Structure

    # -----------------------------------------------------------------------------
    elif type.lower() == 'coulomb' or type.lower() == 'c':
         # return DMA object
         data = open(file)
         line = data.readline()
         querry = " Distributed zeroth-order property"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(3): line = data.readline()
         ZerothMoments = []
         Structure = []
         Origin = []
         while line.split()!=[]:
               ZerothMoments.append( float64( line.split()[0] ) )
                    
               line = data.readline()
         # ----------------------------------
         querry = " Distributed first-order property"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(4): line = data.readline()
         FirstMoments = []
         while line.split()!=[]:
               FirstMoments.append( map(float64, line.split()[:]))
               line = data.readline()
         # ----------------------------------
         querry = " Distributed second-order property"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(4): line = data.readline()
         SecondMoments = []
         while line.split()!=[]:
               SecondMoments.append( map( float64, line.split()[:] ))
               line = data.readline()
         # ----------------------------------
         querry = " Distributed third-order property"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(5): line = data.readline()
         ThirdMoments = []
         while '-----' not in line:
               A = map( float64, line.split()[:] )
               line = data.readline()
               B = map( float64, line.split()[:] )
               ThirdMoments.append( A+B )
               line = data.readline()
               
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
                  atoms.append(Atom(coord[0]))
                  Structure.append( map( float64, coord[1:] ) )
                  line = data.readline()

         Structure = array(Structure,dtype=float64)                  
         
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
                  Origin.append( map( float64, coord[1:] ) )
                  line = data.readline()  
            
            Origin = array(Origin   ,dtype=float64)

         else:
            Origin    = Structure.copy()
          
         return DMA( q=array(ZerothMoments)   ,
                     m=array(FirstMoments )   ,
                     T=array(SecondMoments)   ,
                     O=array(ThirdMoments )   ,
                     atoms=atoms              ,
                     pos=Structure            ,
                     origin=Origin)#,Structure
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
               Structure.append( array(line.split()[-3:],dtype=float64)  )
               line = data.readline()

         # seek for charges!
         querry = " Charge="
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(2): line = data.readline()
         ZerothMoments = []
         for i in range(len(Structure)):
             ZerothMoments.append( float64( line.split()[-1] ) )    
             line = data.readline()
        
         Result = DMA(nfrag=len(Structure))
         Result.pos = array(Structure) * UNITS.AngstromToBohr
         Result.DMA[0] = array(ZerothMoments)
         
         return Result#, array(Structure) * UNITS.AngstromToBohr

def ParseVecFromFchk(file):
    """parse Ci\mu coeeficients from g09 fchk"""
    querry = "Alpha MO coefficients"
    data = open(file)
    line = data.readline()
    while querry not in line:
          line = data.readline()
    NSQ = int(line.split()[-1])
    N   = math.sqrt(NSQ)
    line = data.readline()
    C = []
    g = lambda n: n/5+bool(n%5)
    for i in range(g(NSQ)):
        C+= line.split()
        line = data.readline()
    C = array(C,dtype=float64).reshape(N,N)
    data.close()
    return C

def ParseDmatFromFchk(file,basis_size):
    """parses density matrix from Gaussian fchk file"""
        
    data = open(file)
    line = data.readline()
    querry = "Total SCF Density"
    while 1:
        if querry in line: break
        line = data.readline()
    N = int(line.split()[-1])

    line = data.readline()
    dmat = []
    for i in range(int(ceil(N/5.))): 
        dmat+=[x for x in line.split()] 
        line = data.readline()
    #dmat = array(dmat,dtype=float64)
        
    # construct explicit 2D density matrix
    P = zeros((basis_size,basis_size),dtype=float64)
    #I = 0
    for i in range(basis_size):
        for j in range(i+1):
            P[i,j] = float64(dmat.pop(0))#dmat[I]
            P[j,i] = P[i,j] #dmat[I]
            #I += 1
    data.close()
    return array(P)

def Parse_EDS_InteractionEnergies(file):
    """parses EDS interaction energies from file"""
    
    data = open(file).read()
    #line = data.readline()
    querry = r'.*VARIATIONAL-PERTURBATIONAL DECOMPOSITION.*'
    E = ['DE\\(HL\\)'   , 'E\\(EL,10\\)', 'E\\(EL,M,1\\)',
         'E\\(EL,P,1\\)', 'E\\(EX,HL\\)', 'DE\\(DEL,HF\\)', 
         'DE\\(HF\\)']
         
    for term in E:
        querry+= '\s*%s\s+(%s).*\n' % (term,re_real_e)
    querry = re.compile(querry,re.DOTALL)
    match = re.search(querry,data)
    energies = array(match.groups(),dtype=float64)    
    
    return energies

def CalcStep(step,reduced_mass):
    """return step in Angstroms when given step in normal coordinate unit [Bohr*me-1/2] 
    and reduced mass [AMU]"""
    return step * UNITS.BohrToAngstrom / sqrt(UNITS.AmuToElectronMass*reduced_mass)

def CalculateCAMM(basis='6-311++G**'): 
    """calculates CAMMs from density matrix from GAUSSIAN09
       using COULOMB.py routines. Usage:
       1) gather all the .fchk and and .log files with pop=chelpg
       2) type ./diff 
       3) the files .camm are creating!!! """
       
    from sys import argv
    import os, glob
       
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
        fragment = array(dma.pos)
        
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
            structure.append( (UNITS.atomic_numbers[frag_names[j]],
                                fragment[j]) ) 
        molecule = Molecule('mol',
                            structure,
                            multiplicity=1,
                            charge=0,
                            units='Bohr')
                            
        basis_size = len(Ints.getbasis(molecule,'6-311++G**'))
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
       
        result = DMA(nfrag=len(structure))
        result.DMA[0][:] = CAMM.Mon
        #
        result.DMA[1][:] = CAMM.Dip
        #
        result.DMA[2][:,0] = array(CAMM.Quad)[:,0,0]
        result.DMA[2][:,1] = array(CAMM.Quad)[:,1,1]
        result.DMA[2][:,2] = array(CAMM.Quad)[:,2,2]
        result.DMA[2][:,3] = array(CAMM.Quad)[:,0,1]
        result.DMA[2][:,4] = array(CAMM.Quad)[:,0,2]
        result.DMA[2][:,5] = array(CAMM.Quad)[:,1,2]
        #
        result.DMA[3][:,0] = array(CAMM.Oct)[:,0,0,0]
        result.DMA[3][:,1] = array(CAMM.Oct)[:,1,1,1]
        result.DMA[3][:,2] = array(CAMM.Oct)[:,2,2,2]
        result.DMA[3][:,3] = array(CAMM.Oct)[:,0,0,1]
        result.DMA[3][:,4] = array(CAMM.Oct)[:,0,0,2]
        result.DMA[3][:,5] = array(CAMM.Oct)[:,0,1,1]
        result.DMA[3][:,6] = array(CAMM.Oct)[:,1,1,2]
        result.DMA[3][:,7] = array(CAMM.Oct)[:,0,2,2]
        result.DMA[3][:,8] = array(CAMM.Oct)[:,1,2,2]
        result.DMA[3][:,9] = array(CAMM.Oct)[:,0,1,2]
        #
        #print result
        out = open(file_log[:-4]+'.camm','w')
        out.write(str(result))
        out.close()
        print " Writing file:  :", file_log[:-4]+'.camm'
    print

def RotationMatrix(initial=None,final=None):
    """returns rotation matrix and rms from SVD superposition of two structures.
    The initial structure is rotated into final one"""
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
    charges = zeros((N ,K  ),dtype=float64)
    dipoles = zeros((3 ,N,K),dtype=float64)
    qdrples = zeros((6 ,N,K),dtype=float64)
    octples = zeros((10,N,K),dtype=float64)
    for i,dma in enumerate(dma_list):
        charges[i,:]   = dma.DMA[0]
        dipoles[:,i,:] = transpose(dma.DMA[1])
        qdrples[:,i,:] = transpose(dma.DMA[2])
        octples[:,i,:] = transpose(dma.DMA[3])
 
    ### TRANSFORMATION!    
    charges = dot(matrix,charges)
    dipoles = dot(matrix,dipoles)
    qdrples = dot(matrix,qdrples)
    octples = dot(matrix,octples)
    
    result = []
    for i in range(len(matrix)):
        dma = DMA(nfrag=K)
        dma.DMA[0] = charges[i]
        dma.DMA[1] = transpose(dipoles[i])
        dma.DMA[2] = transpose(qdrples[i])
        dma.DMA[3] = transpose(octples[i])
        dma.pos = positions
        result.append( dma )
        
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
for i in xrange(len(Ra)):
         for j in xrange(len(Rb)):
            R    = Rb[j]-Ra[i]
            Rab=sqrt(sum(R**2,axis=0))
            if (Rab < threshold and Rab !=0):
             qq  +=   qa[i]*qb[j]/Rab                                                               # qa - qb  | R1
             #if not hash:
             qD  +=  -qa[i]*tensordot(Db[j],R,(0,0))/Rab**3                                         # qa - Db  | R2
             Dq  +=  +qb[j]*tensordot(Da[i],R,(0,0))/Rab**3                                         # qb - Da  | R2
             DD  +=-3*tensordot(Da[i],R,(0,0))*tensordot(Db[j],R,(0,0))/Rab**5                      # Da - Db  | R3
             DD  +=   tensordot(Da[i],Db[j],(0,0))/Rab**3                                           # Da - Db  | R3
             qQ  +=   qa[i]*tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))/Rab**5                      # qa - Qb  | R3
             Qq  +=   qb[j]*tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0))/Rab**5                      # qb - Qa  | R3
             DQ  +=-2*tensordot(Da[i],tensordot(Qb[j],R,(0,0)),(0,0))/Rab**5                        # Da - Qb  | R4
             QD  += 2*tensordot(Db[j],tensordot(Qa[i],R,(0,0)),(0,0))/Rab**5                        # Db - Qa  | R4
             DQ  += 5*tensordot(Da[i],R,(0,0))*tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))/Rab**7   # Da - Qb  | R4
             QD  +=-5*tensordot(Db[j],R,(0,0))*tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0))/Rab**7   # Db - Qa  | R4
             qO  +=  -qa[i]*tensordot(R,tensordot(R,tensordot(Ob[j],R,(0,0)),(0,0)),(0,0))/Rab**7   # qa - Ob  | R4
             Oq  +=   qb[j]*tensordot(R,tensordot(R,tensordot(Oa[i],R,(0,0)),(0,0)),(0,0))/Rab**7   # qb - Oa  | R4
             QQ  += (35.)/(3.)* (tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0)) *
                                 tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))  ) / Rab**9            # Qa - Qb  | R5
             OD  +=-7*(tensordot(Db[j],R,(0,0)) *
                       tensordot(R,tensordot(R,tensordot(Oa[i],R,(0,0)),(0,0)),(0,0)) ) / Rab**9    # Db - Oa  | R5
             DO  +=-7*(tensordot(Da[i],R,(0,0)) *
                       tensordot(R,tensordot(R,tensordot(Ob[j],R,(0,0)),(0,0)),(0,0)) ) / Rab**9    # Da - Ob  | R5
             QQ  +=-(20.)/(3.) * tensordot(tensordot(R,Qa[i],(0,0)),
                                           tensordot(R,Qb[j],(0,0)),(0,0)) / Rab**7                 # Qa - Qb  | R5
             QQ  +=(2.)/(3.)  * tensordot(Qa[i],Qb[j])  / Rab**5                                    # Qa - Qb  | R5
             OD  +=3 * tensordot(R,tensordot(R,tensordot(Oa[i],Db[j],(0,0)),(0,0)),(0,0)) / Rab**7  # Db - Oa  | R5
             DO  +=3 * tensordot(R,tensordot(R,tensordot(Ob[j],Da[i],(0,0)),(0,0)),(0,0)) / Rab**7  # Da - Ob  | R5
             ### The remaining terms with hexadecapoles are not implemented yet
             #Eint+= qb[j] * tensordot(R,tensordot(R,tensordot(R,tensordot(R,Ha[i],
             #                (0,0)),(0,0)),(0,0)),(0,0))   / Rab**9                                 # Ha - qb  | R5
             #Eint+= qa[i] * tensordot(R,tensordot(R,tensordot(R,tensordot(R,Hb[j],
             #                (0,0)),(0,0)),(0,0)),(0,0))   / Rab**9                                 # Hb - qj  | R5
             ### these are implemented already !
             #OQ  += 2* tensordot(tensordot(Oa[i],Qb[j],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qb - Oa  | R6
             #QO  +=-2* tensordot(tensordot(Ob[j],Qa[i],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qa - Ob  | R6
             #OQ  +=-14*tensordot(tensordot(R,tensordot(Oa[i],R,(1,0)),(0,0)) ,                      # Qb - Oa  | R6
             #                    tensordot(R,Qb[j],(0,0)) ,(0,0)) / Rab**9                          
             #QO  += 14*tensordot(tensordot(R,tensordot(Ob[j],R,(1,0)),(0,0)) ,                      # Qa - Ob  | R6
             #                    tensordot(R,Qa[i],(0,0)) ,(0,0)) / Rab**9
             #OQ  +=( 21*tensordot(tensordot(R,tensordot(Oa[i],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
             #         * tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))) / Rab**11
             #QO  +=(-21*tensordot(tensordot(R,tensordot(Ob[j],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
             #         * tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0))) / Rab**11   
             #OO  +=(2.)/(5.)*tensordot(Oa[i],Ob[j],((0,1,2),(0,1,2))) / Rab**7                      # Ob - Oa  | R7
             #OO  +=(-42./5.)*tensordot(tensordot(R,Oa[i],(0,0)),
             #                          tensordot(R,Ob[j],(0,0)),
             #                          ((0,1),(0,1))) / Rab**9                                      # Ob - Oa  | R7
             #OO  +=(189.)/(5.)*tensordot(
             #                            tensordot(tensordot(R,Oa[i],(0,0)),R,(0,0)),
             #                            tensordot(tensordot(R,Ob[j],(0,0)),R,(0,0)),
             #                            (0,0)) /Rab**11
             #OO  +=-(231./5.)*(tensordot(tensordot(tensordot(R,Oa[i],(0,0)),R,(0,0)),R,(0,0)) *
             #                  tensordot(tensordot(tensordot(R,Ob[j],(0,0)),R,(0,0)),R,(0,0)) ) /\
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
    Eint,A,B,C,D,E,CC,CD,CQ,CO,DD,DQ=clemtp.clemtp(Ra,qa,Da,Qa,Oa,Rb,qb,Db,Qb,Ob)
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
    
def get_elmtp(DMA1,DMA2):
    """calculates E(EL)MTP from two DMA distributions.
dma1 and dma2 are the objects of the class DMA. Calculations are
in atomic units and a respective interaction energy is in 
a.u. as well. """

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
    for i in xrange(len(Ra)):
         for j in xrange(len(Rb)):
             R    = Rb[j]-Ra[i]
             Rab=sqrt(sum(R**2,axis=0))
             qq  +=   qa[i]*qb[j]/Rab                                                               # qa - qb  | R1
             qD  +=  -qa[i]*tensordot(Db[j],R,(0,0))/Rab**3                                         # qa - Db  | R2
             Dq  +=  +qb[j]*tensordot(Da[i],R,(0,0))/Rab**3                                         # qb - Da  | R2
             DD  +=-3*tensordot(Da[i],R,(0,0))*tensordot(Db[j],R,(0,0))/Rab**5                      # Da - Db  | R3
             DD  +=   tensordot(Da[i],Db[j],(0,0))/Rab**3                                           # Da - Db  | R3
             qQ  +=   qa[i]*tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))/Rab**5                      # qa - Qb  | R3
             Qq  +=   qb[j]*tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0))/Rab**5                      # qb - Qa  | R3
             DQ  +=-2*tensordot(Da[i],tensordot(Qb[j],R,(0,0)),(0,0))/Rab**5                        # Da - Qb  | R4
             QD  += 2*tensordot(Db[j],tensordot(Qa[i],R,(0,0)),(0,0))/Rab**5                        # Db - Qa  | R4
             DQ  += 5*tensordot(Da[i],R,(0,0))*tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))/Rab**7   # Da - Qb  | R4
             QD  +=-5*tensordot(Db[j],R,(0,0))*tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0))/Rab**7   # Db - Qa  | R4
             qO  +=  -qa[i]*tensordot(R,tensordot(R,tensordot(Ob[j],R,(0,0)),(0,0)),(0,0))/Rab**7   # qa - Ob  | R4
             Oq  +=   qb[j]*tensordot(R,tensordot(R,tensordot(Oa[i],R,(0,0)),(0,0)),(0,0))/Rab**7   # qb - Oa  | R4
             #QQ  += (35.)/(3.)* (tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0)) *
             #                    tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))  ) / Rab**9            # Qa - Qb  | R5
             #OD  +=-7*(tensordot(Db[j],R,(0,0)) *
             #          tensordot(R,tensordot(R,tensordot(Oa[i],R,(0,0)),(0,0)),(0,0)) ) / Rab**9    # Db - Oa  | R5
             #DO  +=-7*(tensordot(Da[i],R,(0,0)) *
             #          tensordot(R,tensordot(R,tensordot(Ob[j],R,(0,0)),(0,0)),(0,0)) ) / Rab**9    # Da - Ob  | R5
             #QQ  +=-(20.)/(3.) * tensordot(tensordot(R,Qa[i],(0,0)),
             #                              tensordot(R,Qb[j],(0,0)),(0,0)) / Rab**7                 # Qa - Qb  | R5
             #QQ  +=(2.)/(3.)  * tensordot(Qa[i],Qb[j])  / Rab**5                                    # Qa - Qb  | R5
             #OD  +=3 * tensordot(R,tensordot(R,tensordot(Oa[i],Db[j],(0,0)),(0,0)),(0,0)) / Rab**7  # Db - Oa  | R5
             #DO  +=3 * tensordot(R,tensordot(R,tensordot(Ob[j],Da[i],(0,0)),(0,0)),(0,0)) / Rab**7  # Da - Ob  | R5
             ### The remaining terms with hexadecapoles are not implemented yet
             #Eint+= qb[j] * tensordot(R,tensordot(R,tensordot(R,tensordot(R,Ha[i],
             #                (0,0)),(0,0)),(0,0)),(0,0))   / Rab**9                                 # Ha - qb  | R5
             #Eint+= qa[i] * tensordot(R,tensordot(R,tensordot(R,tensordot(R,Hb[j],
             #                (0,0)),(0,0)),(0,0)),(0,0))   / Rab**9                                 # Hb - qj  | R5
             ### these are implemented already !
             #OQ  += 2* tensordot(tensordot(Oa[i],Qb[j],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qb - Oa  | R6
             #QO  +=-2* tensordot(tensordot(Ob[j],Qa[i],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qa - Ob  | R6
             #OQ  +=-14*tensordot(tensordot(R,tensordot(Oa[i],R,(1,0)),(0,0)) ,                      # Qb - Oa  | R6
             #                    tensordot(R,Qb[j],(0,0)) ,(0,0)) / Rab**9                          
             #QO  += 14*tensordot(tensordot(R,tensordot(Ob[j],R,(1,0)),(0,0)) ,                      # Qa - Ob  | R6
             #                    tensordot(R,Qa[i],(0,0)) ,(0,0)) / Rab**9
             #OQ  +=( 21*tensordot(tensordot(R,tensordot(Oa[i],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
             #         * tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))) / Rab**11
             #QO  +=(-21*tensordot(tensordot(R,tensordot(Ob[j],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
             #         * tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0))) / Rab**11   
             #OO  +=(2.)/(5.)*tensordot(Oa[i],Ob[j],((0,1,2),(0,1,2))) / Rab**7                      # Ob - Oa  | R7
             #OO  +=(-42./5.)*tensordot(tensordot(R,Oa[i],(0,0)),
             #                          tensordot(R,Ob[j],(0,0)),
             #                          ((0,1),(0,1))) / Rab**9                                      # Ob - Oa  | R7
             #OO  +=(189.)/(5.)*tensordot(
             #                            tensordot(tensordot(R,Oa[i],(0,0)),R,(0,0)),
             #                            tensordot(tensordot(R,Ob[j],(0,0)),R,(0,0)),
             #                            (0,0)) /Rab**11                                            # Oa - Ob  | R7
             #OO  +=-(231./5.)*(tensordot(tensordot(tensordot(R,Oa[i],(0,0)),R,(0,0)),R,(0,0)) *
             #                  tensordot(tensordot(tensordot(R,Ob[j],(0,0)),R,(0,0)),R,(0,0)) ) /\
             #                  Rab**13                                                              # Oa - Ob  | R7
             
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
    shift = tensordot(field,tensordot(solpol,field,(0,0)),(0,0))
    shift*= -1./2.
    
    return shift * UNITS.HartreePerHbarToCmRec

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
    result = array([A,B,C,D,E])
    # switch to cm-1
    # result  *= UNITS.HartreePerHbarToCmRec
    return result

class Allign:
    """represents alligning function"""
    def __init__(self,xyz=zeros(3),atid=[],vec=[],axes=(0,1,2),dma=0):
        self.xyz=xyz
        self.atid=atid
        self.vec=vec
        self.axes=axes
        self.__allign() ### ---> init,final
        self.rot,self.rms = RotationMatrix(initial=self.initial,final=self.final)
        if abs(self.rms)>0.0001: print " Warning! Not orthogonal set! (rms=%f)"%self.rms
        if dma: self.__dma_alligned = self.allignDMA(dma); print " DMA is alligned!\n"
        self.xyz=dot(self.xyz,self.rot)

    def allignDMA(self,dma):
        dma_copy=dma.copy()
        dma_copy.MAKE_FULL()
        dma_copy.Rotate(self.rot)
        return dma_copy
        
    def get_transformed(self):
        return self.__dma_alligned, self.xyz
    
    def __allign(self):
        axes=identity(3,dtype=float64)
        if self.vec:
           init = zeros((3,3),dtype=float64)
           for c in [0,1,2]: 
               init[self.axes[c]] = self.vec[c]
               init[self.axes[c]]/=sqrt(sum(init[self.axes[c]]**2))
        
           self.initial=init
        elif self.atid:
           P1 = self.xyz[self.atid[0]-1]
           P2 = self.xyz[self.atid[1]-1]
           P3 = self.xyz[self.atid[2]-1]
           C = P2 - P1
           B = cross(C,P3 - P1)
           A = cross(B,C)
           
           self.initial=array([A,B,C])
           for c in [0,1,2]: 
               self.initial[c]/=sqrt(sum(self.initial[c]**2))
        self.final=axes          

class ModifyStruct(object):
    """structure modifier"""
    def __init__(self,xyz):
        self.xyz = xyz
        self.ring = zeros((1,3),dtype=float64)
        self.n_atoms = len(xyz)
        self.rings = []
        
    def write(self,name,units='angs'):
        ring = self.ring.copy()
        if units=='angs': ring*= UNITS.BohrToAngstrom
        out = open(name,'w')
        out.write('%d\n\n' % len(ring))
        for i in range(len(ring)):
            out.write(" X %13.6f %13.6f %13.6f\n"%tuple(ring[i]))
        out.write('\n')
        return
        
    def makeRing(self,p1,p2,p3,n,r,scale=0):
        """kreuje obwolutek wokół atomu p1 składający się z n punktów
        oddalonych od atomu p1 o odległość r"""
        new, center, rot = self.__makeAxes(p1,p2,p3,scale)
        obw = zeros((n,3),dtype=float64)
        for i in range(n):
            obw[i,0] = r  * cos(2*pi*i/n)
            obw[i,1] = r  * sin(2*pi*i/n)
        obw = dot(obw,rot) + array([center])
        self.ring = concatenate((self.ring,obw),axis=0)
        self.rings.append(obw)
        return 
    
    def makeMidBonds(self,all=True,bonds=None):
        """adds dummy atom in the middle between atoms"""
        midBonds = []
        if all:
           for i in range(self.n_atoms):
               for j in range(i):
                   point = 0.5 * (self.xyz[i]+self.xyz[j])
                   midBonds.append(point)
        else:
            for i in bonds:
                point = 0.5 * (self.xyz[i[0]]+self.xyz[i[1]])
                midBonds.append(point)

        midBonds = array( midBonds,dtype=float64)
        #
        self.ring = concatenate((self.ring, midBonds),axis=0)
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
           ring = [zeros(3,dtype=float64)]
           for group in g:
               average_point = zeros(3,dtype=float64)
               for i in group:
                   average_point+=self.ring[i+1]
               average_point /= len(group)
               ring.append(average_point)
               #ring.append(self.ring[group[0]+1])
           self.ring = array( ring, dtype=float64)
           print "\n %i points deleted for thershold = %.5f a.u.\n" % (n_del,threshold)
           
        return
    
    def reset(self):
        """resets previous changes"""
        self.ring = zeros((1,3),dtype=float64)
        self.rings = []
        return
    
    def __makeAxes(self,p1,p2,p3,scale=0):
        P1,P2,P3 = self.xyz[(p1,p2,p3),]

        c = P2 - P1
        c/= norm(c)
        b = cross(c,P3-P1)
        b/= norm(b)
        a = cross(b,c)
        a/= norm(a)
        
        old = identity(3,dtype=float64)
        new = array([a,b,c],dtype=float64)
        
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
        self.__initial=array(initial)
        self.__final=array(final)
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
                for dma in object.object:
                    dma.pos  =array(self.__initial)
                    dma.origin  = array(self.__initial)
                    dma.MAKE_FULL()
                    dma.Rotate(self.__rot)
              elif object.object.__class__.__name__ == 'ndarray':
              ### rotate the eigenvectors
                  N,M = object.object.shape; N/=3
                  object.object = object.object.reshape(N,3,M)
                  object.object = tensordot(object.object,self.__rot,(1,0))   # dimension: nstat,nmodes,3
                  object.object = transpose(object.object,(0,2,1))            # dimension: nstat,3,nmodes
                  object.object = object.object.reshape(N*3,M)                # dimension: nstat*3,nmodes
              elif object.object.__class__.__name__ == 'DMA':
              ### rotate the DMA object 
                  object.object.MAKE_FULL()
                  object.object.Rotate(self.__rot)
                  
              ### set status to rotated
              object.set_status(rotated=True)
    
class Grid2D:
    """represents 2D-grid of points"""
    def __init__(self,
                 xmin=0, xmax=1, dx=0.5,
                 ymin=0, ymax=1, dy=0.5):
        # coordinates in each space direction
        self.xcoor = seq(xmin, xmax, dx)
        self.ycoor = seq(ymin, ymax, dy)
        
        # store for convenience
        self.dx = dx;  self.dy = dy
        self.nx = self.xcoor.size;  self.ny = self.ycoor.size
        
        # make 2D versions of the coordinate arrays
        # (needed for vectorized  function evaluators)
        self.xcoorv = self.xcoor[:, newaxis]
        self.ycoorv = self.ycoor[newaxis, :]
    
    def vectorized_eval(self,f):
        """Evaluate vectorized function f at each grid point"""
        return f(self.xcoorv,self.ycoorv)

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
    L = len(transpose(M))

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
             for i in range(len(transpose(m))):
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
             for i in range(len(transpose(m))):
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
             for i in range(len(transpose(m))):
               v = "%.5E" % m[u][i]
               print "%15s" % v.rjust(15),
             t3 = "%s" % list3[u]
             print ': %4s' % t3.rjust(4)
           print
           
def P(a):
    """ print helper 1 """
    log=""
    for i in range(len(a)):
        for j in range(len(a[0])):
            log+= "%16.8e " % a[i][j]
            if j==len(a[0])-1: log+='\n'
    return log

def PRINTL(M,list1="",list2=""):
    """ print helper 2 """
    d = 5
    L = len(transpose(M))

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
             for i in range(len(transpose(m))):
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
             for i in range(len(transpose(m))):
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
             for i in range(len(transpose(m))):
               v = "%.6e" % m[u][i]
               print "%15s" % v.rjust(15),
             print
           print


def Histogram(data=[],npoints=100,out="histogram.dat"):
    """make a nice histogram"""
    
    a = min(data)
    b = max(data)
    X = linspace(a,b,npoints)
    spacing = abs(b-a)/(npoints-1)
    
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



if __name__ == '__main__':
   from sys import argv
   a=ParseDMA(argv[1],argv[2])#[0]
   #b=ParseDMA(argv[1][:-4]+'log','gaussian')[0]
   #a.pos = array(b.pos)
   #a.origin = array(b.pos)
   print a
   print a.OverallMoments_old()
