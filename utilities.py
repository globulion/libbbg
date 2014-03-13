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
         'ParseVecFromFchk','interchange','Peak','PUPA','VIB',
         'ParseFCFromFchk','ParseDipoleDerivFromFchk',
         'ParseFockFromGamessLog','lind','order','check_sim','MakeMol',
         'ParseDistributedPolarizabilitiesFromGamessEfpFile','reorder',
         'ParseEFPInteractionEnergies','secant','RungeKutta',
         'numerov1','numerov2','simpson','simpson_nonuniform','fder5pt',
         'QMOscillator','ParseDMAFromGamessEfpFile','dihedral','Peak2DIR',]
         
__version__ = '3.3.0'

import re, gentcf, orbloc, PyQuante, clemtp, \
       scipy.optimize, scipy.integrate
from numpy import transpose, zeros, dot, \
                  float64, shape, array, \
                  sqrt, ceil, tensordot, \
                  cross, sum, where    , \
                  concatenate, average , \
                  exp, linalg, sign    , \
                  arctan2, meshgrid    , \
                  logical_and, fft     , \
                  roll, real, mgrid    , \
                  int64
from math import exp as mexp   ,\
                 sqrt as msqrt ,\
                 pi as mPi
from numpy.linalg import svd, det, norm
from dma   import DMA
from units import *
from re_templates import *
import copy, os, math
#if bool(os.environ.get('__IMPORT_EASYVIZ__')):
from scitools.all import *
from matplotlib.font_manager import FontProperties as FP
from pylab import plt, Line2D, subplots, rcParams
from scitools.numpyutils import seq
from scipy.interpolate import RectBivariateSpline as RBS, \
                              interp2d as I2D
from letters import greek as let_greek

def dihedral(A,unit='radian'):
    """Compute dihedral angle n1-n2-n3-n4. 
Usage: 
dihedral([4x3 array],<unit='radian' or 'deg'>).
Provide real atom numbers. The dihedral evaluated by this code gave opposite signs 
as compared with MOLDEN for a test NMA molecule"""
    P1 = A[1]-A[0]
    P2 = A[2]-A[1]
    P3 = A[3]-A[2]
    N1 = cross(P1,P2)
    N2 = cross(P2,P3)
    N1/= msqrt(sum(N1*N1))
    N2/= msqrt(sum(N2*N2))
    P2/= msqrt(sum(P2*P2))
    M1 = cross(N1,P2)
    x = sum(N1*N2)
    y = sum(M1*N2)
    angle = arctan2(y,x)
    conv = 180./mPi
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
    x = float64(linspace(x0,x0+(npoints-1)*step,npoints))
    y = zeros(npoints,float64)
    y[0] = y0; y[1] = y1
    #
    H = step*step/12.0
    if   isinstance(q,float64): Q = q*ones(npoints,float64)
    elif isinstance(q,ndarray): Q = q
    else:                       Q = q(x,**qarg)
    #
    if   isinstance(s,float64): S = s*ones(npoints,float64)
    elif isinstance(s,ndarray): S = s
    else:                       S = s(x,**sarg)
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
    x = float64(linspace(x0,x0+(npoints-1)*step,npoints))
    y = zeros(npoints,float64)
    y[0] = y0; y[1] = y1
    #
    H = step*step/12.0
    if   isinstance(q,float64): Q = q*ones(npoints,float64)
    elif isinstance(q,ndarray): Q = q
    else:                       Q = q(x,**qarg)
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
    while fabs(x_new-x_old)>delta:
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
   """Runge-Kutta method"""
   def __init__(self):
       self.__gfc    = None
       self.__ndim   = None
       self.__result = None
       pass

   def __call__(self):
       return self.get()

   def set(self,func,tau,init=None):
       """set an array of g-functions!"""
       self.__gfc = func
       self.__tau = tau
       self.set_init(init)
       return

   def set_init(self,init):
       """set the initial variables"""
       self.__init = array(init,float64)
       self.__ndim = len(init)
       return

   def get(self):
       return self.__result

   def eval(self,n_pass,**kwargs):
       """evaluate the dynamic properties as time goes"""
       y = zeros((n_pass,self.__ndim),float64)
       # initial conditions
       y[0] = self.__init
       #
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

class Graph:
    "graph environment"
    fs = 20
    font0 = FP()
    font1 = font0.copy()
    font1.set_family('sans-serif')
    font1.set_weight('bold')
    marker = 'o'
    
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
        self.__x0 = float64(x0)
        self.__xn = float64(xn)
        self.__y0 = 0.000
        self.__y1 = float64(y1)
        self.__np = np
        self.__step = (xn-x0)/np
        self.__x = linspace(x0,xn,np+1)
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
                      
        plt.xticks(fontsize=10, fontproperties=self.font1)
        plt.yticks(fontsize=10, fontproperties=self.font1)
        plt.draw()
        plt.show()
        plt.savefig(name)
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
        self.__wfn = concatenate([yl[:-2],yr[1:]])
        return f

    def _find_xr(self,guess):
        "find the matching point"
        xr = secant(self.__Vxr, 1.0 , guess=guess)
        return xr
    

def check_sim(l):
    """check the sim list"""
    for x,y in l:
        i=0;j=0
        for a,b in l:
            if a==x: i+=1
            if b==y: j+=1
        if (i>1 or j>1): 
            print " --- !ERROR! --- "
            break
        
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

class VIB(UNITS):
    """
Represents vibrational analysis tool
    """
    def __init__(self,mol,hess,weight=True):
        self.hess = hess
        self.mol = mol
        self._prepare()
        if weight: self._weight()

    # public
    
    def eval(self):
        """find normal modes doing tANDr-projection-out"""
        # transformation matrix
        self.__u = self._projout()
        F = dot(transpose(self.__u),dot(self.hess,self.__u))
        #F = dot(self.__u,dot(self.hess,transpose(self.__u)))
        E = linalg.eigh(F)[0]
        # frequencies
        self.__freq = where(E>0.,
            sqrt( E)*self.HartreePerHbarToCmRec,
           -sqrt(-E)*self.HartreePerHbarToCmRec)
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
        M = zeros((N*3,N*3),dtype=float64)
        for i in xrange(N):
            m = 1./sqrt(self.masses[i])
            M[3*i+0,3*i+0] = m
            M[3*i+1,3*i+1] = m
            M[3*i+2,3*i+2] = m
        return M
    def _M1(self):
        """M^-1 matrix"""
        N = len(self.mol.atoms)
        M = zeros((N*3,N*3),dtype=float64)
        for i in xrange(N):
            m = sqrt(self.masses[i])
            M[3*i+0,3*i+0] = m
            M[3*i+1,3*i+1] = m
            M[3*i+2,3*i+2] = m
        return M        
    def _redmass(self):
        """calculate reduced masses and cartesian l matrix"""
        u_cart = dot(self._M(),self.__u)
        redmass = 1./sum(u_cart**2,axis=0)*self.ElectronMassToAmu
        # normalize u_cart
        u_cart = u_cart/sqrt(sum(u_cart**2,axis=0))
        return u_cart, redmass
    
    def _prepare(self):
        """prepare coordinates and masses"""
        coords = []; masses = []
        for atom in self.mol.atoms:
            coords.append(atom.pos())
            masses.append(atom.mass())
        self.masses = array(masses,dtype=float64)*self.AmuToElectronMass
        self.coords = array(coords,dtype=float64)
        return
    def _weight(self):
        """weight hessian"""
        M = self._M()
        self.hess = dot(M,dot(self.hess,M))
        return
    def _weight_old(self):
        """weight Hessian"""
        for i in xrange(len(self.mol.atoms)):
            mi = 1./sqrt(self.masses[i])
            for j in xrange(len(self.mol.atoms)):
                mj = 1./sqrt(self.masses[j])
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
        return sum(self.coords*self.masses[:,newaxis],axis=0)/sum(self.masses)
    def _r(self):
        """translate all coordinates to R_COM"""
        return self.coords - self._rcom()
    def _I(self):
        """moment of inertia tensor"""
        I = zeros((3,3),float64)
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
            m = sqrt(self.masses[i])
            v1 = [m,0,0]
            v2 = [0,m,0]
            v3 = [0,0,m]
            D1+=v1; D2+=v2; D3+=v3
        D1 = self._norm(array(D1,float64))
        D2 = self._norm(array(D2,float64))
        D3 = self._norm(array(D3,float64))
        return D1,D2,D3
    def _vec_r(self):
        """rotational vectors"""
        D4 = []; D5 = []; D6 = []
        X = linalg.eigh(self._I())[1]
        X1,X2,X3 = X
        r = self._r()
        for i in xrange(len(self.coords)):
            m = sqrt(self.masses[i])
            ri = r[i]
            P1 = dot(ri,X1)
            P2 = dot(ri,X2)
            P3 = dot(ri,X3)
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
        D4 = self._norm(array(D4,float64))
        D5 = self._norm(array(D5,float64))
        D6 = self._norm(array(D6,float64))
        return D4,D5,D6
    def _norm(self,u):
        """normalize vector u"""
        return u/sqrt(dot(u,u))
    def _gs(self,u,v):
        """Gram-Schmidt ortogonalization of vector v wrt u"""
        puv = dot(v,u)/dot(u,u)*u
        return v-puv
    def _pre_diag(self):
        """diagonalize Hessian matrix"""
        return linalg.eigh(self.hess)[1]
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
    
class Peak2DIR(UNITS):
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
    def __init__(self,x,y,z,Tw,
                 t_max=1.0,dt=0.008,
                 n_points=2**10,w_cent=1500.0,
                 wx_max=2000.0,wx_min=0.0,
                 wy_max=2000.0,wy_min=0.0):
        
        self.ToCmRec = 1.e-12*self.SpeedOfLight*100.0
        self.FromAngToCmRec = self.ToCmRec*2.*mPi
        dt *= 2.*mPi
        ### grids of data
        self.exp_grid = RBS(y,x,z)
        self.sim_grid = Grid2D(xmin=0.0,xmax=t_max*2.*mPi,
                               ymin=0.0,ymax=t_max*2.*mPi,
                               dx=dt,dy=dt)
        self.freq  = fft.fftshift( fft.fftfreq(n_points,
                                               d=dt*self.ToCmRec) ) + w_cent
        
        ### experimental data
        self.__kx = where(logical_and(self.freq<wx_max, self.freq>wx_min))[0]
        self.__ky = where(logical_and(self.freq<wy_max, self.freq>wy_min))[0]
        self.X = self.freq.copy()[self.__kx]
        self.Y = self.freq.copy()[self.__ky]
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
        
        if func_name=='r':
           if   n==1: self.func = self._r1_no_exch
           elif n==2: self.func = self._r2_no_exch
           
        self.n = n
    
    def get_parameters(self):
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
        return array(peaks,dtype=float64)

    def get_r2(self):
        """return R^2 coefficient of fitting"""
        data_av = average(self.Z)
        sse = sum((self.func(**self.args)-self.Z)**2)
        sst = sum((self.Z-data_av)**2)
        return 1 - sse/sst
        
    def fit(self,opts,method='slsqp',disp=1,bounds=[],ieqcons=[],
            epsilon=1e-06,pgtol=1e-012,factr=100.0,m=8000,
            approx_grad=True,fprime=None,maxfun=10000000,acc=1e-6,):
        """perform fitting using [options] list"""
        self.__set_param(opts)

        if method=='slsqp':
           result  = scipy.optimize.fmin_slsqp(self._residsq,
                                          self.init_list,iter=maxfun,
                                          acc=acc,disp=2,bounds=bounds,
                                          args=(opts,),full_output=True,
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
        return sum((self.Z - self.func(**self.args))**2.)
    
    ### Response Functions

    def __response(self,t3,t1,
                   Tw, w_01, mu_01, mu_12, anh, Delta, tau, T1, T2):
        """3-rd order response function for 3-level system (in time domain)"""
        w_off = -(w_01-self.w_cent) * self.FromAngToCmRec
        anh  *= self.FromAngToCmRec
        # di-Kubo line-shape function
        def g(t):
            G=0.
            for i in range(2):
                G+= tau[i]*Delta[i] * \
          ( exp(-t/tau[i])*tau[i] - tau[i] + t ) 
            G += t/T2
            return G
        #
        M = exp(-Tw/T1) * (1.0 + 0.8 * exp(-Tw/10000000.0))
        RR = 0.0; NR = 0.0
        
        mu_01_2= mu_01*mu_01
        mu_12_2= mu_12*mu_12
        
        gt1     = g(t1)
        gTw     = g(Tw)
        gt3     = g(t3)
        gt1Tw   = g(t1+Tw)
        gTwt3   = g(Tw+t3)
        gt1Twt3 = g(t1+Tw+t3)
        # --- Rephasing
        # R1 and R2
        RR+= mu_01_2*mu_01_2*exp(-1.j*w_off*(-t1+t3))*M*2.0*\
             exp(-gt1+gTw-gt3-gt1Tw-gTwt3+gt1Twt3)
        # R3
        RR-= mu_01_2*mu_12_2*exp(-1.j*(w_off*(-t1+t3)-anh*t3))*M*\
             exp(-gt1+gTw-gt3-gt1Tw-gTwt3+gt1Twt3)
        # --- Nonrephasing
        # R4 and R5
        NR+= mu_01_2*mu_01_2*exp(-1.j*w_off*(t1+t3))*M*2.0*\
             exp(-gt1-gTw-gt3+gt1Tw+gTwt3-gt1Twt3)
        # R6
        NR-= mu_01_2*mu_12_2*exp(-1.j*(w_off*(t1+t3)-anh*t3))*M*\
             exp(-gt1-gTw-gt3+gt1Tw+gTwt3-gt1Twt3)
        #
        return RR, NR

    def _r1_no_exch(self,w_01,anh,delta_1,delta_2,tau_1,tau_2,T1,T2):
        """3-rd order response without exchange and coupling"""
        
        ### signal in time-domain
        rr, nr = self.sim_grid.eval(self.__response, 
                                    # assumed parameters
                                    Tw=self.Tw, mu_01=1., mu_12=msqrt(2.),
                                    # optimizing parameters
                                    w_01=w_01, 
                                    anh=anh, 
                                    Delta=array([delta_1,delta_2]), 
                                    tau=array([tau_1,tau_2]), 
                                    T1=T1, T2=T2, )
        
        ### rephasing and non-rephasing spectras (2D FFT)
        data_rr_f = fft.fftshift( fft.fft2(rr,s=(self.__n_points,self.__n_points)) )
        data_nr_f = fft.fftshift( fft.fft2(nr,s=(self.__n_points,self.__n_points)) )
        data_rr_f = data_rr_f[:,::-1]
        data_rr_f = roll(data_rr_f,1,axis=1)
        
        ### total signal
        data_f = real(data_rr_f + data_nr_f)
        
        data_f = data_f[self.__kx,:]
        data_f = data_f[:,self.__ky]
        data_f/= data_f.max()
        data_f = data_f.transpose()
        
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
              l2 = " %6s"%(('%s'%let_greek.omega+'_01').rjust(6))
              l3 = " %6s"%(let_greek.Delta.rjust(6))
              l4 = " %6s"%('Peak'.rjust(6))
              l5 = " %6s"%('Peak'.rjust(6))
              l6 = " %6s"%('Peak'.rjust(6))
              l7 = " %6s"%('Peak'.rjust(6))
              l8 = " %6s"%('Peak'.rjust(6))
              l9 = " %6s"%('Peak'.rjust(6))
           for line in [l1,l2,l3,l4,l5,l6,l7,l8,l9]:
               log+= line + '\n'
           
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
        return array(fwhm,float64)
               
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
    
    ### Normal distribution

    def _normal(self,xo_1,sigma_1,A_1):
        """single Gaussian distribution"""
        return (A_1/(sigma_1*math.sqrt(2*math.pi)))\
               * exp(-(self.x-xo_1)**2/(2*sigma_1**2))
               
    ### pure Gaussian profiles
                           
    def _gauss1(self,xo_1,sigma_1,A_1):
        """single Gaussian distribution"""
        return A_1/sigma_1 * math.sqrt(4.*math.log(2.)/math.pi) * \
               exp(-4.*math.log(2.)/(sigma_1**2.)*(self.x-xo_1)**2.)

    def _gauss2(self,xo_1,sigma_1,A_1,
                     xo_2,sigma_2,A_2):
        """bimodal gaussian distribution"""
        g1 = A_1/sigma_1 * math.sqrt(4.*math.log(2.)/math.pi) * \
               exp(-4.*math.log(2.)/(sigma_1**2.)*(self.x-xo_1)**2.)
        g2 = A_2/sigma_2 * math.sqrt(4.*math.log(2.)/math.pi) * \
               exp(-4.*math.log(2.)/(sigma_2**2.)*(self.x-xo_2)**2.)
        return A_1*g1 + A_2*g2

    def _gauss3(self,xo_1,sigma_1,A_1,
                     xo_2,sigma_2,A_2,
                     xo_3,sigma_3,A_3):
        """trimodal gaussian distribution"""
        g1 = A_1/sigma_1 * math.sqrt(4.*math.log(2.)/math.pi) * \
               exp(-4.*math.log(2.)/(sigma_1**2.)*(self.x-xo_1)**2.)
        g2 = A_2/sigma_2 * math.sqrt(4.*math.log(2.)/math.pi) * \
               exp(-4.*math.log(2.)/(sigma_2**2.)*(self.x-xo_2)**2.)
        g3 = A_3/sigma_3 * math.sqrt(4.*math.log(2.)/math.pi) * \
               exp(-4.*math.log(2.)/(sigma_3**2.)*(self.x-xo_3)**2.)
        return A_1*g1 + A_2*g2 + A_3*g3

    def _gauss4(self,xo_1,sigma_1,A_1,
                     xo_2,sigma_2,A_2,
                     xo_3,sigma_3,A_3,
                     xo_4,sigma_4,A_4):
        """trimodal gaussian distribution"""
        g1 = A_1/sigma_1 * math.sqrt(4.*math.log(2.)/math.pi) * \
               exp(-4.*math.log(2.)/(sigma_1**2.)*(self.x-xo_1)**2.)
        g2 = A_2/sigma_2 * math.sqrt(4.*math.log(2.)/math.pi) * \
               exp(-4.*math.log(2.)/(sigma_2**2.)*(self.x-xo_2)**2.)
        g3 = A_3/sigma_3 * math.sqrt(4.*math.log(2.)/math.pi) * \
               exp(-4.*math.log(2.)/(sigma_3**2.)*(self.x-xo_3)**2.)
        g4 = A_4/sigma_4 * math.sqrt(4.*math.log(2.)/math.pi) * \
               exp(-4.*math.log(2.)/(sigma_4**2.)*(self.x-xo_4)**2.)

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

    ### pure Voigt profiles
    
    def __v(self,t,x,xc,wl,wg,a):
        A = math.exp(-t**2.)
        B = (math.sqrt(math.log(2.))*wl/wg)**2.
        C = (math.sqrt(4.*math.log(2.)) * (x-xc)/wg - t)**2.
        return a*2.*math.log(2.)/math.pi**(3./2.) * wl/wg**2 * A/(B+C)

    def _voigt1(self,xo_1,sigmaL_1,sigmaG_1,A_1):
        y = zeros(len(self.x),dtype=float64)
        for i in xrange(len(self.x)):
            y[i] = scipy.integrate.quad(self.__v,-inf,inf,full_output=0,args=(self.x[i],xo_1,sigmaL_1,sigmaG_1,A_1))[0]
        return y
    
    def _voigt2(self,xo_1,sigmaL_1,sigmaG_1,A_1,
                     xo_2,sigmaL_2,sigmaG_2,A_2):
        y = zeros(len(self.x),dtype=float64)
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

def reorder(P,sim,axis=0):
    """Reorders the tensor according to <axis> (default is 0). 
<sim> is the list of pairs from 'order' function. 
In normal numbers (starting from 1...)"""
    P_new = zeros(P.shape,dtype=float64)
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
    for i in xrange(len(R)-start):
        s = sum(sign(new_P[i])/sign(R[i]))
        if lprint: print "%10d %f" %(i+1,s)
        r_ = sum(( R[i+start]-new_P[i+start])**2)
        r__= sum((-R[i+start]-new_P[i+start])**2)
        if lprint:
         if s < -154: 
            print "TUTAJ s < -154"
            #new_P[i]*=-1.
         if r__<r_:
            print "TUTAJ r__<r_"
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
           Mol = PyQuante.Molecule(name,Coords,units='Bohr',
                                   multiplicity=mult,charge=charge,
                                   basis=basis,method=method)
    
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
           
       atnos = array(atnos,dtype=int)
       
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
                                   multiplicity=mult,charge=charge,
                                   basis=basis,method=method)
    
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
        field-= 2* tensordot(Qa[i],R,(0,0)) / Rab**5
        
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
    STR = array(STR,float64)
    # MONOPOLES
    l = d.readline();l = d.readline()
    while not l.startswith(' STOP'):
       a,b,c = l.split()
       chg.append(float64(b)+float64(c))
       l = d.readline()
    chg = array(chg,float64)
    # DIPOLES
    l = d.readline();l = d.readline()
    while not l.startswith(' STOP'):
        dip.append(l.split()[1:])
        l = d.readline()
    dip = array(dip,float64)
    # QUADRUPOLES
    l = d.readline();l = d.readline()
    while not l.startswith(' STOP'):
         q = l.split()[1:-1]
         l = d.readline()
         q+= l.split()
         qad.append(q)
         l = d.readline()
    qad = array(qad,float64)
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
    oct = array(oct,float64)
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

    STR=array(STR,float64).reshape(N,3)
    A = array(A,float64).reshape(N,3,3)
    return STR,A

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

    A = array(A,float64).reshape(N,3,3)
    return A

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
    C = array(C,dtype=float64).reshape(N,M)
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
    fock = zeros((nbasis,nbasis),dtype=float64)
    for i in xrange(g(nbasis)):
        line = data.readline()
        line = data.readline()
        nxses= array(line.split(),int)-1
        line = data.readline()

        for j in xrange(nbasis-i*5):
            line = data.readline()
            line = re_dbl_fort_c.sub(r'\1E\2', line)
            ny = int(line.split()[0])-1
            values = line.split()[4:]
            for k in xrange(len(values)):
                nx = nxses[k]
                v  = values[k]
                fock[nx,ny] = v
                fock[ny,nx] = v
    data.close()
    fock = array(fock,dtype=float64)
    return fock

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
    
    FC = array(FC,float64)
    H = zeros((N*3,N*3),dtype=float64)
    I = 0
    for i in xrange(N*3):
        for j in xrange(i+1):
            H[i,j] = FC[I]
            H[j,i] = H[i,j]
            I+=1
    
    return H

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
    fd = array(fd,float64).reshape(M*3,3)
    return fd

def Parse_EDS_InteractionEnergies(file):
    """parses EDS interaction energies from GAMESS log file"""
    
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

def ParseEFPInteractionEnergies(file):
    """parses EFP interaction energies from GAMESS log file"""
    
    dat  = open(file)
    data = dat.read()
    dat.close()
    
    querry = r'.*FRAGMENT-FRAGMENT INTERACTION ENERGIES.*'
    E = ['ELECTROSTATIC ENERGY','REPULSION ENERGY','POLARIZATION ENERGY',
         'DISPERSION ENERGY','CHARGE TRANSFER ENRGY','FINAL EFP ENERGY',]
         
    for term in E:
        querry+= '\s*%s\s+=\s+(%s).*\n' % (term,re_real)
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
             OQ  += 2* tensordot(tensordot(Oa[i],Qb[j],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qb - Oa  | R6
             QO  +=-2* tensordot(tensordot(Ob[j],Qa[i],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qa - Ob  | R6
             OQ  +=-14*tensordot(tensordot(R,tensordot(Oa[i],R,(1,0)),(0,0)) ,                      # Qb - Oa  | R6
                                 tensordot(R,Qb[j],(0,0)) ,(0,0)) / Rab**9                          
             QO  += 14*tensordot(tensordot(R,tensordot(Ob[j],R,(1,0)),(0,0)) ,                      # Qa - Ob  | R6
                                 tensordot(R,Qa[i],(0,0)) ,(0,0)) / Rab**9
             OQ  +=( 21*tensordot(tensordot(R,tensordot(Oa[i],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
                      * tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))) / Rab**11
             QO  +=(-21*tensordot(tensordot(R,tensordot(Ob[j],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
                      * tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0))) / Rab**11   
             OO  +=(2.)/(5.)*tensordot(Oa[i],Ob[j],((0,1,2),(0,1,2))) / Rab**7                      # Ob - Oa  | R7
             OO  +=(-42./5.)*tensordot(tensordot(R,Oa[i],(0,0)),
                                       tensordot(R,Ob[j],(0,0)),
                                       ((0,1),(0,1))) / Rab**9                                      # Ob - Oa  | R7
             OO  +=(189.)/(5.)*tensordot(
                                         tensordot(tensordot(R,Oa[i],(0,0)),R,(0,0)),
                                         tensordot(tensordot(R,Ob[j],(0,0)),R,(0,0)),
                                         (0,0)) /Rab**11                                            # Oa - Ob  | R7
             OO  +=-(231./5.)*(tensordot(tensordot(tensordot(R,Oa[i],(0,0)),R,(0,0)),R,(0,0)) *
                               tensordot(tensordot(tensordot(R,Ob[j],(0,0)),R,(0,0)),R,(0,0)) ) /\
                               Rab**13                                                              # Oa - Ob  | R7
             
             Eint = qq + qD + Dq + qQ + Qq + qO + Oq + DD + DQ + QD + DO + OD + QQ + QO + OQ + OO
             
             ### save the partitioning for current usage
             get_elmtp.qq = qq;get_elmtp.qD = qD;get_elmtp.qQ = qQ;get_elmtp.qO = qO;get_elmtp.QO = QO;
             get_elmtp.DD = DD;get_elmtp.DQ = DQ;get_elmtp.DO = DO;get_elmtp.QQ = QQ;get_elmtp.OO = OO;
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
    log+= " "+"-"*32+":"+"-"*26+"\n"
    log+= "\n"
    get_elmtp.log = log
    
    return get_elmtp.A, get_elmtp.B, get_elmtp.C, get_elmtp.D, get_elmtp.E

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
        self.shape = (self.nx,self.ny)
        # make 2D versions of the coordinate arrays
        # (needed for vectorized  function evaluators)
        #self.xcoorv = self.xcoor[:, newaxis]
        #self.ycoorv = self.ycoor[newaxis, :]
        #self.xcoorv, self.ycoorv = meshgrid(self.xcoor,self.ycoor)
        self.ycoorv, self.xcoorv = meshgrid(self.ycoor,self.xcoor)
    
    def eval(self,f,**kwargs):
        """Evaluate vectorized function f at each grid point"""
        return f(self.xcoorv,self.ycoorv,**kwargs)
    
class Grid3D:
    """represents 3D-grid of points"""
    def __init__(self,
                 xmin=0, xmax=1, dx=0.5,
                 ymin=0, ymax=1, dy=0.5,
                 zmin=0, zmax=1, dz=0.5):
        # coordinates in each space direction
        nx = int64((xmax-xmin)/dx + 1)
        ny = int64((ymax-ymin)/dy + 1)
        nz = int64((zmax-zmin)/dz + 1)
        
        x,y,z = mgrid[0:nx,0:ny,0:nz]
        
        # store for convenience
        self.dx = dx; self.dy = dy; self.dz = dz
        self.nx = nx; self.ny = ny; self.nz = nz
        self.shape = (self.nx,self.ny,self.nz)
        
        # make 3D versions of the coordinate arrays
        # (needed for vectorized  function evaluators)
        self.xcoorv = float64(x)*dx + xmin
        self.ycoorv = float64(y)*dy + ymin
        self.zcoorv = float64(z)*dz + zmin
            
    def eval(self,f,**kwargs):
        """Evaluate vectorized function f at each grid point"""
        return f(self.xcoorv,self.ycoorv,self.zcoorv,**kwargs)

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
