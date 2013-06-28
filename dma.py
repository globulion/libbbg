# ---------------------------------------------- #
#     DISTRIBUTED MULTIPOLE ANALYSIS OBJECT      #
# ---------------------------------------------- #

__all__ = ['DMA']

from numpy import zeros, float64, trace, array,\
                  tensordot, shape, outer, dot,\
                  transpose
from utilities2 import array_outer_product,    \
                       array_outer_product_1_2,\
                       array_outer_product_2_1,\
                       array_outer_product_1_n
import copy
from units import *

class DMA:
    """\
Represents the DMA distribution object. Inputs are q,m,T and O that mean
a set of distributed multipoles on n centers, ie.:                      
   Q = array([Q1,Q2, ... ,Qn])                                          
where Qi is i-th distributed tensor. The DMA class stores also vectors  
of position and origin of each of distributed center in a form of a     
NumPy array. As a zero object you can specify DMA(nfrag=<n>) after which
the DMA of n distributed centers (zero moments at each) will be created.
                                                                        
Usage:                                                                  
                                                                        
DMA(nfrag=<n>)                  return zero DMA object with <n> centers 
<object>.contract(<list>)       contract the DMA <object> using <list>  
                                that specifies the indices of origins to
                                be saved. Others will be removed.       
<object>.MakeUa(<list>,change_origin=<True>,contract=<False>)           
                                create united atoms within the DMA.     
                                unset 'change_origin' option when using 
                                distributed charges                     
<object>.Rotate(<rot>)          rotate the DMA about the rotation matrix
<object>.ChangeOrigin(new_origin_set=<0>,zero=<False>)                  
                                change the origins of distribited moments
                                to the new origin set specified by array
                                (in bohrs) of dimension (<n>,3). If you 
                                want to translate the centers to the    
                                origin of coordinate system specify only
                                zero=True.                              
<object>.ChangeUnits                                                    
<object>.MAKE_FULL()            create full DMA format                  
<object>.MakeTraceless()        transform the DMA quadrupoles and octu- 
                                poles to traceless form. Needs full for-
                                mat created first (above)               
                                                                        
                                                                        
"""
    
    def __init__(self,
                 nfrag=0,
                 name='Untitled',
                 q=0,m=0,T=0,O=0,
                 pos=zeros((1,3),dtype=float64),
                 origin=zeros((1,3),dtype=float64),
                 atoms=[Atom('X')]):
                    
        if nfrag>0:
           q=zeros((nfrag),dtype=float64)
           m=zeros((nfrag,3),dtype=float64)
           T=zeros((nfrag,6),dtype=float64)
           O=zeros((nfrag,10),dtype=float64)
        else:
           pass
        # DMA distribution in reduced format (GAMESS-like)
        self.DMA = [q,m,T,O]
        #self.pos = zeros((nfrag,3),dtype=float64)
        #self.origin = zeros((nfrag,3),dtype=float64)
        # name
        self.name = name
        # origin
        self.origin = origin
        # position array
        self.pos = pos
        # number of distributed sites (fragments)
        self.nfrag = len(self.DMA[0])
        # atomic sites if any
        if len(atoms)>1: self.atoms = atoms
        else: self.atoms = atoms*self.nfrag
        # if DMA FULL formatted memorial were created. Now it is not, so 'False'.
        self.full = False
        # if traceless forms were created. Now it is ordinary (primitive) DMA format so 'False'.
        self.traceless = False
    
    def contract(self,contrlist):
        """shrink the dimensions of DMA object by eliminating rows in attributes"""
        K = len(contrlist)
        self.nfrag = K
        self.DMA[0] = self.DMA[0][contrlist]
        self.DMA[1] = self.DMA[1][contrlist]
        self.DMA[2] = self.DMA[2][contrlist]
        self.DMA[3] = self.DMA[3][contrlist]
        self.origin = self.origin[contrlist]
        return
    
    def set_name(self,name):
        """sets the name for the object""" 
        self.name = name
        return
    
    def set_structure(self,pos=None,origin=None,atoms=None,equal=False):
        """sets new positions or origins and atoms"""
        # update positions
        if pos is not None:
           if len(self.pos) != len(pos): # it means you have to update self.atoms
              self.atoms = [ Atom('X') for x in range(len(pos)) ]
           self.pos = pos.copy()
        # update origins
        if origin is not None:
           self.origin = origin.copy()
        # update atoms
        if atoms is not None:
           self.atoms = [ Atom(x) for x in atoms.split(',') ]
        # equal the positions and origins
        if equal:
           self.origin = pos.copy()
           
    def write(self, file):
        """writes  the DMA distribution in a file"""
        newfile = open(file,'w')
        log = ' SLV output file. All units in AU\n\n'
        log+= self.__repr__()
        log += ' Structure\n'
        log += ' ---------\n\n'
        for i,atom in enumerate(self.atoms):
            log+= " %-10s %14.8f %14.8f %14.8f\n" % (atom.symbol, self.pos[i,0], self.pos[i,1], self.pos[i,2])
        log+='\n'
        
        # print origins in case the where different that structure
        if (self.pos.shape != self.origin.shape or not (self.pos==self.origin).all() ):
           log += ' Origins\n'
           log += ' -------\n\n'
           for i,point in enumerate(self.origin):
               log+= " %-10s %14.8f %14.8f %14.8f\n" % ('X', point[0], point[1], point[2])
           log+='\n'
        
        newfile.write(log)
        newfile.close()
        
    def __getitem__(self,index): 
        return array(self.DMA[index])

    def __setitem__(self,index,value):
        self.DMA[index] = array(value)
    
    def __len__(self):
        return (len(self.DMA[0]))
                
    def copy(self):
        return copy.deepcopy(self)

    def __add__(self,other):
        if (isinstance(self,DMA) and (not isinstance(other,DMA))):
           q = self.DMA[0] + other
           m = self.DMA[1] + other
           T = self.DMA[2] + other
           O = self.DMA[3] + other
        elif ((not isinstance(self,DMA)) and isinstance(other,DMA)):
           q = other.DMA[0] + self
           m = other.DMA[1] + self
           T = other.DMA[2] + self
           O = other.DMA[3] + self
        else:
           q = other.DMA[0] + self.DMA[0]
           m = other.DMA[1] + self.DMA[1]
           T = other.DMA[2] + self.DMA[2]
           O = other.DMA[3] + self.DMA[3]     
        return DMA(q=q,m=m,T=T,O=O)

    def __sub__(self,other):
        q = self.DMA[0] - other.DMA[0]
        m = self.DMA[1] - other.DMA[1]
        T = self.DMA[2] - other.DMA[2]
        O = self.DMA[3] - other.DMA[3]
        return DMA(q=q,m=m,T=T,O=O)
    
    def __div__(self,other):
        if (isinstance(self,DMA) and (not isinstance(other,DMA))):
           q = self.DMA[0] / other
           m = self.DMA[1] / other
           T = self.DMA[2] / other
           O = self.DMA[3] / other
        elif ((not isinstance(self,DMA)) and isinstance(other,DMA)):
           q = other.DMA[0] / self
           m = other.DMA[1] / self
           T = other.DMA[2] / self
           O = other.DMA[3] / self
        else:
           q = self.DMA[0] / other.DMA[0]
           m = self.DMA[1] / other.DMA[1]
           T = self.DMA[2] / other.DMA[2]
           O = self.DMA[3] / other.DMA[3]     

        return DMA(q=q,m=m,T=T,O=O)

    def __mul__(self,other):
        if (isinstance(self,DMA) and (not isinstance(other,DMA))):
           q = self.DMA[0] * other
           m = self.DMA[1] * other
           T = self.DMA[2] * other
           O = self.DMA[3] * other
        elif ((not isinstance(self,DMA)) and isinstance(other,DMA)):
           q = other.DMA[0] * self
           m = other.DMA[1] * self
           T = other.DMA[2] * self
           O = other.DMA[3] * self
        else:
           q = other.DMA[0] * self.DMA[0]
           m = other.DMA[1] * self.DMA[1]
           T = other.DMA[2] * self.DMA[2]
           O = other.DMA[3] * self.DMA[3]            
        return DMA(q=q,m=m,T=T,O=O)

    def __rmul__(self,other):
        if (isinstance(self,DMA) and (not isinstance(other,DMA))):
           q = self.DMA[0] * other
           m = self.DMA[1] * other
           T = self.DMA[2] * other
           O = self.DMA[3] * other
        elif ((not isinstance(self,DMA)) and isinstance(other,DMA)):
           q = other.DMA[0] * self
           m = other.DMA[1] * self
           T = other.DMA[2] * self
           O = other.DMA[3] * self
        return DMA(q=q,m=m,T=T,O=O)

    def __neg__(self):
        q = -self.DMA[0]
        m = -self.DMA[1]
        T = -self.DMA[2]
        O = -self.DMA[3]
        return DMA(q=q,m=m,T=T,O=O)

    def MAKE_FULL(self):
        """creates arrays of ordinary forms of distributed multipoles
        in full format."""
        CHARGES = array(self.DMA[0])
        DIPOLES = array(self.DMA[1])
        QDPOLES = zeros((self.nfrag,3,3),dtype=float64)
        OCTPLES = zeros((self.nfrag,3,3,3),dtype=float64)
        for site in range(self.nfrag):
            # add quadrupoles
            q = array(self.DMA[2])
            QDPOLES[site,0,0] = q[site,0]
            QDPOLES[site,1,1] = q[site,1]
            QDPOLES[site,2,2] = q[site,2]
            QDPOLES[site,0,1] = q[site,3]
            QDPOLES[site,0,2] = q[site,4]
            QDPOLES[site,1,2] = q[site,5]
            #
            QDPOLES[site,1,0] = QDPOLES[site,0,1]
            QDPOLES[site,2,0] = QDPOLES[site,0,2]
            QDPOLES[site,2,1] = QDPOLES[site,1,2]
            
            #add octupoles
            o = array(self.DMA[3])
            OCTPLES[site,0,0,0] = o[site,0]
            OCTPLES[site,1,1,1] = o[site,1]
            OCTPLES[site,2,2,2] = o[site,2]
            OCTPLES[site,0,0,1] = o[site,3]
            OCTPLES[site,0,0,2] = o[site,4]
            OCTPLES[site,0,1,1] = o[site,5] 
            OCTPLES[site,1,1,2] = o[site,6]
            OCTPLES[site,0,2,2] = o[site,7]
            OCTPLES[site,1,2,2] = o[site,8]
            OCTPLES[site,0,1,2] = o[site,9]
            # 
            OCTPLES[site,0,1,0] = OCTPLES[site,0,0,1] 
            OCTPLES[site,1,0,0] = OCTPLES[site,0,0,1] 
            #
            OCTPLES[site,0,2,0] = OCTPLES[site,0,0,2] 
            OCTPLES[site,2,0,0] = OCTPLES[site,0,0,2] 
            #
            OCTPLES[site,1,0,1] = OCTPLES[site,0,1,1] 
            OCTPLES[site,1,1,0] = OCTPLES[site,0,1,1] 
            #
            OCTPLES[site,1,2,1] = OCTPLES[site,1,1,2]
            OCTPLES[site,2,1,1] = OCTPLES[site,1,1,2]
            #
            OCTPLES[site,2,0,2] = OCTPLES[site,0,2,2]
            OCTPLES[site,2,2,0] = OCTPLES[site,0,2,2]
            #
            OCTPLES[site,2,1,2] = OCTPLES[site,1,2,2]
            OCTPLES[site,2,2,1] = OCTPLES[site,1,2,2]
            #
            OCTPLES[site,0,2,1] = OCTPLES[site,0,1,2]
            OCTPLES[site,2,0,1] = OCTPLES[site,0,1,2]
            #
            OCTPLES[site,1,0,2] = OCTPLES[site,0,1,2]
            #
            OCTPLES[site,1,2,0] = OCTPLES[site,0,1,2]
            OCTPLES[site,2,1,0] = OCTPLES[site,0,1,2]
            
        self.DMA_FULL = [ self.pos , CHARGES , DIPOLES , QDPOLES , OCTPLES ]
        self.full = True

    def makeDMAfromFULL(self):
        """saves """
        self.pos = array(self.DMA_FULL[0])
        self.DMA[0] = array(self.DMA_FULL[1])
        self.DMA[1] = array(self.DMA_FULL[2])
        #
        self.DMA[2][:,0] = array(self.DMA_FULL[3][:,0,0])
        self.DMA[2][:,1] = array(self.DMA_FULL[3][:,1,1])
        self.DMA[2][:,2] = array(self.DMA_FULL[3][:,2,2])
        self.DMA[2][:,3] = array(self.DMA_FULL[3][:,0,1])
        self.DMA[2][:,4] = array(self.DMA_FULL[3][:,0,2])
        self.DMA[2][:,5] = array(self.DMA_FULL[3][:,1,2])
        #
        self.DMA[3][:,0] = array(self.DMA_FULL[4][:,0,0,0])
        self.DMA[3][:,1] = array(self.DMA_FULL[4][:,1,1,1])
        self.DMA[3][:,2] = array(self.DMA_FULL[4][:,2,2,2])
        self.DMA[3][:,3] = array(self.DMA_FULL[4][:,0,0,1])
        self.DMA[3][:,4] = array(self.DMA_FULL[4][:,0,0,2])
        self.DMA[3][:,5] = array(self.DMA_FULL[4][:,0,1,1])
        #
        self.DMA[3][:,6] = array(self.DMA_FULL[4][:,1,1,2])
        self.DMA[3][:,7] = array(self.DMA_FULL[4][:,0,2,2])
        self.DMA[3][:,8] = array(self.DMA_FULL[4][:,1,2,2])
        self.DMA[3][:,9] = array(self.DMA_FULL[4][:,0,1,2])
        #
        del self.DMA_FULL
        self.full = False
        
    def MakeTraceless(self):
        """turns ordinary full-formatted DMA_FULL into traceless 
        full-formatted DMA_FULL. The traceless form is taken to be
        following Buckingham /add citation!!!"""
        
        if self.traceless:
           raise  Exception("\nerror: DMA is already traceless!!! quitting ...\n")
           
        if self.full:
           # Quadrupole moment
           for Q in self.DMA_FULL[3]:
               t=trace(Q)
               Q *= (3./2.)
               for i in [0,1,2]: Q[i,i]-=(1./2.)*t
           # Octapole moment
           for O in self.DMA_FULL[4]:
               W= zeros((3,3,3),dtype=float)
               W[:]= O
               O *= (5./2.)
               for i in [0,1,2]:
                   for j in [0,1,2]:
                       for k in [0,1,2]:
                           Wt = trace(W)
                           if   i==j and i==k:
                                O[i,j,k]-= (1./2.) * (Wt[i] + Wt[j] + Wt[k])
                           elif i==j and i!=k:
                                O[i,j,k]-= (1./2.) * Wt[k]
                           elif j!=k and i==k:
                                O[i,j,k]-= (1./2.) * Wt[j]
                           elif i!=j and j==k:
                                O[i,j,k]-= (1./2.) * Wt[i]
                             
           self.traceless = True
           
        else: raise Exception("\nerror: no FULL DMA object created! quitting...\n")

    def ChangeOrigin(self,new_origin_set=0,zero=False):
        """shifts origin(s) of ordinary full-formatted DMA_FULL
           into new one(s)"""
           
        if self.full:
           if zero:
              new_origin_set = zeros((self.nfrag,3),dtype=float64)
           #print new_origin_set
           #new_origin_set*=-1 
           old_origin_set = array(self.origin)#array(self.DMA_FULL[0])
           M0 = array(self.DMA_FULL[1])
           M1 = array(self.DMA_FULL[2]) 
           M2 = array(self.DMA_FULL[3]) 
           M3 = array(self.DMA_FULL[4])

           #A2 = zeros((self.nfrag,3,3),dtype=float64)
           #B2 = zeros((self.nfrag,3,3),dtype=float64)
           A2 = array_outer_product(old_origin_set,old_origin_set)
           B2 = array_outer_product(new_origin_set,new_origin_set)

           D = new_origin_set - old_origin_set
           
           #print old_origin_set
           M1_o = M1 + array_outer_product_1_n(M0, old_origin_set)
           #AM = zeros((self.nfrag,3,3),dtype=float64)
           #MA = zeros((self.nfrag,3,3),dtype=float64)
           AM = array_outer_product(old_origin_set,M1_o)
           MA = array_outer_product(M1_o,old_origin_set)
           
           M2_o = M2 - array_outer_product_1_n(M0, A2) + AM + MA
           
           # First moments
           new_M1 = M1 - array_outer_product_1_n(M0, new_origin_set- old_origin_set)
           # Second moments
           #MD = zeros((self.nfrag,3,3),dtype=float64)
           MD = array_outer_product(M1_o,D)
           #DM = zeros((self.nfrag,3,3),dtype=float64)
           DM =  array_outer_product(D,M1_o)
           new_M2 = M2 + array_outer_product_1_n(M0,(B2-A2)) - MD - DM

           # Third moments
           A3 = array_outer_product_1_2(old_origin_set,A2)#
           B3 = array_outer_product_1_2(new_origin_set,B2)#
           M1A= array_outer_product(M1_o,old_origin_set)#
           M1B= array_outer_product(M1_o,new_origin_set)#
           #BAM2 = zeros((self.nfrag,3,3,3),dtype=float64)
           BAM2  = array_outer_product_1_2(D,M2_o)#
           B2A2M1= array_outer_product_2_1(B2-A2,M1_o)#
           QB3A3 = array_outer_product_1_n(M0,(B3-A3))#
           M2BA  = array_outer_product_2_1(M2_o,D)#
           M1B2A2= array_outer_product_1_2(M1_o,B2-A2)#
           AM1A  = array_outer_product_1_2(old_origin_set,M1A)#
           BM1B  = array_outer_product_1_2(new_origin_set,M1B)#
           #M1AM1 = array_outer_product_2_1(M1A,old_origin_set)
           #M1BM1 = array_outer_product_2_1(M1B,new_origin_set)
           #print shape(array_outer_product_1_2(new_origin_set,M2_o))
           BM2TR = transpose(array_outer_product_1_2(new_origin_set,M2_o),
                             axes=(0,2,1,3))#
           AM2TR = transpose(array_outer_product_1_2(old_origin_set,M2_o),
                             axes=(0,2,1,3))#

           new_M3 = M3   - BAM2 + B2A2M1 - QB3A3 - M2BA + M1B2A2 + \
                    BM1B - AM1A - BM2TR  + AM2TR
                      
           if self.nfrag == 1:
              self.origin = new_origin_set[0]
           
           self.DMA_FULL.append(new_origin_set)
           self.DMA_FULL[2] = new_M1
           self.DMA_FULL[3] = new_M2
           self.DMA_FULL[4] = new_M3
           
           self.makeDMAfromFULL()
           self.origin = array(new_origin_set)
           
        else: raise Exception("\nerror: no FULL DMA object created! quitting...\n")

    def MakeUa(self,ua_list,change_origin=True,contract=False):
        """transforms the object to the contracted DMA form employing united atoms. 
        Usage:
        MakeUa( ua_list )
        
        where 'ua_list' is a list of tuples ti containing integers:
        
        ua+list = [ t1, t2, t3 , ... ], 
        ti = ( A, a1, a2, a3, ... )
        
        For i-th tuple corresponding to i-th UA the first integer A
        is the atom id of UA center atom and the remaining a1, a2 etc relate to 
        the atoms supposed to be contracted within a united atom UA for A-th center"""
        
        origin = self.origin
        if change_origin:
           self.MAKE_FULL()
           self.ChangeOrigin(zero=True)
        for t in ua_list:
            for i in t[1:]:
                for c in [0,1,2,3]:
                    self.DMA[c][t[0]-1]+= self.DMA[c][i-1]
                    if not c: self.DMA[c][i-1] = 0
                    else:     self.DMA[c][i-1].fill(0)
        if change_origin:
           self.MAKE_FULL()
           self.ChangeOrigin(new_origin_set=origin)
      
        if contract:
           contrlist = self.__contrListFromUaList(ua_list)
           self.contract(contrlist)
        
    def __contrListFromUaList(self,ua_list):
       """creates the contraction list from ua_list"""
       return []
    
    def Rotate(self,rotmat):
        """rotates the ordinary full-formatted DMA_FULL
           in molecular space based on rotation matrix
           given"""

        #self.MAKE_FULL()
        if self.full:
           ### transform the dipoles
           #self.DMA[1] = dot(self.DMA[1],transpose(rotmat))
           self.DMA_FULL[2] = dot(self.DMA[1],rotmat)


           ### transform the quadrupoles
           quadrupole = zeros((self.nfrag,3,3),dtype=float64)
           for i in range(self.nfrag):
               quadrupole[i] = dot(transpose(rotmat),dot( self.DMA_FULL[3][i] ,rotmat))
               #quadrupole[i] = tensordot( rot, tensordot( rotmat, self.DMA_FULL[3][i],(1,1)),(1,1))

           #self.DMA[2][:,0] = quadrupole[:,0,0]
           #self.DMA[2][:,1] = quadrupole[:,1,1]
           #self.DMA[2][:,2] = quadrupole[:,2,2]
           #self.DMA[2][:,3] = quadrupole[:,0,1]
           #self.DMA[2][:,4] = quadrupole[:,0,2]
           #self.DMA[2][:,5] = quadrupole[:,1,2]
           self.DMA_FULL[3] = quadrupole
           #del quadrupole
        
           ### transform the octupoles
           octupole = zeros((self.nfrag,3,3,3),dtype=float64)
           for i in range(self.nfrag):
               octupole[i] = tensordot(transpose(rotmat),tensordot(
                                       transpose(rotmat),tensordot(
                                       transpose(rotmat),self.DMA_FULL[4][i] ,(1,2)
                                     ),(1,2)
                                     ),(1,2)
                                     )

           #self.DMA[3][:,0] = octupole[:,0,0,0]
           #self.DMA[3][:,1] = octupole[:,1,1,1]
           #self.DMA[3][:,2] = octupole[:,2,2,2]
           #self.DMA[3][:,3] = octupole[:,0,0,1]
           #self.DMA[3][:,4] = octupole[:,0,0,2]
           #self.DMA[3][:,5] = octupole[:,0,1,1]
           #self.DMA[3][:,6] = octupole[:,1,1,2]
           #self.DMA[3][:,7] = octupole[:,0,2,2]
           #self.DMA[3][:,8] = octupole[:,1,2,2]
           #self.DMA[3][:,9] = octupole[:,0,1,2]
           self.DMA_FULL[4] = octupole
           #del octupole
        
           self.makeDMAfromFULL()
           
        else: raise Exception("\nerror: no FULL DMA object created! quitting...\n")

    def ChangeUnits(self,charges,dipoles,quadrupoles,octupoles):
        """changes the units"""
        self.DMA[0] *= charges
        self.DMA[1] *= dipoles
        self.DMA[2] *= quadrupoles
        self.DMA[3] *= octupoles
        
    def OverallMoments_old(self,origin=zeros(3,dtype=float64)):
        """calculates overall primitive moments from charges.
           This tool has a purpose of testing population analysis obtained by
           fitting to the molecular ab initio potential or other methods"""
        overall = DMA(nfrag=1)
        overall.name = 'Test of reproducing multipoles from charges'
        overall.pos = zeros((1,3),dtype=float64)
        overall.pos[0] = origin
        
        ### make full format of DMA
        self.MAKE_FULL()
        
        ### compute molecular solvatochromic moments
        mu   = zeros((3),dtype=float64)
        quad = zeros((3,3),dtype=float64)
        oct  = zeros((3,3,3),dtype=float64)
        for atom in range(self.nfrag):
            r     = self.pos[atom] - origin
            qatom = self.DMA[0][atom]
            ### calculate dipole moment
            mu   += qatom * r 
            ### calculate quadrupole moment
            quad += qatom * outer (r,r)
            ### calculate octupole moment
            oct  += qatom * outer( r, outer (r,r) ).reshape(3,3,3)

        ### set the molecular moments into the DMA solvent object
        overall.DMA[1][0] = mu
           
        overall.DMA[2][0,0] = quad[0,0]           
        overall.DMA[2][0,1] = quad[1,1]
        overall.DMA[2][0,2] = quad[2,2]
        overall.DMA[2][0,3] = quad[0,1]
        overall.DMA[2][0,4] = quad[0,2]
        overall.DMA[2][0,5] = quad[1,2] 
           
        overall.DMA[3][0,0] = oct[0,0,0]
        overall.DMA[3][0,1] = oct[1,1,1]
        overall.DMA[3][0,2] = oct[2,2,2]
        overall.DMA[3][0,3] = oct[0,0,1]
        overall.DMA[3][0,4] = oct[0,0,2]
        overall.DMA[3][0,5] = oct[0,1,1]
        overall.DMA[3][0,6] = oct[1,1,2]
        overall.DMA[3][0,7] = oct[0,2,2]
        overall.DMA[3][0,8] = oct[1,2,2]
        overall.DMA[3][0,9] = oct[0,1,2]
        
        #overall.DMA[1] *= UNITS.BohrElectronToDebye
        #overall.DMA[2] *= UNITS.BohrElectronToDebye * UNITS.BohrToAngstrom
        #overall.DMA[2] *= UNITS.BohrElectronToDebye * UNITS.BohrToAngstrom**2
          
        return overall
    
    
    def OverallMoments(self,origin=zeros(3,dtype=float64)):
        """calculates overall primitive moments from charges.
           This tool has a purpose of testing population analysis obtained by
           fitting to the molecular ab initio potential or other methods"""
        overall = DMA(nfrag=1)
        #overall.name = 'Test of reproducing multipoles from charges [units: Debyes(*Angstroms^n)]'
        overall.name = 'Test of reproducing overall multipoles from CAMM [units: A.U.]'
        overall.pos = zeros((1,3),dtype=float64)
        overall.pos[0] = origin
        self.MAKE_FULL()
        self.ChangeOrigin(zero=1)
        
        ### make full format of DMA
        self.MAKE_FULL()
        
        ### compute molecular solvatochromic moments
        mu   = zeros((3),dtype=float64)
        quad = zeros((3,3),dtype=float64)
        oct  = zeros((3,3,3),dtype=float64)
        for atom in range(self.nfrag):
            r     = self.origin[atom]### zmiana origin z pos!!!
            qatom = self.DMA[0][atom]
            ### calculate dipole moment
            mu   += self.DMA[1][atom]
            ### calculate quadrupole moment
            quad += self.DMA_FULL[3][atom]
            ### calculate octupole moment
            oct  += self.DMA_FULL[4][atom]

        ### set the molecular moments into the DMA solvent object
        overall.DMA[1][0] = mu
           
        overall.DMA[2][0,0] = quad[0,0]           
        overall.DMA[2][0,1] = quad[1,1]
        overall.DMA[2][0,2] = quad[2,2]
        overall.DMA[2][0,3] = quad[0,1]
        overall.DMA[2][0,4] = quad[0,2]
        overall.DMA[2][0,5] = quad[1,2] 
           
        overall.DMA[3][0,0] = oct[0,0,0]
        overall.DMA[3][0,1] = oct[1,1,1]
        overall.DMA[3][0,2] = oct[2,2,2]
        overall.DMA[3][0,3] = oct[0,0,1]
        overall.DMA[3][0,4] = oct[0,0,2]
        overall.DMA[3][0,5] = oct[0,1,1]
        overall.DMA[3][0,6] = oct[1,1,2]
        overall.DMA[3][0,7] = oct[0,2,2]
        overall.DMA[3][0,8] = oct[1,2,2]
        overall.DMA[3][0,9] = oct[0,1,2]
        
        #overall.DMA[1] *= UNITS.BohrElectronToDebye
        #overall.DMA[2] *= UNITS.BohrElectronToDebye * UNITS.BohrToAngstrom
        #overall.DMA[3] *= UNITS.BohrElectronToDebye * UNITS.BohrToAngstrom**2
          
        return overall
        
    def __repr__(self):
        """PRINTOUT FORM FOR DMA REDUCED FORMAT"""
        log = "\n"
        log+= " --- %s ---\n\n" % (self.name)
        log+= " Distributed zeroth-order property\n"
        log+= " ---------------------------------\n"
        log+= "\n"
        charges = self.DMA[0]
        for xfrag in charges:
            log+= "%16.6E\n"%xfrag
        log+= "\n"

        log+= " Distributed first-order property \n"
        log+= " ---------------------------------\n"
        log+= "\n"
        dipoles = self.DMA[1]
        log+= "%16s %16s %16s\n" % ('X'.rjust(16),'Y'.rjust(16),'Z'.rjust(16))
        for xfrag in dipoles:
            log+= "%16.6E %16.6E %16.6E\n"%(xfrag[0],xfrag[1],xfrag[2])
        log+= "\n"  

        log+= " Distributed second-order property\n"
        log+= " ---------------------------------\n"
        log+= "\n"
        quadrupoles = self.DMA[2]
        log+= "%16s %16s %16s %16s %16s %16s\n" %\
             ('XX'.rjust(16),'YY'.rjust(16),'ZZ'.rjust(16),
              'XY'.rjust(16),'XZ'.rjust(16),'YZ'.rjust(16))
        for xfrag in quadrupoles:
            log+= "%16.6E %16.6E %16.6E %16.6E %16.6E %16.6E\n"%\
             (xfrag[0],xfrag[1],xfrag[2],
              xfrag[3],xfrag[4],xfrag[5])
        log+= "\n"

        log+= " Distributed third-order property \n"
        log+= " ---------------------------------\n"
        log+= "\n"
        octupoles = self.DMA[3]
        log+= "%16s %16s %16s %16s %16s %16s\n" %\
            ('XXX'.rjust(16),'YYY'.rjust(16),'ZZZ'.rjust(16),
             'XXY'.rjust(16),'XXZ'.rjust(16),'XYY'.rjust(16))
        log+= "%16s %16s %16s %16s\n" %\
            ('YYZ'.rjust(16),'XZZ'.rjust(16),'YZZ'.rjust(16),'XYZ'.rjust(16))
        for xfrag in octupoles:
            log+= "%16.6E %16.6E %16.6E %16.6E %16.6E %16.6E\n" %\
                   (xfrag[0],xfrag[1],xfrag[2],
                    xfrag[3],xfrag[4],xfrag[5])
            log+= "%16.6E %16.6E %16.6E %16.6E\n" %\
                   (xfrag[6],xfrag[7],xfrag[8],xfrag[9])
        log+= " "+"-"*100+"\n"
        if self.nfrag == 1:
           if self.traceless: 
              log+= (" traceless form,   origin: %s\n"  % str(self.origin*0.5291772086)[1:-1]).rjust(100)  
           else: 
              log+= (" primitive form,   origin: %s \n" % str(self.origin*0.5291772086)[1:-1]).rjust(100) 
        else:
           if self.traceless: 
              log+= " traceless form\n".rjust(100)  
           else: 
              log+= " primitive form\n".rjust(100) 
        log+= "\n\n\n"
        
        return str(log) 
