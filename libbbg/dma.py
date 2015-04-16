# ---------------------------------------------- #
#     DISTRIBUTED MULTIPOLE ANALYSIS MODULE      #
# ---------------------------------------------- #

__all__ = ['DMA']

from numpy import zeros, float64, trace, array,\
                  tensordot, shape, outer, dot,\
                  transpose, sqrt, sum, linalg
from utilities2 import array_outer_product,    \
                       array_outer_product_1_2,\
                       array_outer_product_2_1,\
                       array_outer_product_1_n,\
                       array_outer_product_2_2,\
                       array_outer_product_1_3,\
                       array_outer_product_3_1
from units import *
import copy

def ParseDMA(file,type='coulomb'):
    """\
>>>>> Copied from LIBBBG.utilities <<<<<
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
                     origin=Origin,
                     is_traceless=is_traceless )
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
      
def interchange(T,ind):
    """\
>>>>> Copied from LIBBBG.utilities <<<<<
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
 
class DMA:
    """\
==============================LIBBBG:DMA================================
Represents the DMA distribution object. Inputs are q,m,T and O that mean
a set of distributed multipoles on n centers, ie.:                      
   Q = array([Q1,Q2, ... ,Qn])                                          
where Qi is i-th distributed tensor. The DMA class stores also vectors  
of position and origin of each of distributed center in a form of a     
NumPy array. As a zero object you can specify DMA(nfrag=<n>) after which
the DMA of n distributed centers (zero moments at each) will be created.
                                                                        
========================================================================
Usage:                                                                  
------------------------------------------------------------------------
DMA(nfrag=<n>)                  return zero DMA object with <n> centers 
<object>.contract(<list>)       contract the DMA <object> using <list>  
                                that specifies the indices of origins to
                                be saved. Others will be removed.       
<object>.MakeUa(<list>,change_origin=<True>,contract=<False>)           
                                create united atoms within the DMA.     
                                unset 'change_origin' option when using 
                                distributed charges                     
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
<object>.Rotate(<rot>)          rotate the DMA about the rotation matrix
                                Needs full format created first (above) 
<object>.copy()                 return a deepcopy of <object>           
<object>.write(<file>)          write DMA distribution to a <file> in cf
<object>.interchange(<ind>)     interchange all memorials using ind list
------------------------------------------------------------------------
"get" methods:                                                          
------------------------------------------------------------------------
<object>.get_pos()              return atomic positions ndarray (in bohr)
<object>.get_origin()           return origins ndarray (in bohr)       
<object>.get_name()             return name of distribution             
<object>.get_nfrags()           return number of distributed sites      
<object>.get_natoms()           return number of atoms                  
------------------------------------------------------------------------
"set" methods:                                                          
------------------------------------------------------------------------
<object>.set_name(<name>)       set the name of distribution            
<object>.set_moments(charges=<None>,dipoles=<None>,                     
                     quadrupoles=<None>,octupoles=<None>)               
                                set moments as numpy arrays in reduced  
                                format                                  
<object>.set_structure(pos=<None>,origin=<None>,                        
                       atoms=<None>,equal=<False>)                      
                                set positions and/or origins and/or atom
                                symbols. If equal=<True> setting only po
                                sitions will set the same origins too   
                                without explicit call for origins       
------------------------------------------------------------------------
mathematical operations:                                                
------------------------------------------------------------------------
1) addition                                                             
   <object> + <number> = adds number to all moments (new <object>)      
   <object> + <object> = adds moments of two objects (new <object>)     
   <object> - <number> = substract number to all moments (new <object>)
   <object> - <object> = substract moments from each other (new <object>)
   <object> * <number> = multiply all moments by number (new <object>)  
   <object> / <number> = divide all moments by number (new <object>)    
   - <object>          = return negative moments (new <object>)         
   len(object)         = return number of distributed sites             
------------------------------------------------------------------------
========================================================================
"""
    
    # ---- CONSTRUCTORS

    def __init__(self,
                 file=None,
                 nfrag=None,
                 name='Untitled',
                 q=None,m=None,T=None,O=None,H=None,
                 pos=zeros((1,3),dtype=float64),
                 origin=zeros((1,3),dtype=float64),
                 atoms=None,
                 is_traceless=False,
                 hexadecapoles=None):

        # save the flag whether this object contains hexadecapoles        
        if hexadecapoles is not None:
           self.has_hexadecapoles = hexadecapoles
        else:
           if H is None: self.has_hexadecapoles = False
           else        : self.has_hexadecapoles = True

        # determine the path of DMA object construction
        if file is not None:
           self.__call__(file)

        elif nfrag is not None:
             q=zeros((nfrag),dtype=float64)
             m=zeros((nfrag,3),dtype=float64)
             T=zeros((nfrag,6),dtype=float64)
             O=zeros((nfrag,10),dtype=float64)
             if self.has_hexadecapoles: H=zeros((nfrag,15),dtype=float64)
             pos = zeros((nfrag,3),dtype=float64)
             origin = zeros((nfrag,3),dtype=float64)

        if file is None:
           # DMA distribution in reduced format (GAMESS-like)
           if self.has_hexadecapoles: self.DMA = [q,m,T,O,H]
           else:                      self.DMA = [q,m,T,O]
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
           if atoms is None:
              #if len(atoms)>1: self.atoms = [Atom('X')]
              #else: self.atoms = [Atom('X')]*self.nfrag
              self.atoms = [Atom('X')]*self.nfrag
           else: self.atoms=atoms
           # if DMA FULL formatted memorial were created. Now it is not, so 'False'.
           self.full = False
           # if traceless forms were created. Now it is ordinary (primitive) DMA format so 'False'.
           self.traceless = is_traceless
 
    def __call__(self,file):
        """open DMA file"""
        dma = ParseDMA(file,'c')
        if dma.has_hexadecapoles:
           self.DMA = [dma[x] for x in [0,1,2,3,4]]
        else:
           self.DMA = [dma[x] for x in [0,1,2,3]]
        self.origin = dma.get_origin()
        self.pos   = dma.get_pos()
        self.name  = dma.get_name()
        self.nfrag = len(dma[0])
        self.atoms = dma.atoms
        self.full  = False
        self.traceless = False
        self.has_hexadecapoles = dma.has_hexadecapoles
        return

    # ---- GET methods

    def get_charges      (self): return self[0]
    def get_dipoles      (self): return self[1]
    def get_quadrupoles  (self): return self[2]
    def get_octupoles    (self): return self[3]
    def get_hexadecapoles(self): return self[4]
    
    def get_pos(self):
        """return positions of atoms"""
        return self.pos

    def get_positions(self): return self.get_pos()
     
    def get_origin(self):
        """return origins of distributed multipoles"""
        return self.origin

    def get_origins(self): return self.get_origin()
    
    def get_name(self):
        """return the name of distribution"""
        return self.name
    
    def get_nfrags(self):
        """return number of distributed sites"""
        return self.nfrag

    def get_natoms(self):
        """return number of atoms"""
        return len(self.pos)

    def get_magnitudes(self):
        """
Return the magnitudes of distributed multipoles. Octupoles and hexadecapoles NOT implemented yet!
To obtain magnitude of quadrupole according to Jackson convention just multiply
the magnitudes by 2.0

Warning: The 0 value is returned for octupoles. Hexadecapole magnitude is not returned at all yet!"""
        # charges
        c1 = self[0].copy()
        # dipoles
        c2 = sqrt ( ((self[1].copy())**2).sum(axis=1) )
        # quadrupoles and octupoles
        temp = self.copy()
        if not temp.traceless:
           temp.MAKE_FULL()       
           temp.MakeTraceless()
           temp.makeDMAfromFULL()
        v3 = temp[2]; v4 = temp[3]
        c3 = zeros(self.nfrag, float64)
        c4 = zeros(self.nfrag, float64)
        for i in range(self.nfrag):
            Q11 = v3[i,0] 
            Q22 = v3[i,1]
            Q33 = v3[i,2]
            Q12 = v3[i,3]
            Q13 = v3[i,4]
            Q23 = v3[i,5]
            v = 1./4. * Q33**2 + 1./6. * (Q12**2 + Q13**2 + Q23**2) + 1./24. * (Q11 - Q22)**2
            c3[i] = sqrt(v)
        return c1, c2, c3, c4

    # ---- SET methods

    def set_name(self,name):
        """sets the name for the object""" 
        self.name = name
        return
    
    def set_structure(self, pos=None, origin=None, atoms=None, equal=False):
        """sets new positions or/and origins and atoms"""
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
           if isinstance(atoms[0], Atom): 
              self.atoms = atoms
           else:
              self.atoms = [ Atom(x) for x in atoms.split(',') ]
        # equal the positions and origins
        if equal:
           self.origin = pos.copy()

    def set_moments(self,charges=None,dipoles=None,
                         quadrupoles=None,octupoles=None,
                         hexadecapoles=None):
        """\
set multipoles given as n-d array
where n is rank of multipole moment"""
        if charges       is not None: self.DMA[0] = charges       .copy()
        if dipoles       is not None: self.DMA[1] = dipoles       .copy()
        if quadrupoles   is not None: self.DMA[2] = quadrupoles   .copy()
        if octupoles     is not None: self.DMA[3] = octupoles     .copy()
        if hexadecapoles is not None: self.DMA[4] = hexadecapoles .copy()
        return
   
    # ---- property descriptor methods
 
    def if_traceless(self):
        """Is the dma object traceless?"""
        return self.traceless

    def if_hexadecapoles(self):
        """Does this object contain hexadecapole moments?"""
        return self.has_hexadecapoles

    # --- auxiliary methods
     
    def interchange(self,ind):
        """
premute all moments, positions and origins using ind list 
(ind has the same meaning as in interchange function"""
        n_range = 4
        if self.has_hexadecapoles: n_range = 5
        for i in range(n_range):
            self.DMA[i]  = interchange(self.DMA[i],ind)
        self.pos = interchange(self.pos,ind)
        self.origin = interchange(self.origin,ind)
        return
        
    def copy(self):
        """creates deepcopy of the object"""
        return copy.deepcopy(self)

    def contract(self,contrlist):
        """shrink the dimensions of DMA object by eliminating rows in attributes"""
        K = len(contrlist)
        self.nfrag = K
        self.DMA[0] = self.DMA[0][contrlist]
        self.DMA[1] = self.DMA[1][contrlist]
        self.DMA[2] = self.DMA[2][contrlist]
        self.DMA[3] = self.DMA[3][contrlist]
        if self.has_hexadecapoles:
           self.DMA[4] = self.DMA[4][contrlist]
        self.origin = self.origin[contrlist]
        return
                     
    def write(self, file, type='c'):
        """writes  the DMA distribution in a file or in a XYZ file (structure)"""
        newfile = open(file,'w')
        if type=='c' or type=='dma':
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
        
        elif type=='xyz':
           log = "  %i\n\n"%len(self.pos)
           pos = self.pos*UNITS.BohrToAngstrom
           for i,atom in enumerate(self.atoms):
              log+= " %-10s %14.8f %14.8f %14.8f\n" % (atom.symbol, pos[i,0], pos[i,1], pos[i,2])
           log+='\n'
           
        newfile.write(log)
        newfile.close() 
        return          

    # PRIVATE METHODS

    def __contrListFromUaList(self,ua_list):
       """creates the contraction list from ua_list"""
       return []
            
    def __getitem__(self,index): 
        return array(self.DMA[index])

    def __setitem__(self,index,value):
        self.DMA[index] = array(value)
    
    def __len__(self):
        """Return the number of distributed sites"""
        n_range = 4
        if self.has_hexadecapoles: n_range = 5
        length = len(self.DMA[0])
        for i in range(n_range-1):
            assert length == len(self.DMA[i]), "The number of distributed %i-th poles is not equal to number of 0-poles (charges)!" % (i+1)
        return(length)

    ### --- mathematical operations
    ###     Warning! If between 2 DMA objects they have to have
    ###     same dimensions!
    
    def __add__(self,other):
        if (isinstance(self,DMA) and (not isinstance(other,DMA))):
           case = 1
           q = self.DMA[0] + other
           m = self.DMA[1] + other
           T = self.DMA[2] + other
           O = self.DMA[3] + other
           if self.has_hexadecapoles: 
              H = self.DMA[4] + other
              case =-1
        elif ((not isinstance(self,DMA)) and isinstance(other,DMA)):
           case = 2
           q = other.DMA[0] + self
           m = other.DMA[1] + self
           T = other.DMA[2] + self
           O = other.DMA[3] + self
           if other.has_hexadecapoles: 
              H = other.DMA[4] + self
              case = -2
        else:
           case = 3
           q = other.DMA[0] + self.DMA[0]
           m = other.DMA[1] + self.DMA[1]
           T = other.DMA[2] + self.DMA[2]
           O = other.DMA[3] + self.DMA[3]     
           if self.has_hexadecapoles and other.has_hexadecapoles: 
              H = other.DMA[4] + self.DMA[4]
              case = -3

        # return DMA without hexadecapoles
        if case > 0:
           return DMA(q=q,m=m,T=T,O=O,
                      pos=self.pos.copy(),origin=self.origin.copy())

        # return DMA with hexadecapoles
        else:
           return DMA(q=q,m=m,T=T,O=O,H=H,
                      pos=self.pos.copy(),origin=self.origin.copy())

    def __sub__(self,other):
        q = self.DMA[0] - other.DMA[0]
        m = self.DMA[1] - other.DMA[1]
        T = self.DMA[2] - other.DMA[2]
        O = self.DMA[3] - other.DMA[3]
        if has.hexadecapoles:
           H = self.DMA[4] - other.DMA[4]
           return DMA(q=q,m=m,T=T,O=O,H=H,
                      pos=self.pos.copy(),origin=self.origin.copy())
        else:
           return DMA(q=q,m=m,T=T,O=O,
                      pos=self.pos.copy(),origin=self.origin.copy())

    
    def __div__(self,other):
        if (isinstance(self,DMA) and (not isinstance(other,DMA))):
           case = 1
           q = self.DMA[0] / other
           m = self.DMA[1] / other
           T = self.DMA[2] / other
           O = self.DMA[3] / other
           if self.has_hexadecapoles:
              H = self.DMA[4] / other 
              case = -1
        elif ((not isinstance(self,DMA)) and isinstance(other,DMA)):
           q = other.DMA[0] / self
           m = other.DMA[1] / self
           T = other.DMA[2] / self
           O = other.DMA[3] / self
           if other.has_hexadecapoles:
              H = other.DMA[4] / self
              case = -2
        else:
           q = self.DMA[0] / other.DMA[0]
           m = self.DMA[1] / other.DMA[1]
           T = self.DMA[2] / other.DMA[2]
           O = self.DMA[3] / other.DMA[3]     
           if self.has_hexadecapoles and other.has_hexadecapoles:
              H = self.DMA[4] / other.DMA[4]
              case = -3

        if case > 0:
           return DMA(q=q,m=m,T=T,O=O,
                      pos=self.pos.copy(),origin=self.origin.copy())
        else:
           return DMA(q=q,m=m,T=T,O=O,H=H,
                      pos=self.pos.copy(),origin=self.origin.copy())


    def __mul__(self,other):
        if (isinstance(self,DMA) and (not isinstance(other,DMA))):
           case = 1
           q = self.DMA[0] * other
           m = self.DMA[1] * other
           T = self.DMA[2] * other
           O = self.DMA[3] * other
           if self.has_hexadecapoles:
              H = self.DMA[4] * other
              case = -1
        elif ((not isinstance(self,DMA)) and isinstance(other,DMA)):
           case = 2
           q = other.DMA[0] * self
           m = other.DMA[1] * self
           T = other.DMA[2] * self
           O = other.DMA[3] * self
           if other.has_hexadecapoles:
              H = other.DMA[4] * self
              case = -2
        else:
           case = 3
           q = self.DMA[0] * other.DMA[0]
           m = self.DMA[1] * other.DMA[1]
           T = self.DMA[2] * other.DMA[2]
           O = self.DMA[3] * other.DMA[3]            
           if self.has_hexadecapoles and other.has_hexadecapoles:
              H = self.DMA[4] * other.DMA[4]
              case = -3

        if case > 0:
           return DMA(q=q,m=m,T=T,O=O,
                      pos=self.pos.copy(),origin=self.origin.copy())
        else:
           return DMA(q=q,m=m,T=T,O=O,H=H,
                      pos=self.pos.copy(),origin=self.origin.copy())

    def __rmul__(self,other):
        if (isinstance(self,DMA) and (not isinstance(other,DMA))):
           case = 1
           q = self.DMA[0] * other
           m = self.DMA[1] * other
           T = self.DMA[2] * other
           O = self.DMA[3] * other
           if self.has_hexadecapoles:
              H = self.DMA[4] * other
              case = -1
        elif ((not isinstance(self,DMA)) and isinstance(other,DMA)):
           case = 2
           q = other.DMA[0] * self
           m = other.DMA[1] * self
           T = other.DMA[2] * self
           O = other.DMA[3] * self
           if other.has_hexadecapoles:
              H = other.DMA[4] * self
              case = -2

        if case > 0:
           return DMA(q=q,m=m,T=T,O=O,
                      pos=self.pos.copy(),origin=self.origin.copy())
        else:
           return DMA(q=q,m=m,T=T,O=O,H=H,
                      pos=self.pos.copy(),origin=self.origin.copy())

    def __neg__(self):
        q = -self.DMA[0]
        m = -self.DMA[1]
        T = -self.DMA[2]
        O = -self.DMA[3]
        if self.has_hexadecapoles:
          H = -self.DMA[4]
          return DMA(q=q,m=m,T=T,O=O,H=H,
                     pos=self.pos.copy(),origin=self.origin.copy())
        else:
          return DMA(q=q,m=m,T=T,O=O,
                     pos=self.pos.copy(),origin=self.origin.copy())


    def MAKE_FULL(self):
        """creates arrays of ordinary forms of distributed multipoles
        in full format."""
        CHARGES = array(self.DMA[0])
        DIPOLES = array(self.DMA[1])
        QDPOLES = zeros((self.nfrag,3,3),dtype=float64)
        OCTPLES = zeros((self.nfrag,3,3,3),dtype=float64)
        if self.has_hexadecapoles:
           HEXPLES = zeros((self.nfrag,3,3,3,3),dtype=float64)
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
            
            # add octupoles
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
            
            # add hexadecapoles
            if self.has_hexadecapoles:
               h = array(self.DMA[4])    # REDUCED FORMAT COMPONENTS

               a = h[site,0]             # XXXX
               HEXPLES[site,0,0,0,0] = a
               
               a = h[site,1]             # YYYY
               HEXPLES[site,1,1,1,1] = a
               
               a = h[site,2]             # ZZZZ
               HEXPLES[site,2,2,2,2] = a
               
               a = h[site,3]             # XXXY
               HEXPLES[site,0,0,0,1] = a
               HEXPLES[site,0,0,1,0] = a
               HEXPLES[site,0,1,0,0] = a
               HEXPLES[site,1,0,0,0] = a
               
               a = h[site,4]             # XXXZ
               HEXPLES[site,0,0,0,2] = a
               HEXPLES[site,0,0,2,0] = a
               HEXPLES[site,0,2,0,0] = a
               HEXPLES[site,2,0,0,0] = a
               
               a = h[site,5]             # YYYX
               HEXPLES[site,0,1,1,1] = a
               HEXPLES[site,1,0,1,1] = a
               HEXPLES[site,1,1,0,1] = a
               HEXPLES[site,1,1,1,0] = a
               
               a = h[site,6]             # YYYZ
               HEXPLES[site,1,1,1,2] = a
               HEXPLES[site,1,1,2,1] = a
               HEXPLES[site,1,2,1,1] = a
               HEXPLES[site,2,1,1,1] = a
               
               a = h[site,7]             # ZZZX
               HEXPLES[site,0,2,2,2] = a
               HEXPLES[site,2,0,2,2] = a
               HEXPLES[site,2,2,0,2] = a
               HEXPLES[site,2,2,2,0] = a
               
               a = h[site,8]             # ZZZY
               HEXPLES[site,1,2,2,2] = a
               HEXPLES[site,2,1,2,2] = a
               HEXPLES[site,2,2,1,2] = a
               HEXPLES[site,2,2,2,1] = a
               
               a = h[site,9]             # XXYY
               HEXPLES[site,0,0,1,1] = a
               HEXPLES[site,0,1,0,1] = a
               HEXPLES[site,0,1,1,0] = a
               HEXPLES[site,1,0,0,1] = a
               HEXPLES[site,1,0,1,0] = a
               HEXPLES[site,1,1,0,0] = a
               
               a = h[site,10]            # XXZZ
               HEXPLES[site,0,0,2,2] = a
               HEXPLES[site,0,2,0,2] = a
               HEXPLES[site,0,2,2,0] = a
               HEXPLES[site,2,0,0,2] = a
               HEXPLES[site,2,0,2,0] = a
               HEXPLES[site,2,2,0,0] = a
               
               a = h[site,11]            # YYZZ
               HEXPLES[site,1,1,2,2] = a
               HEXPLES[site,1,2,1,2] = a
               HEXPLES[site,1,2,2,1] = a
               HEXPLES[site,2,1,1,2] = a
               HEXPLES[site,2,1,2,1] = a
               HEXPLES[site,2,2,1,1] = a
               
               a = h[site,12]            # XXYZ
               HEXPLES[site,0,0,1,2] = a
               HEXPLES[site,0,0,2,1] = a
               HEXPLES[site,0,1,0,2] = a
               HEXPLES[site,0,1,2,0] = a
               HEXPLES[site,0,2,0,1] = a
               HEXPLES[site,0,2,1,0] = a
               HEXPLES[site,1,0,0,2] = a
               HEXPLES[site,1,0,2,0] = a
               HEXPLES[site,1,2,0,0] = a
               HEXPLES[site,2,0,0,1] = a
               HEXPLES[site,2,0,1,0] = a
               HEXPLES[site,2,1,0,0] = a
               
               a = h[site,13]            # YYXZ
               HEXPLES[site,0,1,1,2] = a
               HEXPLES[site,0,1,2,1] = a
               HEXPLES[site,0,2,1,1] = a
               HEXPLES[site,1,0,1,2] = a
               HEXPLES[site,1,0,2,1] = a
               HEXPLES[site,1,1,0,2] = a
               HEXPLES[site,1,1,2,0] = a
               HEXPLES[site,1,2,0,1] = a
               HEXPLES[site,1,2,1,0] = a
               HEXPLES[site,2,0,1,1] = a
               HEXPLES[site,2,1,0,1] = a
               HEXPLES[site,2,1,1,0] = a
               
               a = h[site,14]            # ZZXY
               HEXPLES[site,0,1,2,2] = a
               HEXPLES[site,0,2,1,2] = a
               HEXPLES[site,0,2,2,1] = a
               HEXPLES[site,1,0,2,2] = a
               HEXPLES[site,1,2,0,2] = a
               HEXPLES[site,1,2,2,0] = a
               HEXPLES[site,2,0,1,2] = a
               HEXPLES[site,2,0,2,1] = a
               HEXPLES[site,2,1,0,2] = a
               HEXPLES[site,2,1,2,0] = a
               HEXPLES[site,2,2,0,1] = a
               HEXPLES[site,2,2,1,0] = a
        
        if self.has_hexadecapoles:    
           self.DMA_FULL = [ self.origin , CHARGES , DIPOLES , QDPOLES , OCTPLES , HEXPLES ]
        else:
           self.DMA_FULL = [ self.origin , CHARGES , DIPOLES , QDPOLES , OCTPLES ]
        self.full = True

    def makeDMAfromFULL(self):
        """Saves the reduced-formatted DMA data from full-formatted data"""
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
        if self.has_hexadecapoles:
           self.DMA[4][:, 0]  = array(self.DMA_FULL[5][:,0,0,0,0])  
           self.DMA[4][:, 1]  = array(self.DMA_FULL[5][:,1,1,1,1])
           self.DMA[4][:, 2]  = array(self.DMA_FULL[5][:,2,2,2,2])
           self.DMA[4][:, 3]  = array(self.DMA_FULL[5][:,0,0,0,1])
           self.DMA[4][:, 4]  = array(self.DMA_FULL[5][:,0,0,0,2])
           self.DMA[4][:, 5]  = array(self.DMA_FULL[5][:,1,1,1,0])
           self.DMA[4][:, 6]  = array(self.DMA_FULL[5][:,1,1,1,2])
           self.DMA[4][:, 7]  = array(self.DMA_FULL[5][:,2,2,2,0])
           self.DMA[4][:, 8]  = array(self.DMA_FULL[5][:,2,2,2,1])
           self.DMA[4][:, 9]  = array(self.DMA_FULL[5][:,0,0,1,1])
           self.DMA[4][:,10]  = array(self.DMA_FULL[5][:,0,0,2,2])
           self.DMA[4][:,11]  = array(self.DMA_FULL[5][:,1,1,2,2])
           self.DMA[4][:,12]  = array(self.DMA_FULL[5][:,0,0,1,2])
           self.DMA[4][:,13]  = array(self.DMA_FULL[5][:,1,1,0,2])
           self.DMA[4][:,14]  = array(self.DMA_FULL[5][:,2,2,0,1])
 
        del self.DMA_FULL
        self.full = False
       
    def get_const(self, origin=array([0.,0.,0.]) ):
        """
 Calculate the invariants of total multipole moments
 Returns the tuple of four invariant sets for charge, dipole, quadrupole
 and octupole moments:

 1) charge invariant: 
     total charge

 2) dipole invariant:
     norm of a dipole vector

 3) quadrupole invariants:
     trace(Q)
     0.5 * (trace(Q)^2 - trace(Q^2))
     det(Q)
 
     where Q is quadrupole moment (primitive form)

 4) octupole invariants:
     F1 = E_iik * E_ppk
     F2 = E_ijj * E_iqq
     F3 = E_ijk * E_ijk
     F4 = E_ijk * E_kij
     F5 = 0.5 * (e_kp * E_ijk * E_pji + e_ri * E_ijk * E_kjr)

     where E is octupole moment (primitive form)
     and e is 2-order antisymmetric (permutation) tensor

 Notes:
  1) invariants or 3-rd order Cartesian tensors:
        from. F. Ahmad, Arch. Mech. 63, 4, pp383-392 Warszawa 2011
        formulae on p.390 for F1-F5

  2) no invariants are computed for hexadecapole moments!
     """
        eij = zeros((3, 3), float64)
        eij[0,1] = eij[1,2] = eij[2,0] =  1.0
        eij[1,0] = eij[2,1] = eij[0,2] = -1.0
        
        eijk = zeros((3, 3, 3), float64)
        eijk[0, 1, 2] = eijk[1, 2, 0] = eijk[2, 0, 1] =  1.0
        eijk[0, 2, 1] = eijk[2, 1, 0] = eijk[1, 0, 2] = -1.0
        tot = self.get_mult(origin)
        tot.MAKE_FULL()
        # invariants for dipole moment
        mu  = sqrt(sum(tot.DMA_FULL[2][0]**2))
        # invariants for quadrupole moment
        tr = tot.DMA_FULL[3][0].trace()
        tr2= (dot(tot.DMA_FULL[3][0],tot.DMA_FULL[3][0])).trace()
        qad_1 = tr
        qad_2 = 0.500*(tr*tr-tr2)
        qad_3 = linalg.det(tot.DMA_FULL[3][0])
        # invariants for octupole moment
        # from. F. Ahmad, Arch. Mech. 63, 4, pp383-392 Warszawa 2011
        # formulae from p.390 for F1-F5
        oct = tot.DMA_FULL[4][0]
        oct_1 = dot(oct.trace(axis1=0,axis2=1),oct.trace(axis1=0,axis2=1))
        oct_2 = dot(oct.trace(axis1=1,axis2=2),oct.trace(axis1=1,axis2=2))
        oct_3 = tensordot(oct,oct,((0,1,2),(0,1,2)))
        oct_4 = tensordot(oct,oct,((0,1,2),(1,2,0)))
        A_kp = tensordot(oct,oct,((0,1),(2,1)))
        B_ri = tensordot(oct,oct,((1,2),(1,0)))
        oct_5 = 0.5*(sum(eij*A_kp,axis=None) + sum(eij.transpose()*B_ri,axis=None))
        #oct = 1./6. * (sum(tot.DMA_FULL[4][0]*eijk)) * eijk
        qad = array([qad_1, qad_2, qad_3],float64)
        oct = array([oct_1,oct_2,oct_3,oct_4,oct_5],float64)
        return tot.DMA_FULL[1][0], mu, qad, oct
 
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

           # Hexadecapole moment                                                                                       
           if self.has_hexadecapoles:
              for H in self.DMA_FULL[5]:
                  W = H.copy()
                  H *= (35./8.)
                  Wt = W.trace()
                  Wtt= W.trace().trace()
                  c = (5./8.)
                  for i in [0,1,2]:
                      for j in [0,1,2]:
                          for k in [0,1,2]:
                              for l in [0,1,2]:
                                  if i==j==k==l:                 # AAAA
                                     H -= (15./4.) * Wt[i,i]       - 3.0 * Wtt
                                  elif i==j and k==l and i!=k:   # AABB
                                     H -= c * (Wt[i,i] + Wt[k,k])  - Wtt
                                  elif i==k and j==l and i!=j:   # ABAB
                                     H -= c * (Wt[i,i] + Wt[j,j])  - Wtt
                                  elif i==l and j==k and i!=j:   # ABBA
                                     H -= c * (Wt[i,i] + Wt[j,j])  - Wtt
                                  elif i==j and j!=k and k!=l:   # AABC
                                     H -= c * Wt[k,l]
                                  elif i!=j and j!=k and k==l:   # ABCC
                                     H -= c * Wt[i,j]
                                  elif i!=j and j==k and k!=l:   # ABBC
                                     H -= c * Wt[i,l]
                                  elif i!=j and i==k and k!=l:   # ABAC
                                     H -= c * Wt[j,l]
                                  elif i!=j and j!=k and j==l:   # ABCB
                                     H -= c * Wt[i,k]
                                  elif i==l and i!=j and j!=k:   # ABCA
                                     H -= c * Wt[j,k]
                            
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
           
           # Fourth moments
           if self.has_hexadecapoles:
              M4 = array(self.DMA_FULL[5])   # hexadecapoles centered at old origins 
              #raise NotImplementedError, 'Recentering hexadecapoles is not implemented yet! Quitting...'
              # first find the dipoles, quadrupoles and octupoles at the origin [0,0,0]
              other = self.copy()
              other.has_hexadecapoles = False
              other.MAKE_FULL()
              other.ChangeOrigin(zero=True)
              other.MAKE_FULL()
              M10 = other.DMA_FULL[2]; M20 = other.DMA_FULL[3]; M30 = other.DMA_FULL[4]
              del other
              # next define the permutation sequences for generalized outer products [ M_k x D^(4-k) ] where k is rank of multipole at the origin [0,0,0]
              T1 = [(0,1,2,3,4),(0,2,1,3,4),(0,2,3,1,4),(0,2,3,4,1)]                          # dipole permutations
              T2 = [(0,1,2,3,4),(0,1,3,2,4),(0,3,1,4,2),(0,3,4,1,2),(0,1,3,4,2),(0,3,1,2,4)]  # quadrupole permutations
              T3 = [(0,1,2,3,4),(0,1,2,4,3),(0,1,4,2,3),(0,4,1,2,3)]                          # octupole permutations
              # then construct the basis generalized outer tensor products
              D2 = B2 - A2
              D3 = B3 - A3
              D4 = array_outer_product_2_2(B2,B2) - array_outer_product_2_2(A2,A2)
              Y1 = array_outer_product_1_3(M10,D3)
              Y2 = array_outer_product_2_2(M20,D2)
              Y3 = array_outer_product_3_1(M30,D )
              # finally, compute the new hexadecapoles from zeroth moments
              new_M4 = M4 + array_outer_product_1_n(M0,D4)
              for perm in T1: new_M4 -= Y1.transpose(perm)
              for perm in T2: new_M4 += Y2.transpose(perm)
              for perm in T3: new_M4 -= Y3.transpose(perm)


           if self.nfrag == 1:
              self.origin = new_origin_set[0]
           
           self.DMA_FULL.append(new_origin_set)
           self.DMA_FULL[2] = new_M1
           self.DMA_FULL[3] = new_M2
           self.DMA_FULL[4] = new_M3
           if self.has_hexadecapoles:
              self.DMA_FULL[5] = new_M4
           
           self.makeDMAfromFULL()
           self.origin = array(new_origin_set)
           
        else: raise Exception("\nerror: no FULL DMA object created! quitting...\n")

    def MakeUa(self,ua_list,change_origin=True,contract=False):
        """\
Transforms the object to the contracted DMA form employing united atoms. 
Usage:
MakeUa( ua_list )
        
where 'ua_list' is a list of tuples ti containing integers:
        
ua+list = [ t1, t2, t3 , ... ], 
ti = ( A, a1, a2, a3, ... )
        
For i-th tuple corresponding to i-th UA the first integer A
is the atom id of UA center atom and the remaining a1, a2 etc relate to 
the atoms supposed to be contracted within a united atom UA for A-th center.
The numbers are normal numbers (not in Python convention)."""
       
        if self.has_hexadecapoles: rs = range(5)
        else:                      rs = range(4) 
        origin = self.origin
        if change_origin:
           self.MAKE_FULL()
           self.ChangeOrigin(zero=True)
        for t in ua_list:
            for i in t[1:]:
                for c in rs:
                    self.DMA[c][t[0]-1]+= self.DMA[c][i-1]
                    if not c: self.DMA[c][i-1] = 0
                    else:     self.DMA[c][i-1].fill(0)
        if change_origin:
           self.MAKE_FULL()
           self.ChangeOrigin(new_origin_set=origin)
      
        if contract:
           contrlist = self.__contrListFromUaList(ua_list)
           self.contract(contrlist)

    def translate(self,vector):
        """translate origins and positions by a vector"""
        self.origin += vector
        self.pos    += vector
        return
    
    def Rotate(self,rotmat):
        """rotates the ordinary full-formatted DMA_FULL
           in molecular space based on rotation matrix
           given"""

        #self.MAKE_FULL()
        if self.full:
           ### transform the origins and positions
           for i in xrange(len(self.origin)):
               self.origin[i] = dot(self.origin[i],rotmat)
           for i in xrange(len(self.pos)):
               self.pos[i] = dot(self.pos[i],rotmat)
           ### transform the dipoles
           self.DMA_FULL[2] = dot(self.DMA[1],rotmat)

           ### transform the quadrupoles
           quadrupole = zeros((self.nfrag,3,3),dtype=float64)
           for i in range(self.nfrag):
               quadrupole[i] = dot(rotmat.T,dot( self.DMA_FULL[3][i] ,rotmat))

           self.DMA_FULL[3] = quadrupole
        
           ### transform the octupoles
           octupole = zeros((self.nfrag,3,3,3),dtype=float64)
           for i in range(self.nfrag):
               octupole[i] = tensordot(rotmat.T,tensordot(
                                       rotmat.T,tensordot(
                                       rotmat.T,self.DMA_FULL[4][i] , (1,2)
                                     ),(1,2)
                                     ),(1,2)
                                     )

           self.DMA_FULL[4] = octupole

           ### transform the hexadecapoles
           if self.has_hexadecapoles:
              hexadecapole = zeros((self.nfrag,3,3,3,3),dtype=float64)                       
              for i in range(self.nfrag):
                  hexadecapole[i] = tensordot(rotmat.T,tensordot(
                                              rotmat.T,tensordot(
                                              rotmat.T,tensordot(
                                              rotmat.T,self.DMA_FULL[5][i] , (1,2)
                                            ),(1,2)
                                            ),(1,2)
                                            ),(1,2)
                                            )

              self.DMA_FULL[5] = hexadecapole

           self.makeDMAfromFULL()
           
        else: raise Exception("\nerror: no FULL DMA object created! quitting...\n")

    def ChangeUnits(self,charges,dipoles,quadrupoles,octupoles,hexadecapoles=1.0):
        """changes the units"""
        self.DMA[0] *= charges
        self.DMA[1] *= dipoles
        self.DMA[2] *= quadrupoles
        self.DMA[3] *= octupoles
        if self.has_hexadecapoles:
           self.DMA[4] *= hexadecapoles
        return
        
    def get_mult_c(self, origin=zeros(3,dtype=float64), hexadecapoles=False):
        """Calculates overall primitive moments from charges.
This tool has a purpose of testing population analysis obtained by
fitting to the molecular ab initio potential or other methods.

For an example, refer to the construction of total molecular 
solvatochromic multipoles from fitted distributed charges:

H. Lee, J.-H. Choi and M. Cho, J. Chem. Phys. 137(11), 114307 (2012)
"""
        overall = DMA(nfrag=1)
        overall.name = 'Test of reproducing multipoles from charges'
        overall.pos = zeros((1,3),dtype=float64)
        overall.pos[0] = origin
        
        ### make full format of DMA
        tmp = self.copy()
        tmp.MAKE_FULL()
        
        ### compute molecular solvatochromic moments
        mu   = zeros((3),dtype=float64)
        quad = zeros((3,3),dtype=float64)
        oct  = zeros((3,3,3),dtype=float64)
        hex  = zeros((3,3,3,3),dtype=float64)
        for atom in range(self.nfrag):
            r     = tmp.pos[atom] - origin
            qatom = tmp.DMA[0][atom]
            ### calculate dipole moment
            mu   += qatom * r 
            ### calculate quadrupole moment
            quad += qatom * outer (r,r)
            ### calculate octupole moment
            oct  += qatom * outer( r, outer (r,r) ).reshape(3,3,3)
            ### calculate hexadecapole moment
            if hexadecapoles:
               hex  += qatom * outer( r, outer (r, outer(r,r) ) ).reshape(3,3,3,3)

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

        if hexadecapoles:
           overall.DMA[4][0, 0]  = hex[0,0,0,0]  
           overall.DMA[4][0, 1]  = hex[1,1,1,1]
           overall.DMA[4][0, 2]  = hex[2,2,2,2]
           overall.DMA[4][0, 3]  = hex[0,0,0,1]
           overall.DMA[4][0, 4]  = hex[0,0,0,2]
           overall.DMA[4][0, 5]  = hex[1,1,1,0]
           overall.DMA[4][0, 6]  = hex[1,1,1,2]
           overall.DMA[4][0, 7]  = hex[2,2,2,0]
           overall.DMA[4][0, 8]  = hex[2,2,2,1]
           overall.DMA[4][0, 9]  = hex[0,0,1,1]
           overall.DMA[4][0,10]  = hex[0,0,2,2]
           overall.DMA[4][0,11]  = hex[1,1,2,2]
           overall.DMA[4][0,12]  = hex[0,0,1,2]
           overall.DMA[4][0,13]  = hex[1,1,0,2]
           overall.DMA[4][0,14]  = hex[2,2,0,1]

        #overall.DMA[1] *= UNITS.BohrElectronToDebye
        #overall.DMA[2] *= UNITS.BohrElectronToDebye * UNITS.BohrToAngstrom
        #overall.DMA[2] *= UNITS.BohrElectronToDebye * UNITS.BohrToAngstrom**2
          
        return overall
    
    
    def get_mult(self, origin=zeros(3,dtype=float64), hexadecapoles=False):
        """Calculates overall primitive multipole moments from distributed multipoles.
This can be used for a test of the correctness of multipole analysis. The exact distributed
expansions should sum up to the total molecular expansion with respect to origin privided.
"""
        if hexadecapoles and not self.has_hexadecapoles:
           print " WARNING: This DMA object does not contain hexadecapoles but they were requested. It will be ignored by Coulomb."
           hexadecapoles = False
        overall = DMA(nfrag=1)
        #overall.name = 'Test of reproducing multipoles from charges [units: Debyes(*Angstroms^n)]'
        overall.name = 'Test of reproducing overall multipoles from CAMM/CBAMM/LMTP/DMA [units: A.U.]'
        overall.pos = zeros((1,3),dtype=float64)
        overall.pos[0] = origin
        tmp = self.copy()
        tmp.MAKE_FULL()
        tmp.ChangeOrigin(zero=1)
        
        ### make full format of DMA
        tmp.MAKE_FULL()
        
        ### compute molecular solvatochromic moments
        mu   = zeros((3),dtype=float64)
        quad = zeros((3,3),dtype=float64)
        oct  = zeros((3,3,3),dtype=float64)
        hex  = zeros((3,3,3,3),dtype=float64)
        for atom in range(self.nfrag):
            r     = tmp.origin[atom]### zmiana origin z pos!!!
            qatom = tmp.DMA[0][atom]
            ### calculate dipole moment
            mu   += tmp.DMA[1][atom]
            ### calculate quadrupole moment
            quad += tmp.DMA_FULL[3][atom]
            ### calculate octupole moment
            oct  += tmp.DMA_FULL[4][atom]
            ### calculate hexadecapole moment
            if hexadecapoles:
               hex  += tmp.DMA_FULL[5][atom]

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

        if hexadecapoles:
           overall.DMA[4][0, 0]  = hex[0,0,0,0]
           overall.DMA[4][0, 1]  = hex[1,1,1,1]
           overall.DMA[4][0, 2]  = hex[2,2,2,2]
           overall.DMA[4][0, 3]  = hex[0,0,0,1]
           overall.DMA[4][0, 4]  = hex[0,0,0,2]
           overall.DMA[4][0, 5]  = hex[1,1,1,0]
           overall.DMA[4][0, 6]  = hex[1,1,1,2]
           overall.DMA[4][0, 7]  = hex[2,2,2,0]
           overall.DMA[4][0, 8]  = hex[2,2,2,1]
           overall.DMA[4][0, 9]  = hex[0,0,1,1]
           overall.DMA[4][0,10]  = hex[0,0,2,2]
           overall.DMA[4][0,11]  = hex[1,1,2,2]
           overall.DMA[4][0,12]  = hex[0,0,1,2]
           overall.DMA[4][0,13]  = hex[1,1,0,2]
           overall.DMA[4][0,14]  = hex[2,2,0,1]


        overall.MAKE_FULL()
        overall.ChangeOrigin(new_origin_set=array([origin],float64) )
        
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

        if self.has_hexadecapoles:
           log+= "\n"
           log+= " Distributed fourth-order property\n"              
           log+= " ---------------------------------\n"
           log+= "\n"
           hexadecapoles = self.DMA[4]
           log+= "%16s %16s %16s\n" % ('XXXX'.rjust(16),'YYYY'.rjust(16),'ZZZZ'.rjust(16))
           log+= "%16s %16s %16s %16s %16s %16s\n" %\
                ('XXXY'.rjust(16),'XXXZ'.rjust(16),'YYYX'.rjust(16),
                 'YYYZ'.rjust(16),'ZZZX'.rjust(16),'ZZZY'.rjust(16))
           log+= "%16s %16s %16s %16s %16s %16s\n" %\
                ('XXYY'.rjust(16),'XXZZ'.rjust(16),'YYZZ'.rjust(16),
                 'XXYZ'.rjust(16),'YYXZ'.rjust(16),'ZZXY'.rjust(16))
           for xfrag in hexadecapoles:
               log+= "%16.6E %16.6E %16.6E\n"%(xfrag[0],xfrag[1],xfrag[2])
               log+= "%16.6E %16.6E %16.6E %16.6E %16.6E %16.6E\n"%\
                (xfrag[3],xfrag[4],xfrag[5],
                 xfrag[6],xfrag[7],xfrag[8])
               log+= "%16.6E %16.6E %16.6E %16.6E %16.6E %16.6E\n"%\
                (xfrag[9 ],xfrag[10],xfrag[11],
                 xfrag[12],xfrag[13],xfrag[14])

        log = self.__end_multipole_section(log)
        
        return str(log) 

    def __end_multipole_section(self, log):
       log+= " "+"-"*100+"\n"                                                                               
       if self.nfrag == 1:
          if self.traceless: 
             log+= (" traceless form@,   origin: %s\n"  % str(self.origin*0.5291772086)[1:-1]).rjust(100)  
          else: 
             log+= (" primitive form@,   origin: %s \n" % str(self.origin*0.5291772086)[1:-1]).rjust(100) 
       else:
          if self.traceless: 
             log+= " traceless form@\n".rjust(100)  
          else: 
             log+= " primitive form@\n".rjust(100) 
       log+= "\n\n\n"
       return log

