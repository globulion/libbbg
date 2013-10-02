# --------------------------------------------------------------- #
#       GAUSSIAN FILE WITH ANHARMONIC FREQUENCY ANALYSIS          #
# --------------------------------------------------------------- #

__all__=['FREQ',]
__version__ = '2.2.1'

from units     import *
from numpy     import *
from utilities import *
import copy

class FREQ(UNITS):
      """ represents gaussian log file for anharmonic frequency calculations 
          The job has to have 'nosymm' option specified otherwise the normal
          modes and L matrix will probably have incorrect ordering!!! """
      def __init__(self,file=None):
        if file: 
          self.dir    = '.'                        # working directory
          # file with harmonic frequencies
          self.file   = file
          self.Natoms = self.find_Natoms()       
          self.Nmodes = self.Natoms * 3 - 6        
          self.Atoms  = self.find_Atoms()
          # total number of cartesian coordinates
          self.a      = self.Natoms * 3            
          self.freq   = self.HarmonicFrequencies() # helico
          #
          self.redmass= self.ReducedMasses()       # clio

          self.redmass= self.redmass[::-1]         # helico
          #
          self.L      = self.Trans()               # clio
          self.L_     = self.Trans() 
          self.L      = self.Weight()
          self.L      = self.L[:,::-1]             # helico
          self.L_     = self.L_[:,::-1]
          self.dipole = self.Dipole(file)

          self.K3 = self.K_34()                    # helico [cm-1]
          self._w = False                          # if self.w()

      def copy(self):
          """return deep copy of self"""
          return copy.deepcopy(self)

      def if_w(self):
          """returns the answer whether the object was mass-multiplied"""
          return self._w

      def w(self):
          """return the copy of the original object containing mass-multiplied gijk and L vectors"""
          assert not self._w, 'already mass-multiplied!'
          other = self.copy()
          # L-vectors
          temp = sqrt(self.redmass*self.AmuToElectronMass)[newaxis,:]
          other.L = temp * self.L
          # reduced masses
          other.redmass = self.redmass*self.AmuToElectronMass
          # cubic anharmonic constants
          temp = sqrt(self.redmass)[:,newaxis,newaxis,]
          gijj = temp * self.K3
          temp = sqrt(self.redmass)[newaxis,:,newaxis,]
          gijj = temp * gijj
          temp = sqrt(self.redmass)[newaxis,newaxis,:,]
          gijj = temp * gijj
          other.K3 = gijj
          other._w = True
          return other

      def find_Natoms(self):
          """search for number of atoms"""
          querry = " NAtoms= "
          data = open(self.file)
          line = data.readline()
          while 1:
                 if querry in line: break
                 line = data.readline()
          return int( line.split()[1] )

      def find_Atoms(self):
          """finds atoms and store them in a list"""
          querry = " Symbolic Z-matrix"
          data = open(self.file)
          line = data.readline()
          while 1:
                if querry in line: break
                line = data.readline()
          line = data.readline()
          line = data.readline()
          atoms = []
          for i in range(self.Natoms):
              atoms.append( Atom(line.split()[0]) )
              line = data.readline()
          #for i in atoms: print i
          return atoms
          
      def HarmonicFrequencies(self):
          """derive freqs in helico (anharmonic file)"""
          querry = " Fundamental Bands"
          n = self.Nmodes
          data = open(self.file)
          line = data.readline()

          while 1: 
                if querry in line: break 
                line = data.readline()
          line = data.readline()

          freq = []
          for i in range(n):
              freq.append( float64(line.split()[1]) )
              line = data.readline()

          return array(freq)
      
      def ReducedMasses(self):
          """withdraw reduced masses from Gaussian calculations 
             (anharmonic file) - units: [AMU]"""
          querry = " Reduced masses --- "
          n = self.Nmodes
          data = open(self.file)
          line = data.readline()
          while 1: 
                if querry in line: break 
                line = data.readline()
          
          T = zeros(n)
          for j in range( n/5+bool(n%5) ):
              T[(j*5):j*5+self.dupa(j)] =\
              [ float64(line.replace('D','E').split()[-self.dupa(j):][x])\
                                              for x in range(self.dupa(j)) ]
              for h in range(7+self.a): line = data.readline()

          return array(T,dtype=float64)

      def Trans(self):
          """withdraw transformation matrix"""
          querry = " Coord Atom Element:"
          n = self.Nmodes
          data = open(self.file)
          line = data.readline()

          while 1:
                if querry in line: break
                line = data.readline()
          line = data.readline()
          
          T = zeros((self.a,n))
          for j in range( n/5+bool(n%5) ):
              for i in range(self.a):
                  T[i][(j*5):j*5+self.dupa(j)] =\
                  [ float64(line.replace('D','E').split()[-self.dupa(j):][x])\
                                                  for x in range(self.dupa(j)) ]
                  if (i+1)==self.a:
                     for h in range(8): line = data.readline()
                  else: line = data.readline()
   
          #print "macierz L:"
          #from utilities import PRINTL
          #print PRINTL(array(T,dtype=float64))
          return array(T,dtype=float64)

      def dupa(self,j):
          """some strange but extremely helpful utility:D"""
          if self.Nmodes-j*5 >= 5: return 5
          else                   : return self.Nmodes%5
 
 
      def Weight(self):
          """weight L matrix"""
          querry = " Coord Atom Element:"
          data = open(self.file)
          line = data.readline()

          while 1:
                if querry in line: break
                line = data.readline()
          line = data.readline()
          T = zeros(self.a,dtype=int) # matrix with atom list
          for i in range(self.a):
              T[i] = int( line.split()[1] ) - 1
              line = data.readline()
          
          #for o in self.Atoms : print o
          # weighted L-matrix
          L = zeros((self.a,self.Nmodes))
          for i in range(len(L)):
              for j in range(len(L[0])):
                  sumcia = 0
                  for k in range(len(L)):
                      sumcia += self.Atoms[T[k]].mass * self.AmuToElectronMass * self.L[k][j]**2
                  
                  L[i][j] = self.L[i][j]/ sqrt(sumcia)

          return L

      def Weight_test(self):
          """weight L matrix"""
          querry = " Coord Atom Element:"
          data = open(self.file)
          line = data.readline()
          rrr = array(self.Trans())
          while 1:
                if querry in line: break
                line = data.readline()
          line = data.readline()
          T = zeros(self.a,dtype=int) # matrix with atom list
          for i in range(self.a):
              T[i] = int( line.split()[1] ) - 1
              line = data.readline()
          # weighted L-matrix
          L = zeros((self.a,self.Nmodes))
          for i in range(len(L)):
              for j in range(len(L[0])): 
                  sumcia = 0
                  for k in range(len(L)):
                      sumcia += self.Atoms[T[k]].mass * self.AmuToElectronMass * self.L[k][j]**2
                  L[i][j] = rrr[i][j] /sqrt(sumcia)#/ sqrt(self.mass[T[i]] * self.AmuToElectronMass )

          return L

      def COE(self,structure=[]):
          """calculates the center of squared eigenvector for each mode given the structure [in Bohr]"""
          COEs = []
          for mode in xrange(self.Nmodes):
              vec = self.L[:,mode].reshape(self.Natoms,3)
              r_origin = zeros(3,dtype=float64)       
              for atom in xrange(self.Natoms):
                  r_origin+= sum(vec[atom]**2) * structure[atom] / sum(vec**2) #/ sum(vec**2,axis=0)
              COEs.append(r_origin)
          COEs = array(COEs)
          print " COE for each mode [in Angstrom]\n"
          for i in range(self.Nmodes):
              print (i+1), " mode      ",(COEs[i] * self.BohrToAngstrom)
              
          return array(COEs)
          
      def Dipole(self,file):
          """withdraw dipole moment in cartesian coord in AU!"""

          querry = " Dipole moment (field-independent basis, Debye):"
          data = open(file)
          line = data.readline()

          while 1:
                if querry in line: break
                line = data.readline()
          line = data.readline()
          T = zeros(3,dtype=float64)

          T[0] = line.split()[1]
          T[1] = line.split()[3]
          T[2] = line.split()[5]

          return array(T,dtype=float64) * self.DebyeToBohrElectron
        
      def Quadrupole(self,file):
          """withdraw dipole moment in cartesian coord in AU!"""

          querry = " Quadrupole moment (field-independent basis, Debye-Ang)"
          data = open(file)
          line = data.readline()

          while 1:
                if querry in line: break
                line = data.readline()
          line = data.readline()
          T = zeros(6,dtype=float64)

          T[0] = line.split()[1]
          T[1] = line.split()[3]
          T[2] = line.split()[5]

          line = data.readline()
          
          T[3] = line.split()[1]
          T[4] = line.split()[3]
          T[5] = line.split()[5]          

          return array(T,dtype=float64) * self.DebyeToBohrElectron * self.AngstromToBohr
        
        
      def Octupole(self,file):
          """withdraw octupole moment in cartesian coord in AU!"""

          querry = " Octapole moment (field-independent basis, Debye-Ang**2)"
          data = open(file)
          line = data.readline()

          while 1:
                if querry in line: break
                line = data.readline()
          line = data.readline()
          T = zeros(10,dtype=float64)

          T[0] = line.split()[1]
          T[1] = line.split()[3]
          T[2] = line.split()[5]
          T[5] = line.split()[7]
          
          line = data.readline()
          
          T[3] = line.split()[1]
          T[4] = line.split()[3]
          T[7] = line.split()[5]          
          T[8] = line.split()[7]
          
          line = data.readline()
          
          T[6] = line.split()[1]          
          T[9] = line.split()[3]          
          
          return array(T,dtype=float64) * self.DebyeToBohrElectron * self.AngstromToBohr**2

      def Polarizability(self,file):
          """withdraw polarizability in AU"""

          querry = " Polarizability="
          data = open(file)
          line = data.readline()

          while 1:
                if querry in line: break
                line = data.readline()

          T = zeros((3,3),dtype=float64)
          T[0,0] = line.replace('D','E')[16:31]
          T[0,1] = line.replace('D','E')[31:46]
          T[1,1] = line.replace('D','E')[46:61]
          line = data.readline()
          T[0,2] = line.replace('D','E')[16:31]
          T[1,2] = line.replace('D','E')[31:46]
          T[2,2] = line.replace('D','E')[46:61]
          line = data.readline()
          T[1,0] = T[0,1]
          T[2,0] = T[0,2]
          T[2,1] = T[1,2]

          return array(T,dtype=float64)
        
      # ----------------------------------------------------------------------
      def Dipole_old(self,file):
          """withdraw dipole moment in cartesian coord in AU!. These derivatives 
             are from gaussian log file and hase strange units, unknown to me as for the
             time being! (probably natural units)"""

          querry = " Dipole        ="
          data = open(file)
          line = data.readline()

          while 1:
                if querry in line: break
                line = data.readline()

          T = zeros(3,dtype=float64)
          T[0] = line.replace('D','E')[16:31]
          T[1] = line.replace('D','E')[31:46]
          T[2] = line.replace('D','E')[46:61]

          return array(T,dtype=float64)
        
      def DipoleDeriv(self,file):
          """withdraw first derivs of dipole moment in cartesian coord. These derivatives 
             are from gaussian log file and hase strange units, unknown to me as for the
             time being! (probably natural units)"""

          querry = " DipoleDeriv"
          data = open(file)
          line = data.readline()

          while 1:
                if querry in line: break
                line = data.readline()

          T = zeros((self.a,3),dtype=float64)
          for i in range(self.a):
              T[i][0] = line.replace('D','E')[16:31]
              T[i][1] = line.replace('D','E')[31:46]
              T[i][2] = line.replace('D','E')[46:61]
              line = data.readline()

          return array(T,dtype=float64)

      def FDeriv(self,Print=1,Debye=0,divide=1):
          """ transforms derivatives of dipole moment wrt cartesians to wrt normal modes """
          C = zeros((self.Nmodes,3),dtype=float64)
          LT = transpose( self.Weight() )
          for i in [0,1,2]:
              C[:,i] = dot( LT, self.DipoleDeriv(self.file)[:,i]  )

          # divide by frequencies
          if divide:
             for x in [0,1,2]:
                 for A in range(self.Nmodes):
                     C[A][x] /= sqrt(self.freq[A]  * self.CmRecToHartree )

          if Debye: 
             for x in [0,1,2]:
                 for A in range(self.Nmodes):
                     C[A][x] /= self.BohrElectronToDebye 


          if Print:
             y = arange(self.Nmodes)+1
             print " \n First Derivatives of Dipole Moment \n"
             if Debye:
                       print "                deriv units: Debye "
             else:     print "                deriv inits: AU    "
             print           "                frequencies given in [cm-1]"
             print
             PRINTV(transpose(C),y,self.freq,["x","y","z"])


          return C # C_Aa

      def Intens_1(self,Print,Debye):
          """harmonic intensities of fundamental bands"""
          F = self.FDeriv(0,Debye,1)
          ulazulahyta = zeros((self.Nmodes,3),dtype=float64)
          for cart in [0,1,2]:
              for mode in range(self.Nmodes):
                  ulazulahyta[mode][cart] = F[mode][cart]**2 * sqrt(2*pi) * 1./2. * self.freq[mode] * self.BohrElectronToDebye**2

          intens = sum(ulazulahyta,axis=1)

          if Print:
             y = arange(self.Nmodes)+1
             print " \n Fundamental Harmonic Intensities \n"
             if Debye:
                       print "                deriv units: Debye "
             else:     print "                deriv units: AU    "
             print
             print           "                frequencies [cm-1]"
             PRINT(self.freq)
             print           "                intensities "
             PRINT(intens)

          return intens



      def K_34(self):
          """ derives cubic and quartic anharmonic constants.
              It takes reduced values [cm-1] from gaussian output
              and transforms it to :
                * cubic   force constants : [Hartree*amu(-3/2)*Bohr(-3)]
                * quartic force constants : [Hartree*amu(-2  )*Bohr(-4)]
          """
          n = self.Nmodes
          K3 = zeros((n,n,n)  ,dtype=float64)
          #K4 = zeros((n,n,n,n),dtype=float64)
          data = open(self.file)

          querry = "CUBIC FORCE CONSTANTS IN NORMAL MODES"
          line = data.readline()
          while 1: 
                if querry in line: break 
                line = data.readline()
          for i in range(9): line = data.readline()

          while line.split() != []:
                list = line.split()
                #print list
                i = int(list[0]) - 1
                j = int(list[1]) - 1
                k = int(list[2]) - 1
                d = float64(list[3])

                # transform to [Hartree*amu(-3/2)*Bohr(-3)]
                d *= sqrt(self.freq[i] * self.freq[j] * self.freq[k])
                d *= self.BohrToAngstrom**3
                d /= self.ToRedCubForceConst * self.HartreeToAttoJoule
                # 
                K3[i,j,k]    = d
                K3[i,k,j]    = d
                K3[k,j,i]    = d
                K3[k,i,j]    = d
                K3[j,i,k]    = d
                K3[j,k,i]    = d
                line = data.readline()

          #querry = "QUARTIC FORCE CONSTANTS IN NORMAL MODES"
          #line = data.readline()
          #while 1: 
          #      if querry in line: break 
          #      line = data.readline()
          #for i in range(9): line = data.readline()

          #while line.split() != []:
          #      list = line.split()
          #      i = int(list[0]) - 1
          #      j = int(list[1]) - 1
          #      k = int(list[2]) - 1
          #      l = int(list[3]) - 1
          #      d = float64(list[6])

          #      K4[i][j][k][l] = d
          #      K4[i][j][l][k] = d
          #      K4[i][k][j][l] = d
          #      K4[i][k][l][j] = d
          #      K4[i][l][j][k] = d
          #      K4[i][l][k][j] = d

          #      K4[j][k][l][i] = d
          #      K4[j][k][i][l] = d
          #      K4[j][l][k][i] = d
          #      K4[j][l][i][k] = d
          #      K4[j][i][k][l] = d
          #      K4[j][i][l][k] = d

          #      K4[k][l][i][j] = d
          #      K4[k][l][j][i] = d
          #      K4[k][i][l][j] = d
          #      K4[k][i][j][l] = d
          #      K4[k][j][l][i] = d
          #      K4[k][j][i][l] = d

          #      K4[l][i][j][k] = d
          #      K4[l][i][k][j] = d
          #      K4[l][j][i][k] = d
          #      K4[l][j][k][i] = d
          #      K4[l][k][i][j] = d
          #      K4[l][k][j][i] = d

          #      line = data.readline()


          return K3#,K4
 
