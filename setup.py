#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Standard library for BBG packages 
"""
# -----------------------------------------------
import os, commands
from numpy.distutils.core import setup, Extension
# -----------------------------------------------

# --- Extension modules

# commands to compile and link
f2py_cmd =    'f2py -h gentcf.pyf -m gentcf gentcf.f --overwrite-signature'
f2py_cmd+= '&& f2py -c --fcompiler=gnu95 -m gentcf gentcf.f' 

print "Compiling extension module: ", f2py_cmd
#failure, output = commands.getstatusoutput(f2py_cmd)

# extension module specifications
GENTCF = Extension(name='gentcf',
                   sources=['gentcf.f'],)
ORBLOC = Extension(name='orbloc',
                   sources=['pmloca.f'],)
CLEMTP = Extension(name='clemtp',
                   sources=['clemtp.f'])
CLPLTP = Extension(name='clpltp',
                   sources=['clpltp.f'],
                   library_dirs=['/usr/lib/lapack','/usr/lib/pyshared/python2.7/scipy/lib/lapack'],
                   libraries=['flapack.so',])
#SOLPOL = Extension(name='solpol',
#                   sources=['solpol.f'])

# --- Install libbbg!

setup(name='LIBBBG',
      version='11.03a',
      description='Libraries for BBG packages',
      author='Bartosz Błasiak',
      author_email='globula@o2.pl',
      url='http://www.ex.no/pymod/m1',
      #packages=['solpol'],
      py_modules=['dma','gaussfreq','units',
                  'utilities','utilities2',
                  'dipderiv','re_templates',
                  'mpfit','solpol',],
      ext_modules=[GENTCF,ORBLOC,CLEMTP],
     )
