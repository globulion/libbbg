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
GENTCF = Extension(name='libbbg.qm.gentcf',
                   sources=['libbbg/qm/gentcf.f'],)
CLEMTP = Extension(name='libbbg.qm.clemtp',
                   sources=['libbbg/qm/clemtp.f'])
ORBLOC = Extension(name='libbbg.qm.pmloca',
                   sources=['libbbg/qm/pmloca.f'],)
FT_LIB = Extension(name='libbbg.fourier.ft',
                   sources=['libbbg/fourier/ft.f'],)
EFPROT = Extension(name='libbbg.qm.efprot',
                   sources=['libbbg/qm/efprot.f'],)

#CLPLTP = Extension(name='libbbg.clpltp',
#                   sources=['libbbg/clpltp.f'],
#                   library_dirs=['/usr/lib/lapack','/usr/lib/pyshared/python2.7/scipy/lib/lapack'],
#                   libraries=['flapack.so',])
#SOLPOL = Extension(name='solpol',
#                   sources=['solpol.f'])

# --- Install libbbg!

setup(name='LIBBBG',
      version='1.0.2',
      description='Libraries for BBG packages',
      author='Bartosz Błasiak',
      author_email='blasiak.bartosz@gmail.com',
      url='no-page-yet',
      packages=['libbbg',
                'libbbg.qm',
                'libbbg.letters',
                'libbbg.fourier'],
      package_dir= {'libbbg'         : 'libbbg', 
                    'libbbg.qm'      : 'libbbg/qm', 
                    'libbbg.letters' : 'libbbg/letters',
                    'libbbg.fourier' : 'libbbg/fourier'},
      py_modules=['libbbg.dma',
                  'libbbg.gaussfreq',
                  'libbbg.units',
                  'libbbg.utilities','libbbg.utilities2',
                  'libbbg.dipderiv','libbbg.re_templates',
                  'libbbg.mpfit','libbbg.letters.greek',],
      ext_modules=[GENTCF,ORBLOC,CLEMTP,FT_LIB,EFPROT],
     )

