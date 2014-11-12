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
GENTCF = Extension(name='libbbg.gentcf',
                   sources=['libbbg/gentcf.f'],)
ORBLOC = Extension(name='libbbg.orbloc',
                   sources=['libbbg/pmloca.f'],)
CLEMTP = Extension(name='libbbg.clemtp',
                   sources=['libbbg/clemtp.f'])
#CLPLTP = Extension(name='libbbg.clpltp',
#                   sources=['libbbg/clpltp.f'],
#                   library_dirs=['/usr/lib/lapack','/usr/lib/pyshared/python2.7/scipy/lib/lapack'],
#                   libraries=['flapack.so',])
FT_LIB = Extension(name='libbbg.fourier.ft',
                   sources=['libbbg/fourier/ft.f'],)
#SOLPOL = Extension(name='solpol',
#                   sources=['solpol.f'])

# --- Install libbbg!

setup(name='LIBBBG',
      version='11.03a',
      description='Libraries for BBG packages',
      author='Bartosz Błasiak',
      author_email='blasiak.bartosz@gmail.com',
      url='no-page-yet',
      packages=['libbbg'],
#      packages=['letters','fourier',],
#      py_modules=['dma','gaussfreq','units',
#                  'utilities','utilities2',
#                  'dipderiv','re_templates','fourier',
#                  'mpfit','solpol','letters.greek','fourier.ft'],
      ext_modules=[GENTCF,ORBLOC,CLEMTP,FT_LIB],
     )
