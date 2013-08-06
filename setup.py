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
                   sources=['gentcf.f'],
                   f2py_options=["--fcompiler='gnu95'"])
                   

# --- Install libbbg!

setup(name='LIBBBG',
      version='11.03a',
      description='Libraries for BBG packages',
      author='Bartosz Błasiak',
      author_email='globula@o2.pl',
      url='http://www.ex.no/pymod/m1',
      #packages=['libbbg'],
      py_modules=['dma','gaussfreq','units',
                  'utilities','utilities2',
                  'dipderiv','re_templates'],
      ext_modules=[GENTCF],
     )
