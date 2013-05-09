#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Standard library for BBG packages 
"""

from distutils.core import setup

setup(name='LIBBBG',
      version='1.0',
      description='Libraries for BBG packages',
      author='Bartosz Błasiak',
      author_email='globula@o2.pl',
      url='http://www.ex.no/pymod/m1',
      #packages=['libbbg'],
      py_modules=['dma','gaussfreq','units',
                  'utilities','utilities2',
                  'dipderiv','re_templates']
     )
