### REGULAR EXPRESSIONS TEMPLATES
"""
Reqular expressions templates module
------------------------------------

It contains patterns for re module.
Each patterns can be used two fold:

  1) as a string r'[something]'
  2) as a compiled object (re.compile)

The name for the latter contains end:
_c (compiled object), for example:

re_real_e = r'-?\d\.\d+[Ee][+\-]\d\d?'
                and 
re_real_e_c = re.compile(re_real_e)
"""
import re

__all__ = ['re_real_e'  ,'re_real_e_c',
           're_dbl_fort','re_dbl_fort_c',]

# scientific notation with E/e
re_real_e   = r'-?\d\.\d+[Ee][+\-]\d\d?'
re_real_e_c = re.compile(re_real_e)

# FORTRAN double precision number format
re_dbl_fort   = r'(\d*\.\d+)[dD]([-+]?\d+)'
re_dbl_fort_c = re.compile(re_dbl_fort)