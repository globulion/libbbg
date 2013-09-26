### REGULAR EXPRESSIONS TEMPLATES
import re

__all__ = ['re_real_e','re_dbl_fort']

# scientific notation with E/e
re_real_e=re.compile(r'-?\d\.\d+[Ee][+\-]\d\d?')

# FORTRAN double precision number format
re_dbl_fort = re.compile(r'(\d*\.\d+)[dD]([-+]?\d+)')