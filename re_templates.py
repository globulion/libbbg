### REGULAR EXPRESSIONS TEMPLATES

__all__ = ['re_real_e']

# scientific notation with E/e
re_real_e=r'-?\d\.\d+[Ee][+\-]\d\d?'

# FORTRAN double precision number format
re_dbl_fort = re.compile(r'(\d*\.\d+)[dD]([-+]?\d+)')