libbbg
======

Common Python libraries (LIBBBG) are useful, easy to use Python modules containing functions and classes 
for daily computational chemistry calculation purposes. The following modules are now available:

 - `utilities` which stores various tools (.xyz file handling, rotation matrices etc)
 - `dma` contains DMA class (Distributed Multipole Analysis)
 - `units` contains various unit converters and handling with atom labels
 - `gaussfreq` contains a class for handling Gaussian09 anharmonic analysis file

After installation (see tutorial) you can easily import these modules by typing:

```
$ python
>>> import libbbg.utilities
```
or 
```
>>> from libbbg.utilities import *
```

The tutorial is under preparation.

### Installation ###

Installation prerequisites: 
- Python 2.6 or newer
- NumPy 1.6.1 or newer
- python scitools 0.8 or newer

To install the LIBBBG package you have to make the main directory for LIBBBG,
download the files on this directory and type:

```
$ ./install -p installatio_directory
```

in the main directory. Perhaps you might need super user priviledges if the installation 
path cannot be accessed in other way. Good Luck!

Tutorial
--------

### `units`

Type in your Python console:
```
from libbbg.units import UNITS
```
Now you have an access to the tools of units module.

## Units and standard converters

The class `UNITS` contains all of the necessary converters. Your Python classes can inherit the attributes of `UNITS` class.
To see all of the converters type:
```
t = UNITS()
print t
```
For instance, the converter from Bohrs to Angstroms is `UNITS().BohrToAngstrom`. Therefore, the following example shows how to convert distances by using Python:
```
dist_angs = dist_bohr * t.BohrToAngstrom
```
