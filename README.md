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
>>> import utilities
```
or 
```
>>> from utilities import *
```

The tutorial is under preparation.

### Installation ###

To install the LIBBBG package you have to download the files to your computer and type:

```
$ sudo python setup.py install
```

Good Luck!
