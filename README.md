triplo_gt
=========

Genotypcaller for genotypes of ploidy up to three or more.


---

#### Binaries:

Binaries in this repository were built in Ubuntu Linux.  Use at your own risk.


#### Dependencies:

g++-4.9
The function 'pwin' uses <regex>.  This was not fully supported in earlier versions of g++.  For example, g++ 4.6.4 did not throw a compile time or run time error.  But the resulting program did not return expected results, it did not count many genotypes.  The version g++-4.9 appears to resolve this.  How to get 
[g++-4.9 for Ubuntu](http://askubuntu.com/questions/428198/getting-installing-gcc-g-4-9-on-ubuntu).

[GNU multiple precision arithmetic library](https://gmplib.org/).
Note that [C++ support is not enabled by default by GMP](http://stackoverflow.com/a/22803223).  Please build the GMP with `./configure --prefix=/usr/local --enable-cxx`.

[Boost C++ libraries](http://www.boost.org/)

[SAMtools](http://samtools.sourceforge.net/)


