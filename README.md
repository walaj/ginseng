[![Build Status](https://travis-ci.org/jwalabroad/SwapSV.svg?branch=master)](https://travis-ci.org/jwalabroad/SwapSV)

# SwapSV

Installation
------------

```
###
git clone --recursive https://github.com/jwalabroad/SwapSV

### if on Broad Institute servers, add GCC-4.9
reuse -q GCC-4.9

############## DOWNLOAD BOOST (if not installed) ###############
wget https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.gz
tar -xvzf boost_1_61_0.tar.gz
## we only user header-only libraries, so no compiling of Boost is needed

############### COMPILE AND INSTALL ###############
./configure --with-boost=<path_to_boost>  ## e.g. ~/boost_1_61_0
make
```
