[![Build Status](https://travis-ci.org/walaj/ginseng.svg?branch=master)](https://travis-ci.org/walaj/ginseng)

**License:** [GPLv3][license]

Tools for 1D and 2D genomic significance analysis

Installation
------------
```
## note that GNU scientific library (GSL) and libsql is required. See: http://apophenia.info/setup.html
## If GNU scientific library is not in standard location, set GSL_CONFIG environment
## variable to path to gsl-config. Example setup for GSL below (if not already on system):
##    wget http://www.localmsp.org/gnu/gsl/gsl-2.2.1.tar.gz ## download GSL
##    tar xfzv gsl-2.2.1 && cd gsl-2.2.1
##    ./configure --prefix=`pwd` && make && make install
##    export GSL_HOME=`pwd`
##    export GSL_CONFIG=$GSL_HOME/bin/gsl-config
##    export LD_LIBRARY_PATH=$GSL_HOME/lib:$LD_LIBRARY_PATH

git clone --recursive https://github.com/jwalabroad/ginseng
cd ginseng
make apophenia ## install the apophenia libs
./configure
make
make install
```

Support
-------
This project is being actively developed and maintained by Jeremiah Wala (jwala@broadinstitute.org). 

Attributions
------------
Ginseng is being developed in the Beroukhim Lab of Dana-Farber and the Broad Institute, and the
Imielinski Lab of Cornell University and the New York Genome Center.

Development, support, guidance, testing:
* Steve Schumacher - Computational Biologist, Dana Farber Cancer Institute
* Marcin Imielinski - Assistant Professor, Cornell University
* Rameen Beroukhim - Assistant Professor, Harvard Medical School

[license]: https://github.com/jwalabroad/ginseng/blob/master/LICENSE
