[![Build Status](https://travis-ci.org/walaj/ginseng.svg?branch=master)](https://travis-ci.org/walaj/ginseng)

**License:** [GPLv3][license]

Tools for 1D and 2D genomic significance analysis

Table of contents
=================

  * [Installation](#installation)
  * [Description](#description)
  * [Components](#components)
    * [Swap](#swap)
    * [Sim](#sim)
  * [Example Recipes](#examples-recipes)
  * [Attributions](#attributions)

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

## note that the Apophenia library must be installed after the GSL install
##
##    git clone https://github.com/walaj/apophenia-clone
##    cd apophenia-clone
##    ./configure --prefix=`pwd`
##    make CFLAGS=-I$GSLHOME LIBS="-lgsl -lgslcblas"
##    make install
##    export APOPHOME=`pwd`

git clone --recursive https://github.com/walaj/ginseng
cd ginseng
./configure --with-apophenia=$APOPHOME --with-gsl=$GSLHOME ## must give full (not relative) paths here
make
make install
```

Description
-----------
*ginseng* is set of tools for determing genome-wide signficance for 2-dimensional 
genomic data. The primary use case is for structural variations (aka rearrangements) in cancer, but
the methods in principle could be applied to other data sets (eg HiC interactions).

Components
----------

#### Swap
Swap performs a non-parametric test for enrichment of 2D connections between or within tracks.
2D events are expected in the form of a list of BEDPE files. Tracks are BED tracks of genomic
features. BED tracks are input with a ``bed_list`` file, which has the form ``NAME,PATH_TO_BED`` for each
line.

```
## run the swaps
ginseng swap -n $num_matrices -b $num_bins -k $num_swaps \
	-p $cores -i $bedpe_list -a $my_id -B $track_list -A $animation_step \
	--min-span $min --max-span $max

## plot the results
R/swap-animate.R -a $my_id
```

#### Sim
```
## simulate 1M rearrangements with a 1/L power law distribution using additive generation
COV=data/cov.bed ## some BED file with covered regions (eg sufficient mapping quality)
ginseng sim $COV -m A -N 1000000 > data/model_additive.bed

## simulate 1M rearrangements with a 1/L^2 power law, using multiplicative generation
ginseng sim $COV -m M -N 1000000 > data/model_multiplicative.bed

```


Support
-------
This project is being actively developed and maintained by Jeremiah Wala (jwala@broadinstitute.org). 

Attributions
------------
Ginseng is being developed in the Beroukhim Lab of Dana-Farber and the Broad Institute, and the
Imielinski Lab of Cornell University and the New York Genome Center.

Development, support, guidance, testing:
* Ofer Shapira - Computational Biologist, Dana Farber Cancer Institute
* David Craft - Assistant Professor, Radiation Oncology, Massachusetts General Hospital
* Steve Schumacher - Computational Biologist, Dana Farber Cancer Institute
* Marcin Imielinski - Assistant Professor, Cornell University
* Rameen Beroukhim - Assistant Professor, Harvard Medical School

[license]: https://github.com/walaj/ginseng/blob/master/LICENSE
