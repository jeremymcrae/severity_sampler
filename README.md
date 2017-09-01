[![Build Status](https://travis-ci.org/jeremymcrae/severity_sampler.svg?branch=master)](https://travis-ci.org/jeremymcrae/severity_sampler)

### severity_sampler: functional severity of *de novo* mutations
Calculates the significance of seeing N *de novo* mutations with a summed
severity in a given gene.

#### Install
```sh
# you might need cython installed before this e.g. pip install cython
pip install git+git://github.com/jeremymcrae/severity_sampler.git --user

# Alternatively:
git clone https://github.com/jeremymcrae/severity_sampler.git
cd mupit
python setup.py install --user

```

#### Usage:
```sh
python bin/simulate_severity.py \
    --de-novos PATH_TO_DENOVOS \
    --cadd PATH_TO_SNVS_CADD \
    --constraint PATH_TO_REGIONAL_CONSTRAINT \
    --cache CACHE_FOLDER \
    --genome-build grch37 \
    --output OUTPUT.txt
```

CADD scores for SNVs are available from the
[CADD website](http://cadd.gs.washington.edu/download) (see file for all
possible SNVs of GRCh37/hg19). Regional constraint scores are available for
genes from
[Samocha et al.](http://www.biorxiv.org/content/early/2017/06/12/148353), Table S2.
