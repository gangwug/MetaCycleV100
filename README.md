# MetaCycle
This repository contains files for MetaCycle version 1.0.0 (the first version).

## Introduction
MetaCycle is a newly developed R package for evaluating periodicity in large scale time-series datasets. This package provides 
two functions-**meta2d** and **meta3d**. For analyzing time-series datasets without individual information, **meta2d** is suggested, 
which incorporates [ARSER](https://github.com/cauyrd/ARSER), [JTK_CYCLE](http://openwetware.org/wiki/HughesLab:JTK_Cycle) and
[Lomb-Scargle](http://research.stowers-institute.org/efg/2005/LombScargle/) in the detection of interested rhythms. For analyzing 
time-series datasets with individual information, **meta3d** is suggested, which takes use of Lomb-Scargle to analyze time-series data individual by individual and gives out integrated values based on analysis result of each individual.

## Installation
Use **devtools** to install this version from Github:

```r
# install 'devtools' in R(>3.0.2)
install.packages("devtools")
# install MetaCycle version 1.0.0
devtools::install_github('gangwug/MetaCycleV100')
```

## Usage
```r
library(MetaCycle)
# see detail introduction and associated application examples about meta2d or meta3d
?meta2d
?meta3d
```

## License
This package is free and open source software, licensed under GPL(>= 2).
