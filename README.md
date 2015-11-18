# MetaCycle
This repository contains files for MetaCycle version 1.0.0 (the first version).

MetaCycle is a newly developed R package for evaluating periodicity in large scale time-series datasets. This package provides 
two functions-**meta2d** and **meta3d**. For analyzing time-series datasets without individual information, **meta2d** is suggested, 
which incorporates [ARSER](https://github.com/cauyrd/ARSER), [JTK_CYCLE](http://openwetware.org/wiki/HughesLab:JTK_Cycle) and
[Lomb-Scargle](http://research.stowers-institute.org/efg/2005/LombScargle/) in the detection of intereseted rhythms. For analyzing 
time-series datasets with individual information, **meta3d** is suggested, which takes use of Lomb-Scargle to analyze time-series data 
individual by individual and gives out integrated values based on analysis result of each individual.
