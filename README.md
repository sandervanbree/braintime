[![HitCount](http://hits.dwyl.com/sandervanbree/braintime.svg)](http://hits.dwyl.com/sandervanbree/braintime)
[![Github All Releases](https://img.shields.io/github/downloads/sandervanbree/braintime/total.svg)]()

### Current to do list (for devs) | Official documentation starts below

TO DO:
- sense check clabel (trialinfo)
- in bt_analyzesources, nudge toward centering the specified frequencies exactly on warping frequency
- random gaussian field theory to perform multiple testing correction over power spectra and TGMs?

LOW URGENCY:
- In the paper, remind users to be very careful with high pass filtering; https://doi.org/10.1101/530220

### Start documentation:


# Brain Time Toolbox
<img src="https://i.imgur.com/OAEVqgM.png" width="200">

Warp electrophysiological data from clock time to brain time and analyze the periodicity of cognitive processes. By *Sander van Bree, María Melcón, Maria Wimber & Simon Hanslmayr*.

### Dependencies
- [MATLAB Signal Processing Toolbox](https://uk.mathworks.com/help/signal/getting-started-with-signal-processing-toolbox.html)
- [FieldTrip](http://www.fieldtriptoolbox.org/download/)
- [MVPA Light](https://github.com/treder/MVPA-Light)
- [Several functions](external) included in the toolbox

### Glossary
| **Toolbox term** | **description** |
| --- | --- |
| Clock time (CT) | data ordered by sequences of seconds |
| Brain time (BT) | data transformed according to the dynamics of a warping signal |
| Warping signal | an oscillatory signal that is assumed to clock the cognitive process of interest |
| Warping source | a data structure that contains potential warping signals (e.g., local field potentials, independent component analysis components, virtual channels) |
| Brain time warping | a transformation of clock time data that employs the phase dynamics of the warping signal |
| Periodicity | rhythmic patterns of classification performance |
|             |                                                 |
| **General term** |           |
| Classification timecourse (Diagonal) | training and testing a classifier on one timepoint in the data |
| Temporal generalization matrix (TGM) | training a classifier on one timepoint, and testing it on all other timepoints in the data |
| fast Fourier Transform (FFT)         | an algorithm to go between the frequency and time domain of a signal | 
| Generalized eigendecomposition (GED) | a mathematical approach to combine signals in different warping sources into a single timeseries |
| Independent component analysis (ICA) | a mathematical approach to separating electrophysiological data into additive subcomponents |


## Table of contents
[Installation](#installation)

[Introduction](#introduction)

## Installation
**Option 1:** Download files

1. Download ZIP on the top right of this page (green "code" button)
2. Unzip to preferred folder
3. Run the function [setup_braintime](setup).

**Option 2:** Git client or checkout with SVN

Type:
```java
git clone https://github.com/sandervanbree/braintime.git
```
then run the function [setup_braintime](setup).

## Introduction
Summary of the toolbox here
