|ProjectStatus|_ 

.. |ProjectStatus| image:: http://www.repostatus.org/badges/latest/WIP.svg
.. _ProjectStatus: https://www.repostatus.org/#wip

### Current to do list (for devs) | Official documentation starts below

TO DO:
- in bt_analyzesources, nudge toward centering the specified frequencies exactly on warping frequency
- how to deal with outliers in lvl 2
- random gaussian field theory to perform multiple testing correction over power spectra and TGMs?

LOW URGENCY:
- Fix dots in the asymmetry plot in bt_selectsource
- Add licenses! 
- Update tutorial 4 to 6
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
| Brain Time Toolbox term | description |
| --- | --- |
| Clock time (CT) | data ordered by sequences of seconds |
| Brain time (BT) | data transformed according to the dynamics of a warping signal |
| Warping signal | an oscillatory signal that is assumed to clock the cognitive process of interest |
| Warping source | a data structure that contains potential warping signals (e.g., local field potentials, independent component analysis components, virtual channels) |
| Brain time warping | a transformation of clock time data that employs the phase dynamics of the warping signal |
| Periodicity | rhythmic patterns of classification performance |

| General term | description |
| --- | --- |
| Temporal generalization matrix (TGM) | training a classifier on one timepoint, and testing it on all other timepoints in the data |
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
