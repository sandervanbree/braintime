[![HitCount](http://hits.dwyl.com/sandervanbree/braintime.svg)](http://hits.dwyl.com/sandervanbree/braintime)
[![Github All Releases](https://img.shields.io/github/downloads/sandervanbree/braintime/total.svg)]()

# Brain Time Toolbox
<img src="https://i.imgur.com/OAEVqgM.png" width="200">

Warp electrophysiological data from clock time to brain time and analyze the periodicity of cognitive processes. By *Sander van Bree, María Melcón, Luca Kolibius, Casper Kérren, Maria Wimber, Simon Hanslmayr*.

### Dependencies
- [MATLAB Signal Processing Toolbox](https://uk.mathworks.com/help/signal/getting-started-with-signal-processing-toolbox.html)
- [FieldTrip](http://www.fieldtriptoolbox.org/download/)
- [MVPA Light](https://github.com/treder/MVPA-Light)
- [Several functions](external) included in the toolbox

## Table of contents
[Installation](#installation)

[Glossary](#glossary)

[Introduction](#introduction)

[Brain time warping](#Operation 1: Brain time warping)

## Installation
**Option 1:** Download files

1. [Download ZIP](https://github.com/sandervanbree/braintime/archive/refs/heads/master.zip)
2. Unzip to preferred folder
3. Run the function [setup_braintime](setup).

**Option 2:** Git client or checkout with SVN

Type:
```java
git clone https://github.com/sandervanbree/braintime.git
```
then run the function [setup_braintime](setup).

## Glossary
| **Toolbox term** | **description** |
| --- | --- |
| Brain time warping | a transformation of clock time data based on the phase dynamics of the warping signal |
| Clock time (CT) | time as sequences of seconds
| Brain time (BT) | time as sequences of cycles of an oscillation of interest
| Warping signal | an oscillatory signal that is assumed to clock the cognitive process of interest |
| Warping source | a data structure that contains potential warping signals (e.g., local field potentials, independent component analysis components, virtual channels) |
| Periodicity | rhythmic patterns of classification performance |


## Introduction
What is the premise behind the Brain Time Toolbox (`braintime`)? Insofar a cognitive process is clocked by oscillations, analyses of the process's dynamics require a factoring in of the oscillations' dynamics. To this end, the toolbox implements **brain time warping**, an algorithm to transform electrophysiological data based on relevant oscillations selected by the user—the **warping signal**. This changes the data from **clock time** (seconds), closer to **brain time** (cycles). This completes the first operation of `braintime`: *brain time warping*. 

The second operation is *periodicity analysis*. Here, `braintime` uses multivariate pattern analysis (MVPA) to detect the effects of the applied data transformation. Has evidence of the rhythmicity of the cognitive patterns increased?

Below, we demonstrate how `braintime` achieves both operations step-by-step.

## Operation 1: Brain time warping

**Loading data**

`braintime` works with FieldTrip formatted electrophysiological data. 


