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

**1.1 Loading clock time data**

`braintime` works with FieldTrip formatted electrophysiological data. This can be electroencephalography (EEG), magnetoencephalography (MEG), or intracranial data. The starting data is called clock time data—this will be warped.

**1.2 Loading warping sources data**

Warping sources are the data structure containing the to-be-selected warping signal, which is used to warp clock time data. Warping sources may be obtained separately from clock time data, or extracted from it (using independent component analysis, or the selection of a few channels). Please ensure that your warping source data has 0.5s of additional time extra, before the start and after of your time window of interest, to facilitate step 1.3.

When clock time data and warping sources are dependent, please read [Is it circular to warp to warping sources obtained from my clock time data?](#is it circular to warp to warping sources obtained from my clock time data?).

**1.3 Time frequency analysis of warping sources**

Each warping source contains potential warping signals. [bt_analyzesources](warpingsource/bt_analyzesources.m) performs a time frequency analysis on all warping sources, detecting potential warping signals based on your preferences. These preferences include the frequency range of interest assumed to clock your cognitive process (e.g., 8 to 12 Hz for attention), and the time window of interest that you wish to analyze (e.g. 0 to 1 second, this should match the window you wish to warp). In addition, you can choose one of two methods to cut the data, 'consistenttime', or 'cutartefact'. For details on their relative merit, see [Which cutmethod to choose?](#which cutmethod to choose?)


## Operation 2: Periodicity analysis

Be sure to label your two classes of trials in a field called "clabel" (a vector of 1's and 2's lining up with your trial field).

## Toolbox considerations

### Which cutmethod to choose?

### Is it circular to warp to warping sources obtained from my clock time data?

