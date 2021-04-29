[![HitCount](http://hits.dwyl.com/sandervanbree/braintime.svg)](http://hits.dwyl.com/sandervanbree/braintime)
[![Github All Releases](https://img.shields.io/github/downloads/sandervanbree/braintime/total.svg)]()

# Brain Time Toolbox
<img src="https://i.imgur.com/OAEVqgM.png" width="200">

Warp electrophysiological data from clock time to brain time and analyze the periodicity of cognitive processes. By *Sander van Bree, María Melcón, Luca Kolibius, Casper Kérren, Maria Wimber, Simon Hanslmayr*.

### Dependencies
- [MATLAB Signal Processing Toolbox](https://uk.mathworks.com/help/signal/getting-started-with-signal-processing-toolbox.html)
- [FieldTrip](https://github.com/fieldtrip/fieldtrip/)
- [MVPA Light](https://github.com/treder/MVPA-Light)
- [Several functions](external) included in the toolbox

## Table of contents
[Installation](#installation)

[Glossary](#glossary)

[Introduction](#introduction)

[Brain time warping](#operation-1-brain-time-warping)

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

> :warning: If your clock time data and warping sources are not independent, please read "[is it circular to warp to warping sources obtained from my clock time data?](#is-it-circular-to-warp-to-warping-sources-obtained-from-my-clock-time-data?)".

**1.3 Time frequency analysis of warping sources**

Each warping source contains potential warping signals. [bt_analyzesources](warpingsource/bt_analyzesources.m) performs a time frequency analysis on all warping sources, detecting potential warping signals based on your preferences. These preferences include methods of time frequency analysis, the frequency range of interest assumed to clock your cognitive process (e.g., 8 to 12 Hz for attention), and the time window of interest that you wish to analyze (e.g. 0 to 1 second, this should match the window you wish to warp). In addition, you can choose one of two methods to cut the data, 'consistenttime', or 'cutartefact', both with their own relative merits (see "[which cutmethod to choose?](#which-cutmethod-to-choose)").

[bt_analyzesources](warpingsource/bt_analyzesources.m) collects a whole bunch of information about the warping sources. The phase of the underlying warping signals (using two methods, described in 1.5), its average waveshape and its asymmetry, and a ranking of all warping sources based on your preferences.

> :bulb: You can adjust the ranking of warping sources by their topography by creating a topography using [bt_templatetopo](topography/bt_templatetopo.m) and changing cfg.rankmethod to 'temptopo'.

**1.4 Selecting a warping source and signal**

Now that the toolbox has collected all the relevant time frequency details of warping sources, it's time to make a selection using [bt_selectsource](warpingsource/bt_selectsource.m). If desired, you can input a layout structure to plot each warping source's topography. For each warping source, [bt_selectsource](warpingsource/bt_selectsource.m) plots a time frequency spectrum, a power spectrum, the topography of the warping source, and waveshape details. For each warping source, a warping signal is highlighted, set at the highest power oscillation in the frequency range of interest. Factoring in its characteristics, you can select a warping signal. `braintime` uses the phase of dynamics of this warping signal to brain time warp the clock time data.

**1.5 Warping clock to brain time**

This is where the magic happens. [bt_clocktobrain](bt_clocktobrain/bt_clocktobrain.m) takes the phase of the chosen warping signal and [dynamically time warps](https://en.wikipedia.org/wiki/Dynamic_time_warping) (DTW) it to the phase of a stationary sinusoid of the same frequency. Thereby, DTW enables a readout of where clock time (stationary sinusoid) falls out of tune with brain time (phase of warping signal) by attempting to minimize their difference. That minimization process yields a warping path, which tells `braintime` the samples that need to be repeated for clock and brain time to better align. [bt_clocktobrain](bt_clocktobrain/bt_clocktobrain.m) repeats those samples in the original data, cycle-by-cycle and trial-by-trial, before squeezing things back down to the original data length.

There are two parameters you can change in [bt_clocktobrain](bt_clocktobrain/bt_clocktobrain.m). You can choose to set the clock time signal as a basic stationary sinusoid (warpmethod = 'sinusoid') or a smoothed version of the warping signal's waveshape (warpmethod = 'waveshape').

> :bulb: For asymmetric data, like theta oscillations in intracranial rodent data, using waveshape as a warpmethod is the natural choice

After [bt_clocktobrain](bt_clocktobrain/bt_clocktobrain.m), your data has completed its transformation from clock to brain time. The time axis is now formatted as sequences of cycles, instead of seconds.


## Operation 2: Periodicity analysis

Be sure to label your two classes of trials in a field called "clabel" (a vector of 1's and 2's lining up with your trial field).

## Toolbox considerations

### Which cutmethod to choose?

### Is it circular to warp to warping sources obtained from my clock time data?

