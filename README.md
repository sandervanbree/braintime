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

[Operation 1: Brain time warping](#operation-1-brain-time-warping)

[Operation 2: Periodicity analysis](#operation-2-periodicity-analysis)

[Toolbox considerations](#toolbox-considerations)

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
What's the premise behind the Brain Time Toolbox (`braintime`)? Insofar a cognitive process is clocked by oscillations, analyses of the process's dynamics require a factoring in of the oscillations' dynamics. To this end, the toolbox implements **brain time warping**, an algorithm to transform electrophysiological data based on relevant oscillations selected by the user—the **warping signal**. This changes the data from **clock time** (seconds), closer to **brain time** (cycles). This completes the first operation of `braintime`: *brain time warping*. 

The second operation is *periodicity analysis*. Here, `braintime` uses multivariate pattern analysis (MVPA) to detect the effects of the applied data transformation. Has evidence of the rhythmicity of the cognitive patterns increased?

Below, we explain how to perform both operations step by step.

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

There are two parameters you can change in [bt_clocktobrain](bt_clocktobrain/bt_clocktobrain.m). You can choose to set the clock time signal as a basic stationary sinusoid (``` cfg.warpmethod = 'sinusoid'``` ) or a smoothed version of the warping signal's waveshape (``` cfg.warpmethod = 'waveshape'``` ).

> :bulb: For asymmetric data, such as theta oscillations in intracranial rodent data, using waveshape as a warpmethod is the natural choice.

After [bt_clocktobrain](bt_clocktobrain/bt_clocktobrain.m), your data has completed its transformation from clock to brain time. The time axis is now formatted as sequences of cycles, instead of seconds. You may opt to continue analyses outside of the toolbox, or test for periodic patterns in the data using `braintime`'s second operation described next. In case of the former, please read the previously linked circularity information.


## Operation 2: Periodicity analysis

You may wish to use `braintime` to test whether brain time warping has accentuated dynamic patterns in the data. This is appropriate if your data comprises two underlying conditions which are predicted to yield fluctuating classification evidence. For example, we report a spatial attention dataset where we predict that evidence for motion direction rapidly alternates from low to high as a function of alpha phase (the warping signal). Thus, the second operation requires your brain time warped data to consist of such classes.

Before getting started, ensure the original clock time data as well as the brain time warped data have their classes labeled. Specifically, create a field called "clabel", with a vector of 1's and 2's corresponding to the condition of each trial. If your trial sequence is 'left left right', the clabel field should contain '1 1 2'.

The pipeline is that each participant's data needs to undergo step 2.1 to 2.3, separately for both clock and brain time data. Then, each participant's first level result for clock time must be added as a field to a large group structure, and the same for brain time results. These two group structures (e.g. ```ct_stats1``` and ```bt_stats1```) can separately be used as input for step 2.4 to obtain clock time and brain time group level results, allowing you to compare periodicity in clock time and brain time data.

**2.1 Multivariate pattern analysis**

The first step is to use `MVPA-Light` to classify the data. You may opt to classify across time [mv_classify_across_time.m](https://github.com/treder/MVPA-Light/blob/master/mv_classify_across_time.m). This method tests whether a classifier can separate both classes of data across time in the trials. Alternatively, you can apply temporal generalization by classifying using [mv_classify_timextime.m](https://github.com/treder/MVPA-Light/blob/master/mv_classify_timextime.m). This method tests for temporal generalization of classification. That is, it tests to what extent a classifier trained on one timepoint generalizes its performance to other timepoints.

> :bulb: You can change a variety of parameters when calling `MVPA-Light`, described briefly in [tutorial 2](tutorial/tutorial2_periodicity.m). For more information, check `MVPA-Light`'s awesome [tutorials](https://github.com/treder/MVPA-Light/tree/master/examples).

**2.2 Quantify periodicity in classifier performance**

If the warping signal orchestrates the dynamics of your cognitive function, operationalized by your two classes of data, then the classifier's performance may tap into these dynamics. Thus, [bt_quantify](periodicity/bt_quantify.m) tests and quantifies periodic patterns in the classifier performance, whether that is in the classifier output obtained from [mv_classify_across_time.m](https://github.com/treder/MVPA-Light/blob/master/mv_classify_across_time.m) (1 dimensional; "diagonal"), or [mv_classify_timextime.m](https://github.com/treder/MVPA-Light/blob/master/mv_classify_timextime.m) (2 dimensional; "temporal generalization matrix" (TGM)). In both cases, you get a "periodicity spectrum" which is just a power spectrum of your 1D or 2D classifier performance.

For TGMs, you may opt to perform the periodicity analysis over either the 2 dimensional matrix itself, or its autocorrelation map. Not sure what is better? Check out "[should I perform bt_quantify over the TGM itself, or its autocorrelation map?](#should-i-perform-bt_quantify-over-the-tgm-itself-or-its-autocorrelation-map)".

You also need to specify a range of periodicity frequencies. At which rate do you expect periodic patterns in classifier performance to arise? Finally, you can choose a reference dimension. If you choose ```cfg.refdimension = clocktime```, periodicity in classifier performance will be displayed as a function of cycles per seconds (Hz). Alternatively, if you are quantifying brain time warped data, ```cfg.refdimension = clocktime``` is more appropriate. This references the peroidicity as a function of cycles per second, *normalized* to that participants' warping frequency.

> :bulb: Let's say you warp a participant's data to 11 Hz. Then with ```cfg.refdimension = clocktime```, a peak at the warping frequency will be at 11 Hz, but with ```cfg.refdimension = clocktime``` it will show up at 1 Hz, as the frequencies are normalized to the warped frequency (11/11 = 1).

**2.3 First level statistics**

How does the participant's quantified periodicity compare against the null distribution? [bt_statslevel1](periodicity/bt_statslevel1.m) takes the output from [bt_quantify](periodicity/bt_quantify.m) and performs classification ```cfg.numperms1 = n1``` times over, each time randomly shuffling the classification labels and obtaining x "permuted" periodicity spectra. This provides a null distribution that quantifies how much periodicity is in the data when the class structure is destroyed, setting things up for p-value estimation on the group level.

Now, repeat [bt_statslevel1](periodicity/bt_statslevel1.m), and send its output to a separate field in a group structure (e.g. ```[ct_stats1{subj}] = bt_statslevel1(cfg,data,quant)```). The next step requires this format to perform 2nd level statistics.

**2.4 Second level statistics**

We now have an idea of the periodicity in single participants. But what about the group level? [bt_statslevel2](periodicity/bt_statslevel2.m) employs the following algorithm to generate second level statistics:

1) ```cfg.numperms2 = n2``` times, do the following: randomly grab one of each participants' n1 permuted spectra and average those spectra into one spectrum. This yields n2 spectra.
2) For each frequency, compare the empirical periodicity with the distribution of permuted periodicity to derive a p-value.
3) Plot the average empirical periodicity spectrum and display the p-value for each frequency.

> :bulb: In brain time results, p-values at 0.5 Hz, 1 Hz, and 2 Hz are exempted from multiple testing correction, as periodicity is predicted at one or multiple of these rates depending on the underlying structure'.

> :bulb: `braintime` z-score participants' spectra to enable second level statistics. For details on how it does so, check out "[How does the toolbox z-score periodicity spectra?](#how-does-the-toolbox-z-score-periodicity-spectra)".

## Toolbox considerations

### Which cutmethod to choose?

### Is it circular to warp to warping sources obtained from my clock time data?

### Should I perform bt_quantify over the TGM itself, or its autocorrelation map?

### How does the toolbox z-score periodicity spectra?

