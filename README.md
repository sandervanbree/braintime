<!--- ![Hits](https://hitcounter.pythonanywhere.com/count/tag.svg?url=https%3A%2F%2Fgithub.com%2Fsandervanbree%2Fbraintime) --->

<!--- [![Github All Releases](https://img.shields.io/github/downloads/sandervanbree/braintime/total.svg)]() --->

# Brain Time Toolbox
<img src="https://i.imgur.com/2QP4cnb.png" width="200">

Warp electrophysiological data from clock time to brain time and analyze the periodicity of cognitive processes. A MATLAB toolbox by *Sander van Bree, María Melcón, Luca Kolibius, Casper Kerrén, Maria Wimber, Simon Hanslmayr*. 

Check out the [brain time paper](https://www.biorxiv.org/content/10.1101/2021.06.09.447763v1) for context.

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
| :-- | :-- |
| Brain time warping | a transformation of clock time data based on the phase dynamics of the warping signal |
| Clock time (CT) | time as sequences of seconds |
| Brain time (BT) | time as sequences of cycles of a coordinating brain oscillation |
| Warping signal | an oscillatory signal that is assumed to clock the cognitive process of interest |
| Warping source | data structure containing potential coordinating brain oscillations used for brain time warping (e.g., local field potentials, independent component analysis components, virtual channels) |
| Periodicity | fluctuating patterns of a neural signature |


## Introduction
What's the premise behind the Brain Time Toolbox (`braintime`)? Insofar a cognitive process is clocked by oscillations, analyses of the process's dynamics require a factoring in of the oscillations' dynamics. To this end, the toolbox implements **brain time warping**, an algorithm that transforms electrophysiological data based on relevant oscillations selected by the user — the **warping signal**. This changes the data from **clock time** (seconds), closer to **brain time** (cycles). This completes the first operation of `braintime`: *brain time warping*. 

The second operation is *periodicity analysis*. Here, `braintime` uses multivariate pattern analysis (MVPA) to detect the effects of the applied data transformation. Has evidence of the rhythmicity of the cognitive patterns increased?

Below, we explain how to perform both operations step by step. To see the steps in practice, check out `braintime`'s [tutorials](/tutorial). For extensive methodological documentation, check out the main paper's supplementary material. This figure gives a visual rundown of all steps:

<img src="https://i.imgur.com/uYfJnCm.png" width="1000">

&nbsp;

## Operation 1: Brain time warping

**1.1 Loading clock time data**

`braintime` works with FieldTrip formatted electrophysiological data. This can be electroencephalography (EEG), magnetoencephalography (MEG), or intracranial data. The starting data is called clock time data — this will be warped.

**1.2 Loading warping sources data**

Warping sources are the data structure containing the to-be-selected warping signal, which is used to warp clock time data. Warping sources may be obtained separately from clock time data, or extracted from it (using independent component analysis, or the selection of a few channels). Please ensure that your warping source data has 0.5s of additional time extra (before and after your time window of interest) to facilitate step 1.3.

> :warning: If your clock time data and warping sources are not independent, please read [Is it circular to warp to warping sources obtained from my clock time data?](#is-it-circular-to-warp-to-warping-sources-obtained-from-my-clock-time-data)

**1.3 Time frequency analysis of warping sources**

Each warping source contains potential warping signals. [bt_analyzesources](warpingsource/bt_analyzesources.m) performs a time frequency analysis on all warping sources, detecting potential warping signals based on your preferences. These preferences include methods of time frequency analysis (parameters fed into FieldTrip), the frequency range of interest assumed to clock your cognitive process (e.g., 8 to 12 Hz for attention), and the time window of interest that you wish to analyze (e.g. 0 to 1 second, this should match the window you wish to warp). In addition, you can choose one of two methods to cut the data, ```cfg.warpmethod = 'consistenttime'``` or ```'cutartefact'```, both with their own relative merits (see "[which cutmethod to choose?](#which-cutmethod-to-choose)").

[bt_analyzesources](warpingsource/bt_analyzesources.m) collects a whole bunch of information about the warping sources, such as the phase of the underlying warping signals (obtained using two methods, described in 1.5), its average waveshape and its asymmetry, and a ranking of all warping sources based on your preferences.

> :bulb: You can adjust the ranking of warping sources by their topography by creating a topography using [bt_templatetopo](topography/bt_templatetopo.m) and changing to ```cfg.rankmethod = 'temptopo'```.

**1.4 Selecting a warping source and signal**

Now that the toolbox has collected all the relevant time frequency details of warping sources, it's time to make a selection using [bt_selectsource](warpingsource/bt_selectsource.m). If desired, you can input a layout structure to plot each warping source's topography. For each warping source, [bt_selectsource](warpingsource/bt_selectsource.m) plots a time frequency spectrum, a power spectrum, the topography of the warping source, and waveshape details. For each warping source, a warping signal is highlighted, set at the highest power frequency in the range of interest. Factoring in its characteristics, you can select a warping signal. `braintime` uses the phase dynamics of this warping signal to brain time warp the clock time data in the next step.

**1.5 Warping clock to brain time**

This is where the magic happens. [bt_clocktobrain](bt_clocktobrain/bt_clocktobrain.m) takes the phase of the chosen warping signal and [dynamically time warps](https://en.wikipedia.org/wiki/Dynamic_time_warping) (DTW) it to the phase of a stationary signal of the same frequency. DTW enables a readout of where clock time (stationary signal) falls out of tune with brain time (phase of warping signal) by attempting to minimize their difference. That minimization process yields a warping path, which tells `braintime` the samples that need to be repeated for clock and brain time to better align. [bt_clocktobrain](bt_clocktobrain/bt_clocktobrain.m) repeats those samples in the original data, cycle-by-cycle and trial-by-trial, before squeezing things back down to the original data length.

There are two parameters you can change in [bt_clocktobrain](bt_clocktobrain/bt_clocktobrain.m). You can choose to set the clock time signal as a basic stationary sinusoid (``` cfg.warpmethod = 'sinusoid'``` ) or a smoothed version of the warping signal's waveshape (```cfg.warpmethod = 'waveshape'```).

> :bulb: For data with asymmetric waves such as theta oscillations in intracranial rodent data, using waveshape as a warpmethod is the natural choice.

There are also two "phase methods", which are separate phase estimations of the warping signal both sneakily obtained by [bt_analyzesources](warpingsource/bt_analyzesources.m). One is a regular method of phase estimation (using the Fast Fourier Transform of the warping frequency, selected with ```cfg.phasemethod = 'fft'```, and a special method called [General Eigenvalue Decomposition](http://mikexcohen.com/data/Cohen_STfilter.pdf) (GED), selected with ```cfg.phasemethod = 'ged'```. The former is classical and specific to the warping source, the latter is novel and more holistically estimates phase. The latter is holistic as the phase of the warping frequency is estimated by taking a weighted combination of **all** warping sources. More details on when to use which? Check out "[which phase should I warp to, FFT or GED?](#which-phase-should-i-warp-to-FFT-or-GED)".

After [bt_clocktobrain](bt_clocktobrain/bt_clocktobrain.m), your data has completed its transformation from clock to brain time. The time axis is now formatted as sequences of cycles, instead of seconds. You may opt to continue analyses outside of the toolbox, or test for periodic patterns in the data using `braintime`'s second operation described next. In case of the former, please read the previously linked circularity information.


## Operation 2: Periodicity analysis

You may wish to use `braintime` to test whether brain time warping has accentuated dynamic patterns in the data. This is appropriate if your data comprises two underlying conditions which are predicted to yield fluctuating classification evidence. For example, we report a spatial attention dataset where we predict that evidence for motion direction rapidly alternates from low to high as a function of alpha phase (the warping signal). Thus, the second operation requires your brain time warped data to consist of such classes.

Before getting started, ensure the original clock time data as well as the brain time warped data have their classes labeled. Specifically, create a field called "clabel", with a vector of 1's and 2's corresponding to the condition of each trial. For example, if your trial sequence is 'left left right' conditions, the clabel field should contain ```[1 1 2]```.

The pipeline for this toolbox operation is that each participant's data needs to undergo step 2.1 to 2.3, separately for both clock and brain time data. Then, each participant's output needs to be added as a field to a large group structure, separately for clock and brain time results. These two group structures (e.g. ```ct_stats1``` and ```bt_stats1```) can each be used as input for step 2.4 to obtain group level results, allowing you to compare periodicity in clock time and brain time data.

> :warning: A [recent paper](https://doi.org/10.1016/j.jneumeth.2021.109080) by van Driel et al. (2021) shows that high pass filtering data may cause artefacts through displacement of information. We have found that these artefacts may sometimes look periodic. Thus, when using operation 2 of the toolbox, please ensure no high pass filter was applied to the data, or do so with extreme care. van Driel et al. introduce alternative ways of getting rid of low frequencies (robust detrending using a trial-based mask), and make several suggestions.

**2.1 Multivariate pattern analysis**

The first step is to use `MVPA-Light` to classify the data. You may opt to classify across time [mv_classify_across_time.m](https://github.com/treder/MVPA-Light/blob/master/mv_classify_across_time.m). This method tests whether a classifier can separate both classes of data across time in the trials. Alternatively, you can apply temporal generalization by classifying using [mv_classify_timextime.m](https://github.com/treder/MVPA-Light/blob/master/mv_classify_timextime.m). This method tests for temporal generalization of classification. That is, it tests to what extent a classifier trained on one timepoint generalizes its performance to other timepoints. Not sure whether to classify across time, or test for temporal generalization? Check out "[should I classify across time or test for temporal generalization?](#should-i-classify-across-time-or-test-for-temporal-generalization").

> :bulb: You can change a variety of parameters when calling `MVPA-Light`, described briefly in [tutorial 2](tutorial/tutorial2_periodicity.m). For more information, check `MVPA-Light`'s awesome [tutorials](https://github.com/treder/MVPA-Light/tree/master/examples).

**2.2 Quantify periodicity in classifier performance**

If the warping signal orchestrates the dynamics of your cognitive function, operationalized by your two classes of data, then the classifier's performance may tap into these dynamics. Thus, [bt_quantify](periodicity/bt_quantify.m) tests and quantifies periodic patterns in the classifier performance, whether that is in the classifier output obtained from [mv_classify_across_time.m](https://github.com/treder/MVPA-Light/blob/master/mv_classify_across_time.m) (1 dimensional; "diagonal"), or [mv_classify_timextime.m](https://github.com/treder/MVPA-Light/blob/master/mv_classify_timextime.m) (2 dimensional; "temporal generalization matrix" (TGM)). In both cases, you get a "periodicity spectrum" which is an average power spectrum of all the data in your 1D or 2D classifier performance.

For TGMs, you may opt to perform the periodicity analysis over either the 2 dimensional matrix itself (```cfg.method = 'tgm'```) or its autocorrelation map (```cfg.method = 'ac'```). Not sure what is better? Check out "[should I perform bt_quantify over the TGM itself, or its autocorrelation map?](#should-i-perform-bt_quantify-over-the-tgm-itself-or-its-autocorrelation-map)".

You also need to specify a range of periodicity frequencies. In which range of rates do you expect periodic patterns in classifier performance to arise? Finally, you can choose a reference dimension. If you choose ```cfg.refdimension = clocktime```, periodicity in classifier performance will be displayed as a function of cycles per seconds (Hz). Alternatively, if you are quantifying brain time warped data, ```cfg.refdimension = clocktime``` is more appropriate. This references the peroidicity as a function of cycles per second, *normalized* to that participants' warping frequency.

> :bulb: Let's say you warp a participant's data to an 11 Hz warping signal. In that case, with ```cfg.refdimension = clocktime```, a peak at the warping frequency will be at 11 Hz, but with ```cfg.refdimension = clocktime``` it will show up at 1 Hz, as the frequencies are normalized to the warped frequency (11/11 = 1).

**2.3 First level statistics**

How does the participant's quantified periodicity compare against the null distribution? [bt_statslevel1](periodicity/bt_statslevel1.m) takes the output from [bt_quantify](periodicity/bt_quantify.m) and performs classification ```cfg.numperms1 = n1``` times over, each time randomly shuffling the classification labels and obtaining ```n1``` "permuted" periodicity spectra. This provides a null distribution that quantifies how much periodicity is in the data when the class structure is destroyed, setting things up for p-value estimation on the group level.

Now, repeat [bt_statslevel1](periodicity/bt_statslevel1.m) for every participant, and send its output to a separate field in a group structure (e.g. ```[ct_stats1{subj}] = bt_statslevel1(cfg,data,quant)```). The next step requires this format to perform 2nd level statistics — it loops over the fields.

**2.4 Second level statistics**

We just got an idea of the periodicity in single participants. But what about the group level? [bt_statslevel2](periodicity/bt_statslevel2.m) employs the following algorithm to generate second level statistics:

1) ```cfg.numperms2 = n2``` times, do the following: randomly grab one of each participants' ```n1``` permuted periodicity spectra and average those spectra into one spectrum. This yields ```n2``` spectra.
2) For each frequency, compare the empirical periodicity with the distribution of permuted periodicity to derive a p-value.
3) Plot the average empirical periodicity spectrum and display the p-value for each frequency.

> :bulb: In brain time results, p-values at 0.5 Hz, 1 Hz, and 2 Hz are exempted from multiple testing correction, as periodicity is predicted at one or multiple of these rates depending on the underlying structure'.

> :bulb: `braintime` z-scores participants' spectra to enable second level statistics. For details on how it does this, check out "[How does the toolbox z-score periodicity spectra?](#how-does-the-toolbox-z-score-periodicity-spectra)".

At this point we have tested whether periodicity in classification performance is significant for clock and brain time. However, a separate question is whether classification performance *in itself* is significant, without any question of periodic patterns. [bt_statslevel2](periodicity/bt_statslevel2.m) also tests for this, by calling `MVPA-Light`. For this, the function requires a separate structure called ```cfg_clus```, where you can enter a range of parameters for cluster correction. So, aside from the periodicity results, you also get a display of which datapoints show above chance classification performance.

## Toolbox considerations

### How do I know all this even works?
We provide three sources of evidence for the toolbox. A simulated dataset, a rodent intracranial dataset, and a human EEG dataset. `braintime` includes the script used to generate the simulated dataset. Toying around with it is a great way to both get a feel for the toolbox, and see its effects. Check out [bt_dipsim](dipolesimulation/bt_dipsim.m), change some parameters to your liking, and run it through the `braintime` pipeline to see that the toolbox accounts for clock and brain time disharmony. Details on the rodent and human evidence are described in the [main manuscript](https://www.biorxiv.org/content/10.1101/2021.06.09.447763v1).

### Which cutmethod to choose?
At the start of trials, `braintime` repeats original data samples until the phase of the warping signal aligns with the stationary signal (see paragraph 1.5). This takes a while, depending on their disharmony. This data repetition may cause an artefact in further analyses, called a "first cycle artefact". In `braintime`'s _periodicity analysis_ operation for example, you can sometimes see a stretched out pattern of low or high classification. Most often, this artefact is very small, and unlikely to alter subsequent analyses. Hence, `braintime`'s default option under [bt_analyzesources](warpingsource/bt_analyzesources.m) is ```cfg.cutmethod = 'consistenttime'```. This default method is called ```consistenttime``` because the toolbox warps exactly to the specified time window of interest.

To the extent the first cycle artefact is predicted to be an issue for further analyses, you may decide to remove it by selecting ```cfg.cutmethod = 'cutartefact'``` under [bt_analyzesources](warpingsource/bt_analyzesources.m). With this parameter, `brainime` warps a little bit of time before and after the specified time window of interest before it is removed in later steps. This makes it so the data reptition segment is removed from the final data. The downside to this method is that the warped data will comprise more (or less) than the specified time window of interest.

Finally, there's one more thing to consider here. `braintime` allocates data samples to cycles depending on the warping path. As the two cutmehods comprise warping to a different length of data, there may be variations in which data sample is considered part of which cycle. To check how each method allocates data to cycles, check out [tutorial 7](tutorial/tutorial7_checkallocation.m).

In short, the upside of ```cfg.cutmethod = 'consistenttime'``` is that the warped data will exactly of the specified time window, but its downside is a first cycle artefact. In contrast, the upside of ```cfg.cutmethod = 'curartefact'``` is that the warped data will not contain a first cycle artefact, but its downside is that the data's time may not exactly match the specified window. The methods may vary in which data they allocate to which cycle.

### Is it circular to warp to warping sources obtained from my clock time data?

It depends on subsequent analyses. First, let's consider the concern. The concern is that warping data to the phase dynamics of a warping oscillation trivially introduces patterns in the data oscillating at the warping oscillation's frequency, as you've changed the data according to it.

First, this is only a concern for warping sources obtained from clock time data — if the phase dynamics used for warping are not in the data on which subsequent analyses are done, there's no methodological circularity. Second, even if clock time data and warping sources are dependent, subsequent periodicity analyses _within_ `braintime` are safe, as the analysis inherently accounts for cicularity. Namely, periodicity in classification performance is tested by shuffling class labels, which tests for the relevance of the class structure to the classifier. If the warping trivially introduces periodicity in the classifier, it should do so equally for shuffled results (permuted periodicity spectra) and correctly labeled results (empirical periodicity spectra).

Critically however, if you intend to only use `braintime`'s first function, brain time warping, and continue analyses on the warped data outside the toolbox, it is strongly recommended to remove the warping source from the original clock time data (for example, by removing the independent component or by removing the channel).


### Which phase should I warp to, FFT or GED?

FFT (```cfg.phasemethod = 'fft'```) is more constricted to the warping source of interest, while GED (```cfg.phasemethod = 'ged'```) takes a weighted approach across all warping sources. So the best choice depends on whether you expect the oscillations of interest to be isolated or spread out across the data. To illustrate the best choice depending on circumstance let's look at some examples. If we are studying attention in the parietal cortex, or motor processes in the motor cortex, it may be better to take FFT on the warping source that shows topographies specific to our process of interest. Here, you don't want oscillations from other regions, so GED should be avoided. In contrast, if your warping sources are already restricted to a region, like the local field potential of hippocmapal regions or source localized parietal data, then a holistic estimation using GED is more appropriate as all warping sources are assumed to contain relevant phase dynamics.

### Should I perform bt_quantify over the TGM itself, or its autocorrelation map?

It depends on the uniformity of the periodic patterns of your 2D classifier performance (TGM). In principle, the autocorrelation map (```cfg.method = 'ac'```) is a more powerful approach, as it is known to accentuate periodic patterns in the data. However, when multiple periodicity rates are present in the TGM, the autocorrelation tends to pick up the strongest rate and drown out the weaker one. Similarly, when periodicity is present in only a small part of the TGM, autocorrelation map may not be able to accentuate its pattern. In these cases, it's better to analyze the TGM itself (```cfg.method = 'tgm'```).

Thus, when you expect uniform periodicity patterns in the data, analyzing the autocorrelation map of TGMs is a powerful approach. When the pattern is expected to only be present partially, or when multiple rates are predicted, it is better to analyze the TGM itself.

### How does the toolbox z-score periodicity spectra and implement statistics?

We want to prevent that arbitrary differences in the periodicity power spectra between participants drive a group effect. To this end, `braintime` z-scores each participant's empirical and permuted periodicity spectrum in the following way:

| Step  | *first level* stats ([bt_statslevel1](periodicity/bt_statslevel1.m)) and z-scoring algorithm |
| --- | :-- |
| 1 | Z-score each of ```n1``` permuted periodicity spectra by the mean and standard deviation of that spectrum. |
| 2 | Z-score the single empirical dsitribution by the mean and standard deviation of all ```n1``` permutations. |

| Step  | *second level* stats ([bt_statslevel2](periodicity/bt_statslevel2.m)) |
| --- | :-- |
| 1  | for each second level permutation ```n2```: Grab one (random) permuted spectrum of each participant's ```n1``` permuted spectra and average them into one. This yields ```n2``` permuted averages. |
| 2  | for each frequency, make a p-value from the percentile at which your empirical power lands in the distribution of ```n2``` permuted power values. |

### I have another question!
If it's a technical question about the toolbox, it's likely explained in the [main paper](https://www.biorxiv.org/content/10.1101/2021.06.09.447763v1)'s supplementary material. For other questions, feel free to contact me at sandervanbree@gmail.com and I hope to have some time to answer.
