### Current to do list (for devs) | Official documentation starts below
- random gaussian field theory to perform multiple testing correction over power spectra and TGMs?

LOW URGENCY:
- Fix dots in the asymmetry plot in bt_selectsource
- Add licenses! 
- Update tutorial 4 to 6
- In the paper, remind users to be very careful with high pass filtering; https://doi.org/10.1101/530220

COMPLETED RECENTLY: 
- Implement Cohen/Luca's GED scripts
- For diagonal testing, test significance of diagonal, not whole TGM
- make clear you don't have to get components
- Implement testing of the diagonal only
- Statistical masking requires too much memory
- spectral resolution (alignment warped freq and peak)
- Make lingo consistent within toolbox, in line with SfN poster ("warping signal", "warping source", etc.)
- Resize second level power spectra to the right number of frequency bins, discuss together what that is
- Separately, decide for "shuffled" versus "permuted" and make consistent

### Start documentation:


# Brain Time Toolbox
<img src="https://i.imgur.com/OAEVqgM.png" width="200">

Warp electrophysiological data from clock time (default) to brain time and analyze the recurrence of mental representations. By *Sander van Bree, María Melcón, Maria Wimber & Simon Hanslmayr*.

### Dependencies
- [FieldTrip](http://www.fieldtriptoolbox.org/download/)
- [MVPA Light](https://github.com/treder/MVPA-Light)
- MATLAB Signal Processing Toolbox
- [Several functions](dependencies) included in the toolbox

### Glossary
| Term | Description |
| --- | --- |
| Clock time (CT) | data ordered by sequences of seconds |
| Brain time (BT) | data ordered by the phase of a carrier oscillation |
| Warping source | a data structure that contains the to-be-selected warping signal |
| Warping signal | the signal within your warping source that shows optimal time frequency characteristics |
| Brain time warping | a procedure to nonlinearly rescale clock time data based on the phase dynamics of the warping signal |
| Recurrence | rhythmic patterns of classification performance emerging from data trained on previous time points |
| Temporal Generalization Matrix (TGM): | testing and training a classifier at all time points |

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
