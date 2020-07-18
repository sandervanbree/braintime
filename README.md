# Brain Time Toolbox

![btlogo](https://i.imgur.com/cjhrUnt.png)

Warp electrophysiological data from clock time (default) to brain time and analyze the rate of recurrence of mental representations. By *Sander van Bree, María Melcón, Maria Wimber & Simon Hanslmayr*.

### Dependencies:
- FieldTrip
- MVPA Light
- MATLAB Signal Processing Toolbox
- [Several functions](dependencies) included in the toolbox

### Glossary:
| Term | Description |
| --- | --- |
| Clock time (CT) | data ordered by sequences of seconds |
| Brain time (BT) | data ordered by sequences of neural oscillations |
| Warping | a procedure to nonlinearly rescale time series |
| Carrier | an oscillation to which clock time data can be warped |
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

## Introduction
