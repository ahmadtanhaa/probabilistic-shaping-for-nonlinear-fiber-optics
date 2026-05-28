# Probabilistic Shaping for Nonlinear Fiber Optic Channels

[![IEEE](https://img.shields.io/badge/IEEE-9827181-00629B.svg)](https://ieeexplore.ieee.org/document/9827181/)
[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

## MATLAB Codes for the paper "On the Impact of Probabilistic Shaping on Fiber Nonlinearities".

This repository contains the MATLAB simulation codes used to reproduce the main numerical results of the paper:

**On the Impact of Probabilistic Shaping on Fiber Nonlinearities**  
Ahmad Tanha, Hami Rabbani, and Lotfollah Beygi  

The paper was presented at the 2022 30th International Conference on Electrical Engineering (ICEE).

The paper is available online on IEEE Xplore:

[https://ieeexplore.ieee.org/document/9827181/](https://ieeexplore.ieee.org/document/9827181/)

## Overview

This work studies the impact of probabilistic shaping (PS) on fiber nonlinearities in a polarization-multiplexed 16-QAM wavelength-division multiplexed (WDM) optical fiber system.

The simulations are based on the enhanced Gaussian noise (EGN) model and evaluate the main nonlinear interference (NLI) components, including self-channel interference (SCI), cross-channel interference (XCI), multi-channel interference (MCI), total NLI, effective signal-to-noise ratio (SNR), and symbol error rate (SER).

Two probabilistic shaping scenarios are considered:

1. Probabilistic shaping applied only to the channel of interest (COI).
2. Probabilistic shaping applied to all WDM channels.

The codes in this repository reproduce the numerical plots reported in the paper.

## Contents

The repository includes MATLAB routines for:

- Computing EGN-model interaction coefficients using Monte-Carlo integration,
- Evaluating SCI, XCI, MCI, and total NLI,
- Computing effective SNR over transmission distance,
- Computing SER for uniform and probabilistically shaped 16-QAM,
- Reproducing the results for PS applied only to the channel of interest,
- Reproducing the results for PS applied to all WDM channels,
- Generating paper-style MATLAB figures for Figs. 2--5.

## Files

- `calc_interChannel.m`  
  Computes the EGN-model interaction coefficients using Monte-Carlo integration.

- `run_scenario1_PS_in_COI.m`  
  Reproduces Figs. 2 and 3 of the paper, where probabilistic shaping is applied only to the channel of interest.

- `run_scenario2_PS_all_WDM.m`  
  Reproduces Figs. 4 and 5 of the paper, where probabilistic shaping is applied to all WDM channels.


## Requirements

The codes are written in MATLAB and use standard MATLAB functions.

No additional toolbox is required for the main simulations. A local implementation of the Q-function is included in the scripts for compatibility.

## How to Run

To reproduce the results of Scenario I, where probabilistic shaping is applied only to the channel of interest, run:

```matlab
run_scenario1_PS_in_COI
````

To reproduce the results of Scenario II, where probabilistic shaping is applied to all WDM channels, run:

```matlab
run_scenario2_PS_all_WDM
```

The generated figures are saved automatically in the `results/` folder in both `.fig` and `.png` formats.

## Note

The scripts use Monte-Carlo integration. In the original research code, the number of Monte-Carlo samples is set to:

```matlab
params.N = 1e6;
```

For a quick test, this value may be temporarily reduced, for example:

```matlab
params.N = 1e4;
```

However, reducing the number of Monte-Carlo samples may make the curves less accurate.

## Citation

If you use this code or find it useful in your research, please cite our paper:

```bibtex
@inproceedings{tanha2022impact,
  title={On the Impact of Probabilistic Shaping on Fiber Nonlinearities},
  author={Tanha, Ahmad and Rabbani, Hami and Beygi, Lotfollah},
  booktitle={Proceedings, IEEE 30th International Conference on Electrical Engineering (ICEE)},
  year={2022},
  organization={IEEE}
}
```

The paper is available at:

[https://ieeexplore.ieee.org/document/9827181/](https://ieeexplore.ieee.org/document/9827181/)

## License

This work is licensed under the **Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License**, also known as **CC BY-NC-ND 4.0**.

You are free to share this material, provided that appropriate credit is given, a link to the license is provided, and any changes are indicated.

You may not use the material for commercial purposes.

You may not distribute modified versions of the material.

For more details, see:

[https://creativecommons.org/licenses/by-nc-nd/4.0/](https://creativecommons.org/licenses/by-nc-nd/4.0/)

## Disclaimer

The codes are provided for research and reproducibility purposes. They are intended to accompany the paper and reproduce the numerical results presented therein.
