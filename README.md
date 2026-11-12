# Hook depth and species interactions in the American Samoa longline fishery

This repository contains R code for a depth-based hazard analysis of the Amendment 5 gear requirements in the American Samoa longline fishery. Since September 2011, these requirements have kept baited hooks deeper than about 100 m, mainly through 70 m of mainline without branch lines next to each float and float lines at least 30 m long, to reduce interactions with green sea turtles.

The analysis asks two questions:

1. **Retrospective:** how did interaction rates change after the requirements took effect, and how much of that change came from fishing hooks deeper?
2. **Prospective:** how would interaction rates change if the requirements were removed?

Five species are analyzed: green sea turtle, oceanic whitetip shark, mahimahi, wahoo, and albacore.

## How it works

The depth of every hook in every set is estimated from the recorded gear with a catenary equation, with and without an adjustment for shoaling. For each species, a generalized additive model (`mgcv::bam`) estimates the interaction rate per hook as a smooth function of hook depth. The depth function is shared between periods, so the gear change acts only through where hooks are set. The model also includes a period effect and trip and year random effects. The fitted depth functions are then used to:

- split the observed pre/post change into a gear component and a period component
- project interaction rates under four gear scenarios if the requirements were removed

## Data

The analysis was run on confidential fishery observer data, which are **not** included in this repository.

**If you are not authorized to use the observer data**, generate the simulated dataset and run the analysis on it:

```bash
Rscript R/make_sim_data.R     # writes data/sim_as_data.Rdata
Rscript R/asll_dereg_risk.R   # uses the simulated data by default
```

The simulated data have the same structure and columns as the observer data, with fictitious trip, set, and permit identifiers. Catches are generated from known depth relationships with no period effect, so the code should recover period effects close to one. A few sets have deliberately corrupted gear records so the screening step has something to flag.

The simulated data are for running and testing the code only. They are not calibrated to the fishery, and results from them do not reproduce the results of the analysis and must not be interpreted as estimates for the American Samoa longline fishery.

**If you are authorized to use the observer data**, keep the file outside the repository and point the analysis to it, for example in `~/.Renviron`:

```
ASLL_DATA_FILE=/secure/path/to/as_data_depth.Rdata
```

Do not commit the observer data, or set- or trip-level outputs derived from them (such as `flagged_sets.csv`), to this or any public repository.

## Requirements

R (≥ 4.1) with `dplyr`, `tidyr`, `purrr`, `mgcv`, `ggplot2`, and `here`. Results are written to `output/dereg_risk/`.

## Contact

[Robert Ahrens, NOAA-PIFSC, Robert.Ahrens@noaa.gov]

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
