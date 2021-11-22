Quantifying the probability of detection of wild ungulates using the
‘Judas’ technique
================

## Overview

This repository contains data and code from:

Ramsey, D.S.L., Campbell, K.J., Lavoie, C., MacDonald, N., and Morrison,
S.A. (2022). “Quantifying the probability of detection of wild ungulates
using the ‘Judas’ technique” *Conservation Biology*

### File descriptions:

-   `discrete_survival_analysis.r` performs analysis of association rate
    between Judas dyads using `stan`
-   `discrete_survival_nimble.r` performs analysis of association rate
    using `nimble`
-   `judas_erad_prob.r` Uses estimates of association rate by Judas to
    make inferences about the probability of eradication, given no
    detection of wild individuals, for an area of interest
-   `judas_functions.r` contains various functions required by the above
    scripts

## Prerequisites

The simulation scripts require packages `tidyverse`, `cmdstan`,
`posterior`, `lubridate`, `sf`, `terra`, `simecol`,`ggspatial` and
`gridExtra`.

## Getting started

Data from the manuscript for both Isla Santiago (Judas Goats) and Santa
Cruz Island (Judas Pigs) are available in the `Data/` directory. These
contain data on observed associations and contemporaneous exposure times
for each Judas dyad `*_judas_dyads.csv`, infomation on Judas IDs
including entry and exit times `*_judas_id.csv` and locations of Judas
individuals during routine monitoring `*_judas_locs.csv`.

Start by stepping through the `discrete_survival_analysis.r` script to
generate estimates of the association probabilities for each Judas
individual. Save these to the `out/` directory using the commands
provided. Scripts to produce the various plots in the manuscript are
also provided.

Once estimates are generated, these can then be used to estimate
probabilities of eradication using `judas_erad_prob.r`. Shapefiles of
both Isla Santiago and Santa Cruz Island (zone 1) are provided in the
`Data/` directory for this purpose.
