# TOAST
## Installation
To install the latest version of the TOAST package from GitHub, run the following code in R:
```
install.packages('devtools')
library(devtools)
devtools::install_github('pekjoonwu/TOAST')
```
This command should automatically install all required packages if they are not installed already.

## Quick Start
Our framework includes two parts: Three outcome decision-making for one primary endpoint and for two co-primary endpoints.
1. For one primary endpoint, please use these following commands to see the example usage.
```
library(TOAST)
?three_outcomes_cont ## For continuous endpoint
?three_outcomes_binary ## For binary endpoint
?three_outcomes_count ## For count endpoint using Poisson distribution
?three_outcomes_count_negbin ## For count endpoint using Negative binomial distribution
?three_outcomes_survival ## For time-to-event endpoint
```

2. For co-primary endpoint, please use these following commands to see the example usage.
```
library(TOAST)
?three_outcomes_cont_coprim ## For continuous + continuous endpoints
?three_outcomes_binary_coprim ## For binary + binary endpoints
?three_outcomes_count_coprim ## For count + count endpoints
?three_outcomes_cont_binary_coprimary ## For count + binary endpoints
?three_outcomes_cont_count_coprimary ## For count + count endpoints
?three_outcomes_bin_count_coprimary ## For binary + count endpoints
```
