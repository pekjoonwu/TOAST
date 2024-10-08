% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TOAST.R
\name{three_outcomes_count}
\alias{three_outcomes_count}
\title{TOAST for count type of endpoints}
\usage{
three_outcomes_count(
  shape_1,
  shape_2,
  rate_1,
  rate_2,
  Method = NULL,
  nrep,
  NREP,
  sample_size,
  LRV,
  TV,
  FGR = 0.2,
  FSR = 0.1,
  seed
)
}
\arguments{
\item{shape_1}{Numeric. Shape parameter to generate Poisson rate from control group in the previous studies}

\item{shape_2}{Numeric. Shape parameter to generate Poisson rate from treatment group in the previous studies}

\item{rate_1}{Numeric. Rate parameter to generate Poisson rate from control group in the previous studies}

\item{rate_2}{Numeric. Rate parameter to generate Poisson rate from treatment group in the previous studies}

\item{Method}{Character. Method chosen for confidence interval calculation - empirical or bootstrap}

\item{nrep}{Numeric. Replicates to evaluate the decision once}

\item{NREP}{Numeric. Replicates to calculate the decision proportions}

\item{sample_size}{Numeric. Target Sample size in the new studies}

\item{LRV}{Numeric. Lower reference value (treatment difference) defined prior to the studies}

\item{TV}{Numeric. Target value (treatment difference) defined prior to the studies}

\item{FGR}{Numeric. False go risk level}

\item{FSR}{Numeric. False stop risk level}

\item{seed}{Numeric. Simulation seed for replication purpose}
}
\value{
A list containing the population parameters/ TPP values/ number of simulations/ decision results
}
\description{
TOAST for count type of endpoints
}
\examples{
library(TOAST)
## Example code for one count primary endpoint
count_examples <- three_outcomes_count(shape_1 = 7, shape_2 = 4, rate_1 = 8.8, rate_2 = 2,
   nrep = 500, NREP = 1000, sample_size = 100,
   LRV = 1.0, Method = "empirical",
   FGR = 0.2, FSR = 0.1, TV = 1.2, seed = 2024)
}
