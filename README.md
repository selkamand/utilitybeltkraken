
<!-- README.md is generated from README.Rmd. Please edit that file -->

# utilitybeltkraken

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of utilitybeltkraken is to simplify cohort-analysis of sample
analysed with the kraken2 taxonomic classification system

## Installation

You can install the development version of utilitybeltkraken from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("selkamand/utilitybeltkraken")
```

## Usage

``` r
# Path to folder containing kraken reports 
# (reports must be named with sample_id as prefix)
kreport_dir <- system.file(package="utilitybeltkraken", "simulated_data/simulated_kraken_reports_inc_zero_counts/")

# Create kraken sqlite database
kraken_db <- kraken_reports_parse_to_sqlite_db(kreport_dir)
```
