% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KrakenAnnotateStatistics.R
\name{kraken_info_standardised_metrics}
\alias{kraken_info_standardised_metrics}
\title{Adding Standardised Metrics}
\usage{
kraken_info_standardised_metrics()
}
\description{
\tabular{lll}{
   \strong{Metric} \tab \strong{Description} \tab \strong{Function} \cr
   Zscore \tab Centering = Mean, Scaling = Stdev \tab \code{\link[=kraken_report_add_zscore]{kraken_report_add_zscore()}} \cr
   ZscoreRobust \tab Centering = Median, Scaling = MAD \tab \code{\link[=kraken_report_add_robust_zscore]{kraken_report_add_robust_zscore()}} \cr
   MedianIQR \tab Centering = Median, Scaling = IQR \tab \code{\link[=kraken_report_add_median_centered_iqr_scaled_metric]{kraken_report_add_median_centered_iqr_scaled_metric()}} \cr
   MedianQn \tab Centering = Median, Scaling = Qn \tab \code{\link[=kraken_report_add_median_centered_qn_scaled_metric]{kraken_report_add_median_centered_qn_scaled_metric()}} \cr
   MedianTau2 \tab Centering = Median, Scaling = Tau2 \tab \code{\link[=kraken_report_add_median_centered_tau2_scaled_metric]{kraken_report_add_median_centered_tau2_scaled_metric()}} \cr
}
}
\details{
Run any of the functions listed in description on your report dataframe (generated from kraken_reports_parse or another parsing function provided by this package).
The relevant statistic will be added to the data.table (modified in place, so don't need to save output of functions to variables).

A 'loggable' version of each metric will also be added (instead of centering around zero - centered around a value high enough to ensure even the smallest values are > 0 and therefore loggable)
}
\seealso{
Other standardised metrics: 
\code{\link{kraken_report_add_median_centered_iqr_scaled_metric}()},
\code{\link{kraken_report_add_median_centered_qn_scaled_metric}()},
\code{\link{kraken_report_add_median_centered_tau2_scaled_metric}()},
\code{\link{kraken_report_add_robust_zscore}()},
\code{\link{kraken_report_add_zscore}()}
}
\concept{standardised metrics}
