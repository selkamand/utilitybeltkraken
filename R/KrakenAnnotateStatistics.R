
# Basic Stat Functions ----------------------------------------------------
Zscore <- function(vec){
  z <- (vec-mean(vec))/sd(vec)
  z <- ifelse(is.na(z), yes = 0, no = z) #If denominator is zero AND numerator is zero
  return(z)
}

ZscoreRobust <- function(vec){
  z <- (vec-median(vec))/mad(vec)
  z <- ifelse(is.na(z), yes = 0, no = z)
  return(z)
}

MedianIQR <- function(vec){
  i <- (vec-median(vec))/IQR(vec)
  i <- ifelse(is.na(i), yes = 0, no = i)
  return(i)
}

MedianQn <- function(vec){
  z <- (vec - median(vec))/robustbase::Qn(vec)
  z <- ifelse(is.na(z), yes = 0, no = z)
  return(z)
}

MedianTau2 <- function(vec){
  z <- (vec - median(vec))/robustbase::scaleTau2(vec)
  z <- ifelse(is.na(z), yes = 0, no = z)
  return(z)
}

MakeLoggable <- function(vec){
  return(vec + abs(min(vec)) + 1)
}

# Stats -------------------------------------------------------------------

#' Adding Standardised Metrics
#'
#' | **Metric**   	| **Description**                    	| **Function**                                           	|
#' |--------------	|------------------------------------	|--------------------------------------------------------	|
#' | Zscore       	| Centering = Mean, Scaling = Stdev  	| [kraken_report_add_zscore()]                             	|
#' | ZscoreRobust 	| Centering = Median, Scaling = MAD  	| [kraken_report_add_robust_zscore()]                      	|
#' | MedianIQR    	| Centering = Median, Scaling = IQR  	| [kraken_report_add_median_centered_iqr_scaled_metric()]  	|
#' | MedianQn     	| Centering = Median, Scaling = Qn   	| [kraken_report_add_median_centered_qn_scaled_metric()]   	|
#' | MedianTau2   	| Centering = Median, Scaling = Tau2 	| [kraken_report_add_median_centered_tau2_scaled_metric()] 	|
#'
#' Run any of the functions listed in description on your report dataframe (generated from kraken_reports_parse or another parsing function provided by this package).
#' The relevant statistic will be added to the data.table (modified in place, so don't need to save output of functions to variables).
#'
#' A 'loggable' version of each metric will also be added (instead of centering around zero - centered around a value high enough to ensure even the smallest values are > 0 and therefore loggable)
#'
#' @return NULL
#' @export
#' @family standardised metrics
#'
kraken_info_standardised_metrics <- function(){
  message("Please Run: \n\t?kraken_info_standardised_metrics\n")
  return(invisible(NULL))
}

#' Add Zscore to Report
#'
#' @param kraken_report_df the kraken dataframe produced by \link{kraken_reports_parse} or link{kraken_report_parse}
#'
#' @return Nothing - input data.table is modified in place
#' @export
#' @family standardised metrics
kraken_report_add_zscore <- function(kraken_report_df){
  message("Adding columns: [Zscore & ZscoreLoggable] to kraken report datatable")
  kraken_report_df[, `:=`(Zscore = Zscore(RPM)), by = .(TaxonomyID)]
  kraken_report_df[, `:=`(ZscoreLoggable = MakeLoggable(Zscore))]
  return(invisible(NULL))
}

#' Add Robust Zscore to Report
#'
#' Adds a MAD derived robust zscore to data.table (modifies in place)
#'
#' @inheritParams kraken_report_add_zscore
#' @inherit kraken_report_add_zscore return
#' @export
#' @family standardised metrics
kraken_report_add_robust_zscore <- function(kraken_report_df){
  message("Adding columns: [ZscoreRobust & ZscoreRobustLoggable] to kraken report datatable")
  kraken_report_df[, `:=`(ZscoreRobust = ZscoreRobust(RPM)), by = .(TaxonomyID)]
  kraken_report_df[, `:=`(ZscoreRobustLoggable = MakeLoggable(ZscoreRobust))]

  # browser()
  # kraken_report_df %>%
  #   dtplyr::lazy_dt(immutable = FALSE) %>%
  #   dplyr::group_by(TaxonomyID) %>%
  #     ZscoreRobust =  ZscoreRobust(RPM),
  #     ZscoreRobustLoggable = MakeLoggable(ZscoreRobust)
  #     ) %>%
  #   dplyr::ungroup() %>%
  #   data.table::as.data.table()
  return(invisible(NULL))
}


#' Add Robust Zscore to Report
#'
#' Adds a add median centered IQR scaled  metric to data.table (modifies in place)
#'
#' @inheritParams kraken_report_add_zscore
#' @inherit kraken_report_add_zscore return
#' @export
#' @family standardised metrics
kraken_report_add_median_centered_iqr_scaled_metric <- function(kraken_report_df){
  message("Adding columns: [MedianIQR & MedianIQRLoggable] to kraken report datatable")
  kraken_report_df[, `:=`(MedianIQR = MedianIQR(RPM)), by = .(TaxonomyID)]
  kraken_report_df[, `:=`(MedianIQRLoggable = MakeLoggable(MedianIQR))]
  return(invisible(NULL))
}

#' Add MedianQn to Report
#'
#' Adds a median centered qn scaled metric to data.table (modifies in place)
#'
#' @inheritParams kraken_report_add_zscore
#' @inherit kraken_report_add_zscore return
#' @export
#' @family standardised metrics
kraken_report_add_median_centered_qn_scaled_metric <- function(kraken_report_df){
  message("Adding columns: [MedianQn & MedianQnLoggable] to kraken report datatable")
  kraken_report_df[, `:=`(MedianQn = MedianQn(RPM)), by = .(TaxonomyID)]
  kraken_report_df[, `:=`(MedianQnLoggable = MakeLoggable(MedianQn)), by = .(TaxonomyID)]
  return(invisible(NULL))
}


#' Add MedianTau2 to Report
#'
#' Adds a median centered Tau2 scaled metric to data.table (modifies in place)
#'
#' @inheritParams kraken_report_add_zscore
#' @inherit kraken_report_add_zscore return
#' @export
#' @family standardised metrics
kraken_report_add_median_centered_tau2_scaled_metric <- function(kraken_report_df){
  message("Adding columns: [MedianTau2 & MedianTau2Loggable] to kraken report datatable")
  kraken_report_df[, `:=`(MedianTau2 = MedianTau2(RPM)), by = .(TaxonomyID)]
  kraken_report_df[, `:=`(MedianTau2Loggable = MakeLoggable(MedianTau2))]
  return(invisible(NULL))
}





