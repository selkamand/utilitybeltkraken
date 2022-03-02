
#' Visualise Distribution
#'
#' @inheritParams kraken_report_add_zscore
#' @param taxids_of_interest
#'
#' @return
#' @export
#'
kraken_report_visualise_distributions <- function(kraken_report_df, taxids_of_interest, metric = "ZscoreRobust" ,show_sample_vlines = FALSE, supporting_reads_threshold = 50, sample_of_interest = NA, use_loggable_version = TRUE, seed = 1, pointsize = c(1, 2), ...){

  # Defensive Assertions
  standardised_metrics = c("Zscore", "ZscoreRobust", "MedianIQR", "MedianQn", "MedianTau2")
  assertthat::assert_that(assertthat::is.string(metric))
  assertthat::assert_that(metric %in% standardised_metrics, msg = paste0("metric (", metric,") must be one of: [", paste0(standardised_metrics, collapse = ", "), "]"))
  assertthat::assert_that(is.data.frame(kraken_report_df))
  if(use_loggable_version) metric <- paste0(metric, "Loggable")
  assertthat::assert_that(metric %in% colnames(kraken_report_df), msg = paste0("kraken_report_df must contain the column: [", metric, "]. Please see ?kraken_info_standardised_metrics for information on how to add standard metrics to your dataframe."))
  columns = c("SampleID", "TaxonomyID", "ScientificName", "ReadsCoveredByClade", "RPM", metric, "ZscoreRobust")
  assertthat::assert_that(all(columns %in% colnames(kraken_report_df)), msg = paste0("kraken_report_df must contain the columns [", paste0(columns, collapse=", "), "]"))
  assertthat::assert_that(length(pointsize) == 2, msg = "pointsize must a vector of length 2. First number is size of most samples and second number for size of 'samples of interest' ")
  #assertthat::assert_that(color_aesthetic %in% colnames(kraken_report_df), msg = paste0("kraken_report_df must contain the column [", paste0(color_aesthetic, collapse=", "), "], to use it as a color aesthetic"))


  # Subset the data
  subset_df = kraken_report_df[, (columns), with=FALSE][TaxonomyID %in% taxids_of_interest][, `:=` (Infinite = is.infinite(get(metric)))]
  subset_df = subset_df[, `:=`(SampleOfInterest = SampleID %in% sample_of_interest)]
  subset_df = subset_df[, `:=` (Shape = data.table::fcase(SampleOfInterest, "Sample of Interest", is.infinite(get(metric)), "Infinite", default = "General Sample"))]

  #subset_df = subset_df[, `:=`(ReadsCoveredByCladeBinned = data.table::fcase(ReadsCoveredByClade < 50, "< 50", ReadsCoveredByClade < 1e3, "< 1000", ReadsCoveredByClade < 1e6, "< 1,000,000", default = ">= 1,000,000")) ]

  set.seed(seed)

  # d = subset_df[SampleID == sample_of_interest,]

  pos_jitter = ggplot2::position_jitter(width = 0, height = 0.5, seed = seed)
#browser()
  plot = subset_df %>%
    ggplot2::ggplot(ggplot2::aes_string(x = metric, y= "ScientificName", SampleID = "SampleID", ReadsCoveredByClade = "ReadsCoveredByClade", RPM = "RPM", Metric = metric, ZscoreRobust = "ZscoreRobust", customdata = "SampleID", ...), alpha = 0.5) +
    ggplot2::geom_jitter(position = pos_jitter, ggplot2::aes(shape = Shape, color =  ReadsCoveredByClade > 50, size = SampleOfInterest)) +
    #ggplot2::geom_jitter(position = pos_jitter, ggplot2::aes(shape = Infinite, color =  ReadsCoveredByClade > 50), stroke = 2) +
    ggplot2::scale_x_continuous(oob = scales::oob_squish_infinite, trans="log10") +
    ggplot2::facet_wrap(~ScientificName, ncol = 1, scales = "free_y") +
    #ggplot2::theme(strip.text.x = ggplot2::element_blank()) +
    ggplot2::scale_shape_manual(values = c("General Sample" = 21, "Infinite" = 24, "Sample of Interest" = 17)) +
    #ggplot2::scale_alpha_manual(values = c(FALSE=0.7, TRUE=1)) +
    ggplot2::scale_size_manual(values = pointsize) +
    #utilitybeltgg::theme_fivethirtyeight_two() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), axis.title.y  = ggplot2::element_blank()) +
    utilitybeltgg::theme_legend_right() +
    ggplot2::scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue"))

  if(show_sample_vlines){
    plot = plot + ggplot2::geom_vline(data = function(df) {df[SampleID %in% sample_of_interest,]} ,ggplot2::aes_string(xintercept = metric, linetype = "SampleID"))
  }

  return(plot)
}

kraken_report_visualise_signal_focus <- function(kraken_report_df, parent_taxids){

  lapply(
    parent_taxids, kraken_calculate_proportion_of_signal_explained_by_n_strongest_taxids(kraken_report_df = as.data.frame(kraken_report_df), n = 1, parent_taxid)
    )

  return(results)
}


#' Title
#'
#' @param kraken_report_df
#' @param samples_of_interest
#' @param rank
#' @param taxonomy_id_to_highlight
#' @param read_threshold
#'
#' @return
#' @export
#'
kraken_report_visualise_single_sample <- function(kraken_report_df, samples_of_interest, rank = "S",taxonomy_id_to_highlight = NULL, read_threshold = 50){
  assertthat::assert_that("ZscoreRobustLoggable" %in% colnames(kraken_report_df), msg = "Report dataframe must include the column [ZscoreRobustLoggable]")

  taxonomy_name_to_highlight <- kraken_report_df[["ScientificName"]][match(taxonomy_id_to_highlight, kraken_report_df[["TaxonomyID"]])]

  #browser()
  gg = kraken_report_df[SampleID %in% samples_of_interest & ReadsCoveredByClade > 0 & Rank == rank,] %>%
    ggplot2::ggplot(
      ggplot2::aes(
        color = ReadsCoveredByClade > 50,
        ReadsCoveredByClade = ReadsCoveredByClade,
        x=ZscoreRobustLoggable,
        ZscoreRobust = ZscoreRobust,
        y=RPM,
        ScientificName = ScientificName,
        TaxonomyID = TaxonomyID,
        #ReadsCoveredByClade = ReadsCoveredByClade,
        shape = is.infinite(ZscoreRobust)
        #color = TaxonomyID %in% taxonomy_id_to_highlight,
        )) +
    ggplot2::scale_shape_manual(values = c("FALSE" = 21,"TRUE" = 24)) +
    ggplot2::scale_x_continuous(trans = "log10", oob = scales::oob_squish_infinite) +
    ggplot2::scale_y_continuous(trans = "log10", oob = scales::oob_squish_infinite) +
    ggplot2::ggtitle(paste0(samples_of_interest, collapse = ",")) +
    ggplot2::guides(shape=ggplot2::guide_legend("Infinite")) +
    ggplot2::geom_point() +
    #ggthemes::theme_fivethirtyeight() +
    #ggthemes::theme_calc() +
    utilitybeltgg::theme_common_adjustments() +
    #utilitybeltgg::theme_fivethirtyeight_two() +
    ggplot2::scale_color_brewer(palette = "Set2")

  if (is.null(taxonomy_id_to_highlight)){
    #gg = gg + ggplot2::guides(color="none")
  }
  else
  {
    gg = gg + ggrepel::geom_text_repel(data = function(df) {df[(TaxonomyID %in% taxonomy_id_to_highlight),]}, ggplot2::aes(label=ScientificName), color="black", min.segment.length = 0, nudge_y=3)
    #gg = gg + ggplot2::guides(color = ggplot2::guide_legend(taxonomy_name_to_highlight))
  }

  return(gg)

}

#' Interactive Single Sample Visualisation
#'
#'
#'
#' @param ... arguments passed to \strong{kraken_report_visualise_single_sample}
#'
#' @return plotly graph
#' @export
#'
kraken_report_visualise_single_sample_interactive <- function(tooltip = c("ScientificName", "TaxonomyID", "ReadsCoveredByClade", "ZscoreRobust", "RPM"), ...){
  plotly::ggplotly(kraken_report_visualise_single_sample(...), tooltip = tooltip)
}

kraken_report_visualise_two_group_comparision <- function(kraken_report_df, taxonomy1 = 2, taxonomy2 = 10239, sample_of_interest = NA){

  subset = kraken_report_df[(TaxonomyID %in% c(taxonomy1, taxonomy2)),]
  taxname1 = subset[["ScientificName"]][match(taxonomy1, subset[["TaxonomyID"]])]
  taxname2 = subset[["ScientificName"]][match(taxonomy2, subset[["TaxonomyID"]])]

  gg = subset %>%
    tidyr::pivot_wider(id_cols = SampleID, names_from = ScientificName, values_from = ZscoreRobustLoggable) %>%
    ggplot2::ggplot(ggplot2::aes_string(x=taxname1, y=taxname2, SampleID = "SampleID")) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "log10", oob = scales::oob_squish_infinite)+
    ggplot2::scale_y_continuous(trans = "log10", oob = scales::oob_squish_infinite) +
    utilitybeltgg::theme_fivethirtyeight_two()

  if(!is.na(sample_of_interest)){
    gg = gg + ggplot2::geom_point(data = function(df) {dplyr::filter(df, SampleID %in% sample_of_interest)} ,ggplot2::aes(color = "Sample of Interest"))
  }

  return(gg)
}


kraken_report_tsne <- function(kraken_report_df){
  #assertthat::assert_that("ZscoreRobust" %in% colnames(kraken_report_df), msg = "Report dataframe must include the column [ZscoreRobustLoggable]")
  assertthat::assert_that("ZscoreRobust" %in% colnames(kraken_report_df), msg = "Report dataframe must include the column [ZscoreRobustLoggable]")
  assertthat::assert_that(nchar(system.file(package = "Rtsne")) > 0, msg = "Tsne plots require that the package `Rtsne` is installed. Run install.packages('Rtsne')")

  message("Using ZscoreRobust column in TSNE dimensionality reduction")

  #Filter for rank of interest
  tsne_input_dt = kraken_report_df[Rank=="S", .(SampleID, ZscoreRobust, ScientificName)]

  #Filter zerovariance groups
  tsne_input_dt[, `:=`(novariance = data.table::uniqueN(RPM) == 1), by = ScientificName]
  tsne_input_dt = tsne_input_dt[novariance == FALSE]

  # Remove samples with infinite values ((alternative - maybe turn those inf values into the max Zscore observed in DF))
  tsne_input_noinfinites_dt = tsne_input_dt %>%
    #dtplyr::lazy_dt(immutable = FALSE) %>%
    dplyr::group_by(ScientificName) %>%
    dplyr::filter(!any(is.infinite(ZscoreRobust))) %>%
    data.table::as.data.table()

  # convert to wide form
  tsne_input_wide_dt = data.table::dcast(tsne_input_dt, SampleID ~ ScientificName, value.var = "ZscoreRobust")


  # Convert to matrix
  tsne_input_wide_mx = as.matrix(tsne_input_wide_dt[,!"SampleID"])

  # Add SampleID as rownmaes
  rownames(tsne_input_wide_mx) <- tsne_input_wide_dt$SampleID
  #return(tsne_input_dt)
}
