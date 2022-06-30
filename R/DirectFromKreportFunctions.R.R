#testkreport_path = system.file(package = "utilitybeltkraken", "DRR002014_1.100sequences.minikraken8GB.kreport")

#' Get Child Taxids
#'
#' Scan a kraken report to find what taxids kraken considered to be the chldren of others
#'
#' @param path_to_kraken_report path to kraken report file.
#' @param taxid taxonomy ID of species you want to find  the relatives of (int)
#' @param inclusive should we include the user-specified taxid in the result (flag)
#'
#' @return taxids that the kraken report says are children of the user-supplied taxid (int)
#' @export
#'
kraken_report_get_child_taxids <- function(path_to_kraken_report, taxid, inclusive = FALSE){

  assertthat::assert_that(assertthat::is.flag(inclusive))

  kraken_report_df = kraken_report_parse(path_to_kraken_report, verbose = FALSE)
  row_containing_taxid = match(taxid, kraken_report_df[["TaxonomyID"]])

  #If taxid isnt in kraken report
  if(is.na(row_containing_taxid)){
    stop(paste0("Taxonomy ID [",taxid , "] was not present in Kraken Report.... Exiting"))
  }


  level_of_specified_taxid = kraken_report_df[[row_containing_taxid, "Level"]]
  rows_with_levels_equivalent_or_lower_than_taxid <- which(kraken_report_df[["Level"]] <= level_of_specified_taxid)

  first_row_since_taxid_thats_not_its_child = rows_with_levels_equivalent_or_lower_than_taxid[rows_with_levels_equivalent_or_lower_than_taxid > row_containing_taxid][1]

  if(is.na(first_row_since_taxid_thats_not_its_child)){
    first_row_since_taxid_thats_not_its_child = nrow(kraken_report_df)+1
  }

  child_taxid_rows = row_containing_taxid:(first_row_since_taxid_thats_not_its_child-1)
  child_taxids = as.numeric(unlist(kraken_report_df[child_taxid_rows, "TaxonomyID"]))

  # IF we only find the taxid specified by the user, return just that taxid if inclusive = TRUE, If inclusive = FALSE, throw error
  if(length(child_taxids) == 1){
    if(inclusive)
      return(child_taxids)
    else
      stop("Couldn't find any descendants of taxid [", taxid, "]")
  }

  # If we find any child taxids we return them (and optionally also the parent taxid if `inclusive` flag is set)
  if(inclusive)
    return(child_taxids)
  else
    return(child_taxids[-1])
}


