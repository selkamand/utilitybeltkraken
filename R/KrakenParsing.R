#' kraken_report_parse
#'
#' If you want to import multiple reports use \strong{kraken_report_parses}
#'
#' @param kraken2report path to kraken directory (string)
#' @param sample_id (string)
#'
#' @return dataframe describing the kraken report
#' @export
#'
kraken_report_parse2 <- function(kraken2report, sample_id=NA){
  # Read Data
  kreport_headings <- c("PercentReadsCoveredByCladeLowResolution","ReadsCoveredByClade", "ReadsDirectlyAssigned", "Rank", "TaxonomyID", "ScientificName")
  kraken_df <- read.csv(file = kraken2report, header = FALSE, sep = "\t", col.names = kreport_headings)

  #Use number of indents in scientific name to extrapolate depth
  kraken_df$Level <- stringr::str_count(kraken_df$ScientificName, pattern = "  ")
  kraken_df$LevelDrops <- c(kraken_df$Level[-length(kraken_df$Level)]-kraken_df$Level[-1] < 0, FALSE)
  # Strip tabs from scientific name column
  kraken_df$ScientificName <- sub(x = kraken_df$ScientificName, pattern = "^ +", replacement = "")

  # Simplify Rank code - We don't care about differentiating between S, S1, or S2 ranks - theyre all species level
  kraken_df$RankSimple <- gsub(x = kraken_df$Rank, pattern = "[0-9]", replacement = "")

  # Recalculate Percent Reads Covered as normalised score. Due to no. of sci figures supported by kraken2 - clades with 100,000 reads classifed will still come out as '0%'. We fix this below
  total_reads = sum(kraken_df$ReadsDirectlyAssigned)
  kraken_df$PercentReadsCoveredByClade <- kraken_df$ReadsCoveredByClade*100/total_reads

  # Calculate RPM (reads covered by clade per million total reads
  kraken_df$RPM <- kraken_df$ReadsCoveredByClade*1e6/total_reads

  #Add SampleID column derived from report filename
  if(is.na(sample_id))
    kraken_df$SampleID = sub(x=basename(kraken2report), pattern = '\\..*', replacement = "")
  else
    kraken_df$SampleID = sample_id
  return(kraken_df)
}

#' Parse Kraken Reports
#'
#' @param kraken2directory path to a directory filled with ONLY kraken2 reports
#'
#' @return a dataframe describing all samples kraken reports
#' @export
kraken_reports_parse2 <- function(kraken2directory){
  krakenreportpaths <- dir(path = kraken2directory, full.names = TRUE, recursive = FALSE, include.dirs = FALSE)
  #browser()
  kraken_repots_df <- purrr::map_dfr(
    krakenreportpaths,
    ~ kraken_report_parse(
        kraken2report = .x
      )
    )
  # Maybe add code to make sure all samples have the same number of rows
}


#' Parse Kraken Reports
#'
#' @param kraken2directory path to a directory filled with ONLY kraken2 reports
#'
#' @return a dataframe describing all samples kraken reports
#' @export
kraken_reports_parse <- function(kraken2directory){

  assertthat::assert_that(dir.exists(kraken2directory), msg = paste0("Could not find directory: ", kraken2directory))
  kreport_headings <- c("PercentReadsCoveredByCladeLowResolution","ReadsCoveredByClade", "ReadsDirectlyAssigned", "Rank", "TaxonomyID", "ScientificName")
  krakenreportpaths <- dir(path = kraken2directory, full.names = TRUE, recursive = FALSE, include.dirs = FALSE)

  # Read Files into one huge df
  kraken_reports_df <- data.table::rbindlist(
    lapply(krakenreportpaths, data.table::fread, col.names = kreport_headings, sep="\t", strip.white = FALSE) %>%
      magrittr::set_names(basename(krakenreportpaths)),
    idcol = "Filename")

  #browser()
  #Use number of indents in scientific name to extrapolate depth
  kraken_reports_df[, `:=` (Level = stringr::str_count(ScientificName, pattern = "  "))]
  #kraken_reports_df$Level <- stringr::str_count(kraken_df$ScientificName, pattern = "  ")

  #kraken_reports_df[, `:=` (LevelDrops = c(kraken_df$Level[-length(kraken_df$Level)]-kraken_df$Level[-1] < 0, FALSE))]
  #kraken_df$LevelDrops <- c(kraken_df$Level[-length(kraken_df$Level)]-kraken_df$Level[-1] < 0, FALSE)

  # Strip tabs from scientific name column
  kraken_reports_df[, `:=` (ScientificName = stringr::str_replace(string = ScientificName, pattern = "^ +", replacement = ""))]
  #kraken_df$ScientificName <- sub(x = kraken_df$ScientificName, pattern = "^ +", replacement = "")

  # Simplify Rank code - In case we don't care about differentiating between S, S1, or S2 ranks - theyre all species level -- most of the time original rankings should be fine
  kraken_reports_df[, `:=` (RankSimple = stringr::str_replace_all(string =  Rank, pattern = "[0-9]", replacement = ""))]
  #kraken_df$RankSimple <- gsub(x = kraken_df$Rank, pattern = "[0-9]", replacement = "")

  # Recalculate Percent Reads Covered as normalised score. Due to no. of sci figures supported by kraken2 - clades with 100,000 reads classifed will still come out as '0%'. We fix this below
  total_reads = sum(kraken_reports_df$ReadsDirectlyAssigned)
  kraken_reports_df[, `:=`(PercentReadsCoveredByClade = ReadsCoveredByClade*100/total_reads) ]

  # total_reads = sum(kraken_df$ReadsDirectlyAssigned)
  # kraken_df$PercentReadsCoveredByClade <- kraken_df$ReadsCoveredByClade*100/total_reads

  # Calculate RPM (reads covered by clade per million total reads
  kraken_reports_df[, `:=`(RPM = ReadsCoveredByClade*1e6/total_reads)]
  #kraken_df$RPM <- kraken_df$ReadsCoveredByClade*1e6/total_reads

  #Add SampleID column derived from report filename
  kraken_reports_df[, `:=`(SampleID = stringr::str_replace(string = Filename, pattern = '\\..*', replacement = ""))]
  #kraken_df$SampleID = sub(x=basename(kraken2report), pattern = '\\..*', replacement = "")
  return(kraken_reports_df)

  # Maybe add code to make sure all samples have the same number of rows
}

#' Identify Taxid Descendancy Status
#'
#' Identify taxids in your \strong{kraken_report_df} that are either equal to the user specified taxonomy ID or one of its descendants/children
#' This is seful for telling what species are bacterial / viral / belong to a particular genus etc.
#' \strong{WARNING: } this function relies on inherent properties of kraken reports - make sure you run this ONLY on unsorted dataframes produced by \link{kraken_reports_parse} or link{kraken_report_parse}
#'
#'
#' @param kraken_report_df an UNSORTED kraken report dataframe produced by \link{kraken_reports_parse} or link{kraken_report_parse}
#' @param taxid ncbi taxonomic id (integer)
#' @param careful should we do a bunch of quality assertions (slight speed overhead) (flag)
#' @param columname name of the column to add to the dataframe
#'
#'
#' @return This function returns the same dataframe from \strong{kraken_report_df} produced by \link{kraken_reports_parse} or link{kraken_report_parse}  but with  a new columns describing inclusive descendancey status to a particular taxid
#' @export
#'
kraken_report_add_descendancy_status <- function(kraken_report_df, taxid, columname=NA, careful=TRUE, verbose = T){
  assertthat::assert_that(assertthat::is.flag(careful))

  if(careful){
    assertthat::assert_that(assertthat::is.number(taxid))
    assertthat::assert_that(is.data.frame(kraken_report_df))
    assertthat::assert_that(assertthat::has_name(kraken_report_df, c("TaxonomyID", "Level", "SampleID", "ReadsCoveredByClade", "ScientificName")))
    assertthat::assert_that(taxid %in% kraken_report_df[["TaxonomyID"]])
    assertthat::assert_that(dplyr::n_distinct(table(kraken_report_df[["TaxonomyID"]])) == 1, msg = "Some Taxonomy Ids appear more than others. Its likely you haven't used --report-zero-counts. I would advise you add this flag as it makes it very easy to catch problems arising from different databases. All the functions in this package will still work fine - just run this code again with the option careful=FALSE")
    assertthat::assert_that(!is.unsorted(kraken_report_df[["SampleID"]]))

  }

  taxid_index = which(taxid == kraken_report_df$TaxonomyID)
  taxid_level = unique(kraken_report_df$Level[taxid_index])
  if(length(taxid_level) > 1) stop("Found inconsistencies in krakenreport formating -- same taxid is at different levels for different samples. Are you sure you ran all samples against the same database")
  terminating_indexes = which(kraken_report_df$Level <= taxid_level)
  terminating_indexes = terminating_indexes[!terminating_indexes %in% taxid_index]-1
  lastSampleEntry=which(c(kraken_report_df$SampleID[-length(kraken_report_df$SampleID)] != kraken_report_df$SampleID[-1], TRUE))
  terminating_indexes = sort.int(c(terminating_indexes, lastSampleEntry), decreasing = FALSE)

  matching_terminating_indexes = purrr::map_int(taxid_index, ~as.integer(terminating_indexes[terminating_indexes >= .x])[1])

  indexes <- as.vector(sapply(seq_along(taxid_index), FUN = function(i){taxid_index[i]:matching_terminating_indexes[i]}))

  if(is.na(columname))
    columname = paste0("DescendancyFrom", gsub(x=kraken_report_df$ScientificName[match(taxid, kraken_report_df$TaxonomyID)], pattern = " ", replacement = "_"))

  if (verbose)
    message("Adding Column: ", columname)

  kraken_report_df[columname] = seq_along(kraken_report_df$TaxonomyID) %in% indexes
  #browser()
  return(kraken_report_df)
}



#' Quantifying Signal Spread Across Taxid
#'
#' This function tells you the proportion of reads belonging to a taxid (or its descendents) whose classification within the larger group is driven by classification into the most commonly classified species/genus/family (user can choose).
#' Say we notice a high number of bacterial reads and we want to know if this phenomina is driven largely by a single species or a more diffuse binning of reads into lots of different bacterial species.
#' The low spread binning (most in single species) would increase our confidence that the sample actually has the in it. If binning is more diffuse its more likely to have some other explanation.
#'
#' @param parent_taxid the taxid representing a larger group within which we want to know whether binning of reads is diffuse or focused (integer)
#' @param n  use the n most binned species/genus'/etc to calculate 'proportion of reads in parent_taxid lineage explained by a focused group' (integer)
#' @param rank one of: (S, G, F, O, C, P, K). S=speies, G=genus, etc. Determines which rank in calculation of binning 'diffuseness'. You will almost always want to use species (S), maybe genus.
#' @inheritParams kraken_report_add_zscore
#'
#' @return Dataframe describing the percentage of reads in parent clade that are explained by the n most frequently classified taxids
#' @export
kraken_calculate_proportion_of_signal_explained_by_n_strongest_taxids <- function(kraken_report_df,parent_taxid, n=1, rank = "S"){
  assertthat::assert_that(assertthat::is.number(parent_taxid))
  assertthat::assert_that(is.data.frame(kraken_report_df))
  assertthat::assert_that(assertthat::has_name(kraken_report_df, c("TaxonomyID", "Level", "SampleID", "ReadsCoveredByClade", "ScientificName")))
  assertthat::assert_that(parent_taxid %in% kraken_report_df[["TaxonomyID"]])
  assertthat::assert_that(dplyr::n_distinct(table(kraken_report_df[["TaxonomyID"]])) == 1, msg = "Some Taxonomy Ids appear more than others. Its likely you haven't used --report-zero-counts. I would advise you add this flag as it makes it very easy to catch problems arising from different databases. All the functions in this package will still work fine - just run this code again with the option careful=FALSE")

  taxid_name = kraken_report_df[["ScientificName"]][match(parent_taxid, kraken_report_df[["TaxonomyID"]])]
  #reads_covering_parent_clade = kraken_report_df[["ReadsCoveredByClade"]][taxid %in% ]
  #browser()
  kraken_report_with_descendancy_status = kraken_report_add_descendancy_status(kraken_report_df = kraken_report_df, taxid = parent_taxid, columname = "InParentClade", verbose=F)
  reads_covered_by_clade <- kraken_report_with_descendancy_status %>%
    dplyr::filter(TaxonomyID == parent_taxid) %>%
    dplyr::distinct(SampleID, ReadsCoveredByClade,.keep_all = TRUE) %>%
    dplyr::select(SampleID, ReadsCoveredByParentClade = ReadsCoveredByClade)

  #browser()

  kraken_report_with_descendancy_status %>%
    dplyr::filter(InParentClade==TRUE,Rank == rank) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(ContributionRank = bit64::rank(ties.method="first", -ReadsCoveredByClade)) %>%
    dplyr::summarise(
      TopRankClade = paste0(ScientificName[ContributionRank<=n], collapse = ","),
      TopRankCladeTaxid = paste0(TaxonomyID[ContributionRank<=n], collapse = ","),
      ReadsCoveredByTopRankClade = sum(ReadsCoveredByClade[ContributionRank<=n]),
      #ReadsCoveredByParentClade = reads_covered_by_clade$ReadsCoveredByClade[match(SampleID, reads_covered_by_clade$SampleID)],
      #PercentageOfReadsInParentCladeCoveredByTopRankClade = ReadsCoveredByTopRankClade*100/ReadsCoveredByParentClade
      ) %>%
    dplyr::left_join(reads_covered_by_clade, by="SampleID") %>%
    dplyr::mutate(
      ParentCladeQueried = parent_taxid,
      Nused = n,
      PercentageOfReadsInParentCladeCoveredByTopRankClade = ReadsCoveredByTopRankClade*100/ReadsCoveredByParentClade,
      RPM = ReadsCoveredByTopRankClade*1000000/ReadsCoveredByParentClade
      )


}

#' Taxids of Interest
#'
#'
#' @param purpose one of ("general" or "humaninfectiondetection", or "pedcanviral")
#'
#' @return Named numeric vector with some commonly used taxids (returns invisible vector)
#' @export
#'
#' @examples
#' kraken_info_taxids_of_interest()
kraken_info_taxids_of_interest <- function(purpose = "info"){
  purpose_possible_values = c("info","general", "Detecting Infection from Human Samples", "Detecting Viruses Potentially Related to Paediatric Cancer")
  assertthat::assert_that(purpose %in% purpose_possible_values, msg = paste0("[purpose] must be one of [", paste0(purpose_possible_values, collapse = ", "), "]"))

  if(purpose == "info"){
    message("Please run one of the following based on your goal: \n\n",
            'taxids <- kraken_info_taxids_of_interest("general")\n\n',
            'taxids <- kraken_info_taxids_of_interest("Detecting Infection from Human Samples")\n\n',
            'taxids <- kraken_info_taxids_of_interest("Detecting Viruses Potentially Related to Paediatric Cancer")\n\n'
            )
  }
  else if(purpose == "general"){
    message("Listing General Purpose Taxids")
    general_vector <- c(
      "unclassified" = 0,
      "root (a.k.a classified)" = 1,
      "Bacteria" = 2,
      "Viruses" = 10239,
      "Fungi" = 4751,
      "Archaea" = 2157,
      "SAR" = 2698737,
      "Homo sapiens" = 9606,
      "Homo" = 9605,
      "Hominidae" = 9604,
      "Chordata" = 7711,
      "Other Sequences" = 28384,
      "Escherichia virus phiX174" = 10847
      #"Eukaryota" = 2759,
      )

    message(paste0("[",general_vector, "]", "\t\t",names(general_vector), collapse = "\n"))
    message("\nNotes:\n> SAR is a subdomain of eukaryotic microbes that includes many eukaryotic microbes including protists")
    message("> Bacteria, Viruses, and Archaea are Domains, while Fungi is a kingdom within the eukarya domain")
    message("> Chordata is a phylum containing humans and many other animals - no microbes")
    message("> If Using PlusPF DB: Other sequences refers to UniVecCore sequences (common artificial contaminants, plasmids, etc)\n  and includes PhiX")
    return(invisible(general_vector))
  }
  else if(purpose == "Detecting Viruses Potentially Related to Paediatric Cancer"){
    pedcanviral = c(
      "Human gammaherpesvirus 4 (EBV)" = 10376,
      "Human betaherpesvirus 6B (HHV-6B)" = 32604,
      "Hepatitis B virus (Hep B)" = 10407,
      "Torque teno virus (TTV)" = 68887,
      "Human papillomavirus (HPV)" = 10566,
      "Cytomegalovirus (CMV)" = 10358,
      "Macaca mulatta polyomavirus 1 (SV40)" = 1891767,
      "Human polyomavirus 1 (BKV)" = 1891762
      )
    message(paste0("[",pedcanviral, "]", "\t\t",names(pedcanviral), collapse = "\n"))
    return(invisible(pedcanviral))
  }
  else if (purpose == "Detecting Infection from Human Samples"){
    humaninfectiondetection = c(
      "unclassified" = 0,
      "Homo" = 9605,
      "Viruses" = 10239,
      "Bacteria"  = 2,
      "Fungi" = 4751,
      "SAR" = 2698737,
      message("\nNotes:\n> SAR is a subdomain of eukaryotic microbes that includes many eukaryotic microbes including protists")
      )
    message(paste0("[",humaninfectiondetection, "]", "\t\t",names(humaninfectiondetection), collapse = "\n"))
    return(invisible(humaninfectiondetection))
  }
  else
    stop("purpose must be one of: ",paste0(purpose_possible_values, collapse=", "))
}

KrakenParsingSingleFileAddLineage <- function(krakenreport){
  kreport_headings <- c("PercentReadsCoveredByCladeLowResolution","ReadsCoveredByClade", "ReadsDirectlyAssigned", "Rank", "TaxonomyID", "ScientificName")
  data.table::fread(col.names = kreport_headings, sep="\t", strip.white = FALSE) %>%
  magrittr::set_names(basename(krakenreport))
}


