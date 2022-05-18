

#' Convert Kraken Report Dataframe to sqlite database
#'
#' @inherit kraken_report_add_descendancy_status
#' @param database_name (If unsure leave defaults) name of sqlite database to create.
#' @param table_name (If unsure leave defaults) name of kraken report dataframe table in sqlite database
#'
#' @return Name of the Database Created (string). Primarily Run for its side effects
#' @export
#'
kraken_report_to_sqlite_db = function(kraken_report_df, database_name = paste0(getwd(),"/kraken_report_database.sqlite"), table_name = "kreports"){

  if(file.exists(database_name)){
    message("Already found a database file with the given name [", tools::file_path_as_absolute(database_name), "]")
   overwrite = askYesNo(msg = paste0("Are you sure you want to overwrite the existing database"))

    if(!overwrite){
      message("Leaving database untouched")
      return()
    }
   else message("Overwriting existing database ... ")
  }

  message("Creating sqlite database [",tools::file_path_as_absolute(database_name),"]")
  kreport_sqlite_db <- DBI::dbConnect(RSQLite::SQLite(), database_name)

  message("Converting kraken_report_df to sqlite database table [",table_name,"]")
  DBI::dbWriteTable(kreport_sqlite_db, name = table_name, value = kraken_report_df, overwrite = TRUE)

  message("Creating SampleID TaxonomyID multicolumn index:")
  DBI::dbSendQuery(conn = kreport_sqlite_db, statement = paste0(
    "CREATE INDEX sample_taxonomy_index
    ON ",table_name," (SampleID, TaxonomyID);
    "
    ))

  message("Creating TaxonomyID SampleID index:")
  DBI::dbSendQuery(conn = kreport_sqlite_db, statement = paste0(
    "CREATE INDEX taxonomy_sample_index
    ON ",table_name," (TaxonomyID, SampleID);
    "))


  message("Creating TaxonomyID to ScientificName Mapping Views:")
  DBI::dbSendQuery(conn = kreport_sqlite_db, statement = paste0(
    "CREATE VIEW species AS
    SELECT DISTINCT(TaxonomyID, ScientificName) FROM ",table_name, ";"))

  message("Creating TaxonomyID, ScientificName Index:")
  DBI::dbSendQuery(conn = kreport_sqlite_db, statement = paste0(
    "CREATE INDEX taxonomy_sample_index
    ON ",species," (TaxonomyID, ScientificName);
    "))

  DBI::dbDisconnect(kreport_sqlite_db)

  message("")
  message("All Finished!

          Connect to the database using the following code snippet:
          DBI::dbConnect(
            RSQLite::SQLite(),
            '",database_name,"'
          )"
          )
  return(database_name)
}


#' kraken_reports_parse_to_sqlite_db
#'
#' Read kraken reports from a directory and convert to sql database
#'
#' @inheritDotParams kraken_report_to_sqlite_db
#' @inheritParams kraken_reports_parse
#' @param ... additional paramaaters to pass to [utilitybeltkraken::kraken_report_to_sqlite_db()]
#'
#' @return database name
#' @export
#'
kraken_reports_parse_to_sqlite_db <- function(krakenreport_directory, ...){
  message("Parsing kraken reports")
  kraken_df <- kraken_reports_parse(krakenreport_directory)
  kraken_report_add_robust_zscore(kraken_df)

  message("Starting Conversion to sqlite database")
  dbname = kraken_report_to_sqlite_db(kraken_report_df = kraken_df, ...)
  return(dbname)
}


