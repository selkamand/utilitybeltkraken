

#' Convert Kraken Report Dataframe to sqlite database
#'
#' @inherit kraken_report_add_descendancy_status
#' @param database_name (If unsure leave defaults) name of sqlite database to create.
#' @param table_name (If unsure leave defaults) name of kraken report dataframe table in sqlite database
#'
#' @return Name of the Database Created (string). Primarily Run for its side effects
#' @export
#'
kraken_report_to_sqlite_db = function(kraken_report_df, database_name = paste0(getwd(),"/kraken_report_database.sqlite"), table_name = "kreports", overwrite = TRUE){

  if(file.exists(database_name)){
    message("Already found a database file with the given name [", tools::file_path_as_absolute(database_name), "]")
   #overwrite = askYesNo(msg = paste0("Are you sure you want to overwrite the existing database [",database_name,"]"))

    if(!overwrite){
      message("Leaving database untouched")
      return()
    }
   else {
     message("Overwriting existing database ... ")
     file.remove(tools::file_path_as_absolute(database_name))
   }

  }

  message("Creating sqlite database [",database_name,"]")
  kreport_sqlite_db <- DBI::dbConnect(RSQLite::SQLite(), database_name)

  message("Converting kraken_report_df to sqlite database table [",table_name,"]")
  DBI::dbWriteTable(kreport_sqlite_db, name = table_name, value = kraken_report_df, overwrite = TRUE)


  # Indexing ----------------------------------------------------------------
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


  message("Creating TaxonomyID, ScientificName Index:")
  DBI::dbSendQuery(conn = kreport_sqlite_db, statement = paste0(
    "CREATE INDEX taxonomy_scientific_name_index
    ON ",table_name," (TaxonomyID, ScientificName);
    "))

  message("Creating Rank, Zscore, RPM Index")
  DBI::dbSendQuery(conn = kreport_sqlite_db, statement = paste0(
    "CREATE INDEX rank_zscore_rpm
    ON ",table_name," (Rank, ZscoreRobust, RPM);
    "))

  # message("Creating Rank Index")
  # DBI::dbSendQuery(conn = kreport_sqlite_db, statement = paste0(
  #   "CREATE INDEX rank
  #   ON ",table_name," (Rank);
  #   "))

  # Create Views ----------------------------------------------------------------
  message("Creating TaxonomyID to ScientificName Mapping Views:")
  DBI::dbSendQuery(conn = kreport_sqlite_db, statement = paste0(
    "CREATE VIEW species AS
    SELECT DISTINCT(TaxonomyID, ScientificName) FROM ",table_name, ";"))

  DBI::dbDisconnect(kreport_sqlite_db)

  message("All Finished!

          Database files can be found at: ",
          message(tools::file_path_as_absolute(database_name))
          ,"\n

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


#' Connect to Kraken Database
#'
#' @param path_to_sqlite_db path to sqlite db file produced by [kraken_reports_parse_to_sqlite_db()]
#'
#' @return dbplyr object
#' @export
#'
kraken_database_connect <- function(path_to_sqlite_db){
  assertthat::assert_that(file.exists(path_to_sqlite_db), msg = paste0("Cant find database file at [", path_to_sqlite_db,"]"))
  connection <- DBI::dbConnect(RSQLite::SQLite(), path_to_sqlite_db)
  kreports_db <- dplyr::tbl(connection, "kreports")
  return(kreports_db)
}

#' Disconnect from Kraken Database
#'
#' @param kreports_db object produced by [kraken_database_connect()]
#'
#' @return run for its side effects
#' @export
#'
#' @examples
#' \dontrun{
#' # open connection to database
#' kreports_db = kraken_database_connect("path_to_database")
#'
#' # do required visualisations/analyses
#'
#' # close connection when done
#' kraken_database_close_connection(kreports_db)
#' }
kraken_database_close_connection <- function(kreports_db){
  connnection=kreports_db[["src"]][["con"]]
  DBI::dbDisconnect(connnection)
}


