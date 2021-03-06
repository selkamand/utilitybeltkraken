% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KrakenParsing.R
\name{kraken_info_taxids_of_interest}
\alias{kraken_info_taxids_of_interest}
\title{Taxids of Interest}
\usage{
kraken_info_taxids_of_interest(purpose = "info")
}
\arguments{
\item{purpose}{one of ("info","general", "Detecting Infection from Human Samples", "Detecting Viruses Potentially Related to Paediatric Cancer")}
}
\value{
Named numeric vector with some commonly used taxids (returns invisible vector)
}
\description{
Taxids of Interest
}
\examples{
kraken_info_taxids_of_interest()
}
