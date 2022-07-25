#' Sequence types of filepaths
#'
#' @param filepaths paths of sequence file (fastq/fasta/other)
#'
#' @return  'fastq', 'fasta' or 'Other' depending on the filetype of the paths supplied (character)
#' @export
#'
#' @examples
#' filepath_get_sequence_type("file.fastq.gz")
filepath_get_sequence_type <- function(filepaths){
  input_seq_class = dplyr::case_when(
    grepl(filepaths, pattern = "\\.(fastq|fq)(.gz)?$") ~ "fastq",
    grepl(filepaths, pattern = "\\.(fasta|fa|fna)(.gz)?$") ~ "fasta",
    TRUE ~ "other"
  )
  return(input_seq_class)
}
