test_that("filepath_get_sequence_type", {
  expect_equal(filepath_get_sequence_type("folder/subfolder/file.fastq.gz"), "fastq")
  expect_equal(filepath_get_sequence_type("file.fastq.gz"), "fastq")
  expect_equal(filepath_get_sequence_type("file.fastq"), "fastq")
  expect_equal(filepath_get_sequence_type("file.fq.gz"), "fastq")
  expect_equal(filepath_get_sequence_type("file.fq"), "fastq")

  expect_equal(filepath_get_sequence_type("folder/subfolder/file.fasta.gz"), "fasta")
  expect_equal(filepath_get_sequence_type("file.fasta.gz"), "fasta")
  expect_equal(filepath_get_sequence_type("file.fasta"), "fasta")
  expect_equal(filepath_get_sequence_type("file.fa.gz"), "fasta")
  expect_equal(filepath_get_sequence_type("file.fa"), "fasta")
  expect_equal(filepath_get_sequence_type("file.fna.gz"), "fasta")
  expect_equal(filepath_get_sequence_type("file.fna"), "fasta")

  expect_equal(filepath_get_sequence_type("asdasdasd"), "other")
  expect_equal(filepath_get_sequence_type("folder/file.alfalfa"), "other")


  expect_equal(
    filepath_get_sequence_type(c("file.fna", "file.fastq.gz", "file.fastq", "file.alfalfa")),
    c("fasta", "fastq", "fastq", "other")
    )
})
