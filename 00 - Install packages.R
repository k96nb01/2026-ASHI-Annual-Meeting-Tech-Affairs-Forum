# Install the packages used in this session.
# On Posit Cloud these are already installed in the shared project, so you only
# need this if you are running the code on your own computer. Run it once; it
# takes a few minutes.
install.packages("pak")
pak::pkg_install(c(
  "tidyverse",     # data tables: select(), filter(), mutate() and friends
  "immunogenetr",  # HLA typing: GL strings, mismatches, match grades
  "HLAtools",      # the IPD-IMGT/HLA Database: alignments and allele names
  "httr2",         # talking to web services, like the matching API
  "shiny",         # web apps
  "bslib",         # page layouts for Shiny apps
  "rmarkdown"      # rendering the .qmd files into documents
))
