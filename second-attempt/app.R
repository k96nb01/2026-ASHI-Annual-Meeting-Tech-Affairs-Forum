# =============================================================================
# HLA-DPB1 epitope viewer
#
# A small Shiny app that shows the three major antibody epitopes of any
# HLA-DPB1 allele (Morris et al., PMID 39277973):
#
#   Epitope   Amino acid positions   Common sequences
#   55        55, 56, 57             AAE, DEE, DED, EAE
#   84        84, 85, 86, 87         DEAV, GGPM, VGPM
#   96/170    96 and 170             K/I, R/T
#
# The protein sequences come from the IPD-IMGT/HLA Database, downloaded with
# the HLAtools package. Allele names are handled with the immunogenetr package.
#
# To run it: open this file in RStudio and click "Run App", or run
#   shiny::runApp("app")
# =============================================================================


# ---- Packages ---------------------------------------------------------------
# Each library() call loads a package: a bundle of ready-made functions.
library(shiny)        # builds the web app
library(bslib)        # gives the app a clean, modern layout
library(dplyr)        # tools for working with data tables
library(stringr)      # tools for working with text (used only for residues)
library(HLAtools)     # downloads alignments from the IPD-IMGT/HLA Database
library(immunogenetr) # understands HLA allele names


# ---- Settings ---------------------------------------------------------------
# The IPD-IMGT/HLA Database releases a new version every few months, and new
# versions add alleles. We "pin" one version here so that the same allele
# always gives the same answer, no matter when the app is run. To use a newer
# release, change this one line.
IPD_IMGT_HLA_VERSION <- "3.65.0"

# The allele shown when the app first opens, so the demo starts with a result.
DEFAULT_ALLELE <- "DPB1*04:01:01:01"

# The common sequences for each epitope (Morris et al., PMID 39277973).
COMMON_55     <- c("AAE", "DEE", "DED", "EAE")
COMMON_84     <- c("DEAV", "GGPM", "VGPM")
COMMON_96_170 <- c("K/I", "R/T")

# The amino acid positions the epitopes are built from.
EPITOPE_POSITIONS <- c("55", "56", "57", "84", "85", "86", "87", "96", "170")

# Symbols in the alignment that are NOT amino acids:
#   "*"  the residue is not known (that part was never sequenced)
#   "."  there is no residue here (a gap, or the protein has already ended)
#   "X"  a stop codon (the protein ends here)
# An epitope that contains any of these is reported as "not known".
NOT_AMINO_ACIDS <- c("*", ".", "X")

# Plain-language labels used when one of those symbols is shown in the table,
# so the audience never mistakes them for amino acids.
NOT_AMINO_ACID_LABELS <- c(
  "*" = "not sequenced",
  "." = "no residue",
  "X" = "stop codon"
)


# ---- Load the alignment (once, when the app starts) -------------------------
# alignmentFull() downloads the DPB1 protein alignment from the IPD-IMGT/HLA
# Database. It returns a list; the table we want is inside it at
# [["prot"]][["DPB1"]]. The table has one row per allele and one column per
# amino acid position ("55", "84", "170", ...).
#
# This code sits outside the server function, so it runs only once when the
# app starts, not every time someone picks an allele.
dpb1_alignment <- alignmentFull(
  loci      = "DPB1",
  alignType = "prot",
  version   = IPD_IMGT_HLA_VERSION
)[["prot"]][["DPB1"]]

# Keep only the columns we need: the full allele name and the nine positions.
dpb1 <- dpb1_alignment |>
  select(allele_name, all_of(EPITOPE_POSITIONS))

# Every full allele name, in alignment order, for the drop-down list.
all_alleles <- dpb1$allele_name


# ---- Helper: clean up a typed allele name -----------------------------------
# People type allele names in many ways: " 04:01 ", "DPB1*04:01",
# "HLA-DPB1*04:01". This turns all of them into one standard form,
# "DPB1*04:01", using immunogenetr functions only.
# It returns NA if the text is not a valid allele name, or is for another gene.
clean_allele_name <- function(typed) {
  # HLA_validate() removes stray spaces and returns NA for text that is not
  # an HLA name at all (for example "abc" or "").
  valid <- HLA_validate(typed)
  if (is.na(valid)) return(NA_character_)

  # HLA_prefix_remove() strips "HLA-" and the gene, leaving just the fields
  # ("04:01"). HLA_prefix_add() then puts back "DPB1*", the form the
  # alignment uses.
  fields_only <- HLA_prefix_remove(valid)
  cleaned     <- HLA_prefix_add(fields_only, "DPB1*")

  # If the user typed a gene, it must be DPB1: we compare the name with the
  # gene kept (keep_locus = TRUE) to the cleaned name. A name typed without
  # any gene ("04:01") passes this check too.
  with_gene <- HLA_prefix_remove(valid, keep_locus = TRUE)
  if (with_gene != fields_only && with_gene != cleaned) return(NA_character_)

  cleaned
}


# ---- Helper: how many fields does a name have? ------------------------------
# "DPB1*04:01" has 2 fields, "DPB1*04:01:01" has 3, and so on. Instead of
# counting colons ourselves, we ask HLA_truncate(): the smallest number of
# fields that leaves the name unchanged is the number of fields it has.
count_fields <- function(allele) {
  for (n in 1:4) {
    if (HLA_truncate(allele, fields = n) == allele) return(n)
  }
  4
}


# ---- Helper: is this a null allele? ------------------------------------------
# A null allele carries the expression suffix "N" (for example
# DPB1*61:01N): it is not expressed on the cell surface.
# GLstring_drop_non_expressed() removes alleles that are not expressed; given
# a single null allele, nothing is left, so it returns NA. We tell it to treat
# only "N" as non-expressed, so alleles with other suffixes (such as "Q",
# "questionable" expression) are not called null.
is_null_allele <- function(allele) {
  is.na(GLstring_drop_non_expressed(allele, suffixes = "N"))
}


# ---- Helper: build one epitope and decide its status ------------------------
# Takes the residues for one epitope (for example "D", "E", "A", "V") and the
# list of common sequences. Returns the epitope text and its status.
make_epitope <- function(residues, common, separator = "") {
  # If any residue is missing or is not an amino acid, we cannot say what the
  # epitope is.
  if (any(is.na(residues)) || any(residues %in% NOT_AMINO_ACIDS)) {
    return(list(sequence = "not known", status = "not known"))
  }
  # Otherwise join the residues together: "DEAV", or "K/I" for 96/170.
  sequence <- str_c(residues, collapse = separator)
  status   <- if (sequence %in% common) "common" else "uncommon"
  list(sequence = sequence, status = status)
}


# ---- The epitope function ----------------------------------------------------
# This is the heart of the app. Give it an allele name, and it returns the
# epitopes. It is a plain R function, so you can try it in the console:
#
#   dpb1_epitopes("04:01")
#   dpb1_epitopes("DPB1*02:01:02:01")
#
# It returns a list with:
#   found      TRUE if a matching allele was found
#   message    a plain-language message (why nothing was found, if so)
#   typed      what was typed
#   shown      the full allele name whose epitopes are shown
#   label      a short heading, e.g. "DPB1*04:01 (shown: DPB1*04:01:01:01)"
#   is_null    TRUE if the shown allele is a null allele
#   epitopes   a table: epitope, sequence, status
#   residues   a table: the residue at each of the nine positions
dpb1_epitopes <- function(allele) {

  # Step 1: clean up what was typed.
  query <- clean_allele_name(allele)
  if (is.na(query)) {
    return(list(
      found   = FALSE,
      typed   = allele,
      message = paste0("\"", allele, "\" is not a valid HLA-DPB1 allele name.")
    ))
  }

  # Step 2: find matching alleles. We shorten every allele name in the
  # alignment to the same number of fields as the query, then compare whole
  # names. HLA_truncate() keeps expression suffixes, so "DPB1*04:01:01:24N"
  # becomes "DPB1*04:01N", which does NOT equal "DPB1*04:01". That is how a
  # null allele is never returned for a name that did not ask for one.
  n_fields <- count_fields(query)
  matches  <- dpb1 |>
    filter(HLA_truncate(allele_name, fields = n_fields, keep_suffix = TRUE) == query)

  if (nrow(matches) == 0) {
    return(list(
      found   = FALSE,
      typed   = allele,
      message = paste0(
        query, " was not found in IPD-IMGT/HLA Database version ",
        IPD_IMGT_HLA_VERSION, "."
      )
    ))
  }

  # Step 3: when several alleles match (for example "04:01" matches
  # DPB1*04:01:01:01, DPB1*04:01:01:02, ...), use the first one in the
  # alignment, and say which one is shown.
  row   <- matches |> slice(1)
  shown <- row$allele_name
  label <- if (shown == query) shown else paste0(query, " (shown: ", shown, ")")

  # Step 4: pull out the residue at each position as a plain named vector,
  # e.g. r["55"] is the amino acid at position 55.
  r <- unlist(row[EPITOPE_POSITIONS])

  # Step 5: build the three epitopes.
  e55  <- make_epitope(r[c("55", "56", "57")],       COMMON_55)
  e84  <- make_epitope(r[c("84", "85", "86", "87")], COMMON_84)
  e96  <- make_epitope(r[c("96", "170")],            COMMON_96_170, separator = "/")

  epitopes <- tibble(
    epitope  = c("55", "84", "96/170"),
    sequence = c(e55$sequence, e84$sequence, e96$sequence),
    status   = c(e55$status,   e84$status,   e96$status)
  )

  # Step 6: the residue table. Symbols that are not amino acids are replaced
  # with words ("not sequenced", "no residue", "stop codon").
  residue_text <- if_else(
    r %in% NOT_AMINO_ACIDS,
    NOT_AMINO_ACID_LABELS[r],
    r
  )
  residues <- tibble(position = EPITOPE_POSITIONS, residue = unname(residue_text))

  list(
    found    = TRUE,
    message  = "",
    typed    = allele,
    shown    = shown,
    label    = label,
    is_null  = is_null_allele(shown),
    epitopes = epitopes,
    residues = residues
  )
}


# ---- User interface -----------------------------------------------------------
# page_sidebar() gives a sidebar on the left (for choosing the allele) and a
# main area on the right (for the results).
ui <- page_sidebar(
  title = "HLA-DPB1 antibody epitopes",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  sidebar = sidebar(
    width = 320,
    # The choices are filled in by the server (see below), because sending
    # ~3,000 names to the browser at once would be slow.
    selectizeInput(
      "allele",
      "Choose a DPB1 allele (type to search, e.g. 04:01)",
      choices = NULL
    ),
    p(
      class = "text-muted small",
      "Protein sequences: IPD-IMGT/HLA Database version ",
      IPD_IMGT_HLA_VERSION, ". Epitopes: Morris et al., PMID 39277973."
    )
  ),

  # Card 1: the allele name and its three epitopes.
  card(
    card_header("Epitopes"),
    uiOutput("summary"),
    tableOutput("epitope_table")
  ),

  # Card 2: the individual residues the epitopes are built from.
  card(
    card_header("Residues at each position"),
    p(class = "text-muted small",
      "Epitope 55 = positions 55, 56, 57. Epitope 84 = positions 84 to 87. ",
      "Epitope 96/170 = positions 96 and 170."),
    tableOutput("residue_table")
  )
)


# ---- Server -------------------------------------------------------------------
# The server is the part that reacts when the user picks an allele.
server <- function(input, output, session) {

  # Fill the drop-down with every allele. server = TRUE means the browser
  # only receives the names that match what the user is typing.
  updateSelectizeInput(
    session, "allele",
    choices  = all_alleles,
    selected = DEFAULT_ALLELE,
    server   = TRUE
  )

  # Run the epitope function whenever the chosen allele changes.
  result <- reactive({
    req(input$allele)          # wait until an allele is chosen
    dpb1_epitopes(input$allele)
  })

  # The heading: allele name, plus a note for null alleles.
  output$summary <- renderUI({
    res <- result()
    if (!res$found) return(div(class = "alert alert-warning", res$message))
    tagList(
      h4(res$label),
      if (res$is_null) {
        div(class = "alert alert-danger",
            "This is a null allele: it is not expressed on the cell surface, ",
            "so these epitopes are not presented to antibodies.")
      }
    )
  })

  # The epitope table, with friendlier column names.
  output$epitope_table <- renderTable({
    res <- result()
    req(res$found)
    res$epitopes |>
      rename(Epitope = epitope, Sequence = sequence, Status = status)
  })

  # The residue table: one column per position, so it reads left to right.
  output$residue_table <- renderTable({
    res <- result()
    req(res$found)
    residues <- res$residues
    wide <- as.data.frame(t(residues$residue))
    names(wide) <- residues$position
    wide
  })
}


# ---- Start the app -------------------------------------------------------------
shinyApp(ui, server)
