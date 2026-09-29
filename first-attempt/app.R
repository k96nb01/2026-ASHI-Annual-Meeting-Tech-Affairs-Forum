# =============================================================================
# HLA-DPB1 epitope viewer
#
# Pick any HLA-DPB1 allele and see its three major antibody epitopes
# (55, 84 and 96/170; Morris et al., PMID 39277973).
#
# The protein sequences come from the IPD-IMGT/HLA Database, downloaded with
# the HLAtools package when the app starts.
# =============================================================================

# --- Packages -----------------------------------------------------------------
# Each library() call loads a set of ready-made tools.
library(shiny)         # the web-app framework
library(bslib)         # a clean, modern page layout
library(dplyr)         # tools for working with tables
library(stringr)       # tools for working with text
library(HLAtools)      # downloads alignments from the IPD-IMGT/HLA Database
library(immunogenetr)  # HLA helpers, such as adding or removing "DPB1*"
# Note: HLAtools must be loaded with library(), not only called as
# HLAtools::alignmentFull(). Its download function looks up a built-in table
# of loci that is only found once the package is loaded.


# --- Settings -----------------------------------------------------------------
# The IPD-IMGT/HLA Database is updated every three months. Fixing (pinning)
# the version here means the same allele always gives the same answer, no
# matter when the app is run. Change this one line to use a newer release.
IPD_IMGT_HLA_VERSION <- "3.65.0"

# The allele shown when the app first opens.
DEFAULT_ALLELE <- "DPB1*04:01:01:01"


# --- The epitope definitions --------------------------------------------------
# Each epitope is built from the amino acids at a few positions in the
# protein. "common" lists the sequences most often seen for each epitope.
# 96/170 is written as two residues separated by a slash, e.g. "K/I".
EPITOPES <- list(
  "55"     = list(positions = c("55", "56", "57"),       separator = "",  common = c("AAE", "DEE", "DED")),
  "84"     = list(positions = c("84", "85", "86", "87"), separator = "",  common = c("DEAV", "GGPM", "VGPM")),
  "96/170" = list(positions = c("96", "170"),            separator = "/", common = c("K/I", "R/T"))
)

# Every position we need, in order: 55, 56, 57, 84, 85, 86, 87, 96, 170.
ALL_POSITIONS <- unlist(lapply(EPITOPES, function(e) e$positions), use.names = FALSE)

# Symbols in the alignment that are NOT real amino acids:
#   "*" = the residue is not known (that part of the sequence was never determined)
#   "." = an insertion/deletion placeholder (no residue at this position)
#   "X" = a stop codon; the protein ends here (mostly seen in null alleles)
# An epitope that contains any of these is reported as "not known".
NOT_AMINO_ACIDS <- c("*", ".", "X")


# --- Load the data (once, when the app starts) --------------------------------
# This downloads the DPB1 protein alignment from the IPD-IMGT/HLA Database
# (about 2 seconds). It gives one row per allele and one column per amino
# acid position, named by position number ("55", "84", "170", ...).
# Because this code sits outside the server function, it runs only once,
# not every time someone picks an allele.
dpb1_alignment <- alignmentFull(
  loci = "DPB1", alignType = "prot", version = IPD_IMGT_HLA_VERSION
)[["prot"]][["DPB1"]] |>
  select(allele_name, all_of(ALL_POSITIONS))

# The full list of allele names for the search box, e.g. "DPB1*04:01:01:01".
all_alleles <- dpb1_alignment$allele_name


# --- The epitope function -----------------------------------------------------
# This small function does all of the HLA work. It lives outside the app so
# you can try it on its own in the R console, for example:
#   dpb1_epitopes("DPB1*04:01:01:01")
#   dpb1_epitopes("04:01:01:01")          # the "DPB1*" prefix is optional
#   dpb1_epitopes("04:01")                # a shorter name uses the first
#                                         # matching allele, DPB1*04:01:01:01
dpb1_epitopes <- function(allele) {

  # Tidy up the name: remove any "HLA-" or "DPB1*" prefix, then add "DPB1*"
  # back, so "04:01:01:01", "DPB1*04:01:01:01" and "HLA-DPB1*04:01:01:01"
  # all become "DPB1*04:01:01:01", matching the alignment.
  allele <- allele |>
    HLA_prefix_remove() |>
    HLA_prefix_add(prefix = "DPB1*")

  # Find this allele's row in the alignment. First look for an exact match.
  row <- dpb1_alignment |> filter(allele_name == allele)

  # No exact match? People often type a shorter name, such as "04:01" (two
  # fields) or "04:01:01" (three fields). In that case, use the first allele
  # in the alignment whose name starts with that name followed by ":".
  # The ":" matters: "04:01" should find "DPB1*04:01:01:01", not "DPB1*04:010".
  if (nrow(row) == 0) {
    row <- dpb1_alignment |>
      filter(str_starts(allele_name, fixed(str_c(allele, ":")))) |>
      slice(1)
  }

  # Still nothing? Then the name is not in the alignment at all.
  if (nrow(row) == 0) {
    stop("Allele not found in the IPD-IMGT/HLA Database ",
         IPD_IMGT_HLA_VERSION, " alignment: ", allele)
  }

  # The full allele name actually used, e.g. "DPB1*04:01:01:01".
  allele_used <- row$allele_name[1]

  # The residue (amino acid letter) at each position we care about.
  residues <- unlist(row[1, ALL_POSITIONS])

  # Build each epitope and decide whether it is common, uncommon or not known.
  epitopes <- lapply(names(EPITOPES), function(name) {
    e <- EPITOPES[[name]]
    letters_here <- residues[e$positions]

    if (any(letters_here %in% NOT_AMINO_ACIDS)) {
      # At least one position is not a real amino acid: don't guess.
      sequence <- "not known"
      status   <- "not known"
    } else {
      # Join the letters, e.g. D + E + A + V -> "DEAV", or K + I -> "K/I".
      sequence <- str_c(letters_here, collapse = e$separator)
      status   <- if_else(sequence %in% e$common, "common", "uncommon")
    }

    data.frame(epitope = name, sequence = sequence, status = status)
  }) |>
    bind_rows()

  # The residue table: replace the non-amino-acid symbols with "?" so they
  # can never be mistaken for real amino acids.
  residue_table <- data.frame(
    position = names(residues),
    residue  = if_else(residues %in% NOT_AMINO_ACIDS, "?", residues)
  )

  # Send everything back as one list.
  # A label for display: just the name for an exact match, otherwise the
  # name typed plus the full allele shown, e.g.
  # "DPB1*04:01 (shown: DPB1*04:01:01:01)".
  label <- if_else(
    allele_used == allele, allele,
    str_c(allele, " (shown: ", allele_used, ")")
  )

  list(
    allele      = allele,        # the name asked for (with "DPB1*" added)
    allele_used = allele_used,   # the full allele whose sequence is shown
    label       = label,
    is_null     = str_detect(allele_used, "N$"),  # null alleles end in "N"
    epitopes    = epitopes,
    residues    = residue_table
  )
}


# --- The user interface (what people see) -------------------------------------
ui <- page_sidebar(
  title = "HLA-DPB1 epitope viewer",

  # The left-hand panel: choose an allele.
  sidebar = sidebar(
    # The list of ~3,000 alleles is filled in by the server (see below),
    # which keeps the page fast.
    selectizeInput(
      "allele", "Choose a DPB1 allele",
      choices = NULL,
      options = list(placeholder = "Type e.g. 04:01")
    ),
    p("Type part of a name, e.g. ", code("04:01"),
      "; there is no need to type ", code("DPB1*"), "."),
    hr(),
    p(class = "small text-muted",
      "Data: IPD-IMGT/HLA Database version ", IPD_IMGT_HLA_VERSION, ".", br(),
      "Epitopes: Morris et al., PMID 39277973.")
  ),

  # The main area: the results.
  card(
    card_header("Allele"),
    h3(textOutput("allele_name")),
    uiOutput("null_note")
  ),
  card(
    card_header("Major antibody epitopes"),
    tableOutput("epitope_table")
  ),
  card(
    card_header("Residues at each position"),
    tableOutput("residue_table"),
    p(class = "small text-muted",
      "? = not known in the IPD-IMGT/HLA Database alignment ",
      "(sequence not determined, a gap, or after a stop codon).")
  )
)


# --- The server (what the app does) -------------------------------------------
server <- function(input, output, session) {

  # Fill the allele list on the server side ("server = TRUE"). The browser
  # only receives the names that match what the user is typing.
  updateSelectizeInput(
    session, "allele",
    choices = all_alleles, selected = DEFAULT_ALLELE, server = TRUE
  )

  # Re-run the epitope function whenever a new allele is picked.
  result <- reactive({
    req(input$allele)  # wait until an allele has been chosen
    dpb1_epitopes(input$allele)
  })

  output$allele_name <- renderText(result()$label)

  # Show a warning only for null alleles.
  output$null_note <- renderUI({
    if (result()$is_null) {
      div(class = "alert alert-warning mb-0",
          "This is a null allele (its name ends in N): the protein is not ",
          "expressed on the cell surface, so its epitopes are not expected ",
          "to be seen by antibodies.")
    }
  })

  output$epitope_table <- renderTable(
    result()$epitopes |>
      rename(Epitope = epitope, Sequence = sequence, Status = status)
  )

  # Show the residues as one row: positions across the top.
  output$residue_table <- renderTable({
    r <- result()$residues
    as.data.frame(as.list(setNames(r$residue, r$position)), check.names = FALSE)
  })
}


# --- Start the app ------------------------------------------------------------
shinyApp(ui, server)
