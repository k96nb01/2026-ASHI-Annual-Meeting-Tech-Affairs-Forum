# =============================================================================
# HLA-DPB1 epitope viewer
#
# Pick any HLA-DPB1 allele and see its three major antibody epitopes
# (55, 84 and 96/170), as defined by Morris et al. (PMID 39277973).
#
# The protein sequences come from the IPD-IMGT/HLA Database, downloaded with
# the HLAtools package. Allele names are handled with the immunogenetr package.
#
# To run: open this file in RStudio and click "Run App", or run
#   shiny::runApp("app")
# =============================================================================


# --- 1. Packages -------------------------------------------------------------
# Each library() call loads a package: a bundle of ready-made functions.
library(shiny)         # builds the web app
library(bslib)         # gives the app a clean, modern layout
library(dplyr)         # tools for working with data tables
library(stringr)       # tools for working with text
library(HLAtools)      # downloads HLA alignments from the IPD-IMGT/HLA Database
library(immunogenetr)  # understands HLA allele names


# --- 2. Settings --------------------------------------------------------------
# The IPD-IMGT/HLA Database releases a new version every few months, and
# sequences can change between versions. Fixing the version here means the
# same allele always gives the same answer. To use a newer release, change
# this one line.
ipd_imgt_hla_version <- "3.65.0"

# The allele shown when the app first opens, so the demo starts with a result.
default_allele <- "DPB1*04:01"

# The common sequence of each epitope (Morris et al., PMID 39277973).
# Anything else we find is called "uncommon".
common_epitopes <- list(
  "55"     = c("AAE", "DEE", "DED", "EAE"),
  "84"     = c("DEAV", "GGPM", "VGPM"),
  "96/170" = c("K/I", "R/T")
)

# The individual amino acid positions that make up the three epitopes.
epitope_positions <- c("55", "56", "57", "84", "85", "86", "87", "96", "170")

# Symbols in the alignment that are NOT amino acids:
#   "*" = the residue is not known (that part of the sequence was never determined)
#   "." = there is no residue here (a gap, or the protein has already ended)
#   "X" = a stop codon (the protein ends here)
# We never show these as if they were amino acids.
not_amino_acids <- c("*", ".", "X")


# --- 3. Load the data (once, when the app starts) ----------------------------
# alignmentFull() downloads the DPB1 protein alignment (about 2 seconds).
# It returns a list; the table we want is inside it at [["prot"]][["DPB1"]].
# The table has one row per allele, and one column per amino acid position.
# Because this runs outside the server function below, it happens only once,
# not every time someone picks an allele.
message("Downloading the DPB1 protein alignment from the IPD-IMGT/HLA Database, version ",
        ipd_imgt_hla_version, " ...")
dpb1_alignment <- alignmentFull(loci = "DPB1", alignType = "prot",
                                version = ipd_imgt_hla_version)[["prot"]][["DPB1"]]

# The list of alleles to choose from. Many alleles differ only in the third or
# fourth field (changes that don't alter the protein), so we shorten every name
# to two fields with HLA_truncate() and keep one of each with distinct().
# keep_suffix = FALSE drops the N/Q suffix, so "04:01" is listed once rather
# than also as "04:01N"; the epitope function below still works out whether a
# name is a null allele.
dpb1_choices <- dpb1_alignment |>
  mutate(two_field = HLA_truncate(allele_name, keep_suffix = FALSE)) |>
  distinct(two_field) |>
  pull(two_field)


# --- 4. The epitope function -------------------------------------------------
# This is the heart of the app. It is an ordinary R function, so you can try
# it on its own in the console, for example:
#   dpb1_epitopes("04:01")
#
# It takes an allele name as typed by a person and returns a list with:
#   typed     - the cleaned-up name we searched for
#   wrong_gene - TRUE if the name was typed with a gene other than DPB1
#   found     - TRUE if any allele in the alignment matched
#   allele    - the full name of the allele shown
#   is_null   - TRUE if that allele is a null (not expressed) allele
#   epitopes  - a table of the three epitopes (empty for a null allele)
#   residues  - a table of the amino acid at each position (empty for a null allele)
dpb1_epitopes <- function(allele_typed, alignment = dpb1_alignment) {

  # Step 1: clean up what the user typed, using immunogenetr.
  #   HLA_validate() removes stray spaces and returns NA for text that is not
  #     an allele name at all.
  allele_valid <- HLA_validate(allele_typed)

  # This is the empty answer we give when nothing matches.
  not_found <- list(typed = allele_typed, wrong_gene = FALSE, found = FALSE,
                    allele = NA_character_, is_null = NA,
                    epitopes = NULL, residues = NULL)
  if (is.na(allele_valid)) {
    return(not_found)
  }

  # Step 1b: make sure the name is for the DPB1 gene.
  #   HLA_prefix_remove() strips the "HLA-" prefix and the gene, leaving only
  #     the fields ("04:01"). With keep_locus = TRUE it strips only "HLA-" and
  #     keeps the gene ("DPB1*04:01", "A*02:01", or just "04:01" if the user
  #     typed no gene at all).
  #   HLA_prefix_add() puts "DPB1*" in front of the fields, giving the name
  #     written the same way as in the alignment.
  # If the user typed a gene, the name with its gene must be the same as the
  # DPB1 name; otherwise it belongs to another gene (such as A*02:01 or
  # DQB1*04:01) and we must not look it up as if it were DPB1.
  allele_fields    <- HLA_prefix_remove(allele_valid)
  allele_with_gene <- HLA_prefix_remove(allele_valid, keep_locus = TRUE)
  allele_clean     <- HLA_prefix_add(allele_fields, "DPB1*")

  typed_a_gene <- allele_with_gene != allele_fields
  if (typed_a_gene && allele_with_gene != allele_clean) {
    not_found$typed <- allele_with_gene
    not_found$wrong_gene <- TRUE
    return(not_found)
  }
  not_found$typed <- allele_clean

  # Step 2: find every allele in the alignment that matches the name.
  # GLstring_regex() builds a search pattern that matches the allele and any
  # longer name that extends it: "04:01" matches "04:01:01:01", but not
  # "04:010:01". No part of the name is cut off or ignored.
  pattern <- GLstring_regex(allele_clean)
  matches <- alignment |>
    filter(str_detect(allele_name, pattern))

  if (nrow(matches) == 0) {
    return(not_found)
  }

  # Step 3: pick which matching allele to show.
  # A name without a suffix also matches null alleles ("61:01" matches
  # "61:01N"). We prefer the first EXPRESSED allele in the alignment.
  # GLstring_drop_non_expressed() removes null ("N") alleles: it returns NA
  # for a null allele and the name itself for an expressed one. (Q alleles,
  # "questionable" expression, are kept, so they count as expressed.)
  # Only if every match is null do we show the first null allele.
  matches <- matches |>
    mutate(is_null = is.na(GLstring_drop_non_expressed(allele_name, suffixes = "N")))
  expressed <- matches |> filter(!is_null)
  chosen <- if (nrow(expressed) > 0) expressed[1, ] else matches[1, ]

  answer <- list(typed = allele_clean, wrong_gene = FALSE, found = TRUE,
                 allele = chosen$allele_name, is_null = chosen$is_null,
                 epitopes = NULL, residues = NULL)

  # Step 4: a null allele makes no protein, so antibodies have nothing to
  # see. We stop here and report no epitopes.
  if (chosen$is_null) {
    return(answer)
  }

  # Step 5: read the amino acid at each epitope position.
  # Any "*", "." or "X" is turned into NA, which R uses for "not known".
  residues <- unlist(chosen[1, epitope_positions])
  residues[residues %in% not_amino_acids] <- NA

  # Step 6: build the three epitope sequences.
  # 55 and 84 are the residues run together ("DEAV"); 96/170 is two residues
  # with a slash between them ("K/I"). If any residue is not known, the whole
  # epitope is not known.
  build <- function(positions, separator) {
    if (any(is.na(residues[positions]))) return(NA_character_)
    str_c(residues[positions], collapse = separator)
  }
  sequences <- c(
    "55"     = build(c("55", "56", "57"), ""),
    "84"     = build(c("84", "85", "86", "87"), ""),
    "96/170" = build(c("96", "170"), "/")
  )

  # Step 7: label each epitope as common, uncommon or not known, by
  # comparing it with the list of common sequences at the top of the file.
  label <- function(epitope) {
    sequence <- sequences[[epitope]]
    if (is.na(sequence)) return("Not known")
    if (sequence %in% common_epitopes[[epitope]]) return("Common")
    "Uncommon"
  }
  answer$epitopes <- tibble(
    Epitope  = names(sequences),
    Sequence = if_else(is.na(sequences), "not known", sequences),
    Type     = c(label("55"), label("84"), label("96/170"))
  )

  # Step 8: a table of the single residues, so you can see how the epitope
  # strings above are put together.
  answer$residues <- tibble(
    Position   = names(residues),
    `Amino acid` = if_else(is.na(residues), "not known", residues)
  )

  answer
}


# --- 5. The user interface (what the audience sees) ---------------------------
ui <- page_sidebar(
  title = "HLA-DPB1 antibody epitopes",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  # The sidebar holds the allele picker and some background notes.
  sidebar = sidebar(
    width = 320,
    # The choices are filled in by the server (see below), because a list of
    # about 2,000 alleles is too long to send to the browser all at once.
    selectizeInput("allele", "Choose a DPB1 allele",
                   choices = NULL,
                   options = list(
                     placeholder = "Type e.g. 04:01",
                     # create = TRUE lets the user look up a name that is not
                     # in the list as typed, such as "04:01" or "61:01". The
                     # list then offers "Add 04:01..." as its first line.
                     create = TRUE,
                     # addPrecedence = TRUE makes that "Add ..." line the one
                     # Enter picks, instead of the first allele in the list.
                     addPrecedence = TRUE,
                     # One case addPrecedence misses: if the allele already
                     # shown (say DPB1*04:01:01:01) is in the list, the box
                     # highlights it and Enter would pick it again. This small
                     # piece of JavaScript (the language of web pages) runs
                     # once when the box is set up. Every time the list is
                     # redrawn, it moves the highlight back to the "Add ..."
                     # line, so Enter always looks up exactly what was typed.
                     onInitialize = I("function() {
                       var box = this, redraw = box.refreshOptions;
                       box.refreshOptions = function() {
                         redraw.apply(box, arguments);
                         var typed = box.$dropdown_content.find('.create');
                         if (typed.length) box.setActiveOption(typed);
                       };
                     }")
                   )),
    helpText(paste0("Type a name, for example 04:01 (the DPB1* prefix is ",
                    "optional). Press Enter to look up exactly what you ",
                    "typed, or click an allele in the list.")),
    hr(),
    p(strong("Epitopes"), " (Morris et al., PMID 39277973):"),
    tags$ul(
      tags$li("55: positions 55, 56, 57"),
      tags$li("84: positions 84, 85, 86, 87"),
      tags$li("96/170: positions 96 and 170")
    ),
    # paste0() glues the pieces into one sentence with no extra spaces.
    p(class = "small text-muted",
      paste0("Sequences from the IPD-IMGT/HLA Database, version ",
             ipd_imgt_hla_version, "."))
  ),

  # The main panel shows the results.
  card(
    card_header("Selected allele"),
    uiOutput("allele_summary")
  ),
  uiOutput("result_cards")
)


# --- 6. The server (what the app does when the user acts) ---------------------
server <- function(input, output, session) {

  # Fill the allele list. server = TRUE means the browser only receives the
  # alleles that match what the user is typing ("server-side selectize").
  updateSelectizeInput(session, "allele",
                       choices = dpb1_choices,
                       selected = default_allele,
                       server = TRUE)

  # Run the epitope function whenever the chosen allele changes.
  result <- reactive({
    req(input$allele)
    dpb1_epitopes(input$allele)
  })

  # The top card: which allele is shown, or a clear message if none is.
  output$allele_summary <- renderUI({
    r <- result()
    if (r$wrong_gene) {
      return(p(class = "text-danger",
               strong(r$typed), " is not an HLA-DPB1 allele name. ",
               "Type a DPB1 allele, such as 04:01 or DPB1*04:01."))
    }
    if (!r$found) {
      return(p(class = "text-danger",
               "No DPB1 allele matching ", strong(r$typed),
               paste0(" was found in the IPD-IMGT/HLA Database, version ",
                      ipd_imgt_hla_version, ".")))
    }
    tagList(
      h3(r$allele),
      p(class = "text-muted",
        paste0("You chose ", r$typed, " (shown: ", r$allele, ")")),
      if (r$is_null) {
        div(class = "alert alert-warning",
            strong("This is a null allele. "),
            "It is not expressed: no protein reaches the cell surface, ",
            "so there are no epitopes for antibodies to recognize.")
      }
    )
  })

  # The epitope and residue tables: only for expressed alleles.
  output$result_cards <- renderUI({
    r <- result()
    if (!r$found || r$is_null) return(NULL)
    layout_columns(
      col_widths = c(7, 5),
      card(card_header("Epitopes"), tableOutput("epitope_table")),
      card(card_header("Amino acid at each position"), tableOutput("residue_table"))
    )
  })

  output$epitope_table <- renderTable(result()$epitopes)
  output$residue_table <- renderTable(result()$residues)
}


# --- 7. Start the app ---------------------------------------------------------
shinyApp(ui, server)
