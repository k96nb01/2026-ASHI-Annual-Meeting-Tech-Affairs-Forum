# Prompt used to build the DPB1 epitope app

*First attempt, kept for comparison. The app it produced is [`app.R`](app.R) in this folder; the current prompt is [`../PROMPT.md`](../PROMPT.md). Links below refer to the original layout.*

This is the prompt given to an AI coding agent (Claude Code) to build the Shiny
app in [`app/`](app/). The back-and-forth that followed is recorded in
[`agent-log.md`](agent-log.md).

Notice what a good prompt contains: **what** to build, **who** it is for,
**where the data come from**, the **rules** the answer must follow, and **how to
check** that it works. The more of this you write down, the less the agent has
to guess.

---

Build a small R Shiny app that shows the major antibody epitopes of any
HLA-DPB1 allele. It will be used as a live demo in a 40-minute introductory
talk to HLA laboratory technologists who have little or no coding experience,
and run on Posit Cloud.

## Where to put it

- Write a single file: `app/app.R` in this repository.
- Do not create any other files, and do not commit or push anything.

## Background

DPB1 has three variable regions that define the major functional antibody
epitopes (Morris et al., PMID 39277973):

| Epitope | Amino acid positions | Common sequences |
|---|---|---|
| 55 | 55, 56, 57 | AAE, DEE, DED |
| 84 | 84, 85, 86, 87 | DEAV, GGPM, VGPM |
| 96/170 | 96 and 170 | K/I, R/T |

Write the 55 and 84 epitopes as the residues run together (for example
`DEAV`) and the 96/170 epitope as the two residues separated by a slash (for
example `K/I`).

## Data

- Get the DPB1 protein alignment from the IPD-IMGT/HLA Database with the
  HLAtools package:
  `HLAtools::alignmentFull(loci = "DPB1", alignType = "prot", version = "3.65.0")`.
  This downloads the alignment (about 2 seconds) and returns a list; the data
  frame is at `[["prot"]][["DPB1"]]`. It has one row per allele, the columns
  `locus`, `allele`, `trimmed_allele` (two-field name with the `DPB1*` prefix),
  and `allele_name` (full name), followed by one column per amino acid
  position, named by position number (`"55"`, `"84"`, `"170"`, ...).
- Pin the database version (3.65.0) in one clearly named variable at the top
  of the file, and show the version in the app, so the same input always gives
  the same answer.
- Load the alignment once, when the app starts, not on every click.
- In the alignment, `*` means the residue is not known (the sequence was never
  determined) and `.` means an insertion/deletion placeholder. Never show these
  as if they were real amino acids: an epitope containing either is "not known".

## What the app does

1. The user picks a DPB1 allele from a searchable list. The list should hold
   every allele in the alignment. Typing `04:01` should find it; the user
   should not have to type the `DPB1*` prefix. (With ~3,000 alleles, use
   server-side selectize.)
2. The app shows, for the selected allele:
   - the full allele name,
   - the three epitopes (55, 84, 96/170),
   - for each epitope, whether it is one of the common sequences in the table
     above, is an uncommon sequence, or is not known,
   - a note if the allele is a null allele (name ends in `N`).
3. Below that, a small table of the residues at each individual position used
   (55, 56, 57, 84, 85, 86, 87, 96, 170), so the audience can see how the
   epitope strings are built.
4. Start with a sensible default allele selected (DPB1\*04:01:01:01) so the
   demo shows something immediately.

## Rules for the code

- The audience is new to coding. Comment generously, in plain language:
  explain what each block does and why, not only how.
- Use the base R pipe `|>`, not `%>%`.
- Use only these packages: shiny, bslib, dplyr, stringr, HLAtools,
  immunogenetr. Use immunogenetr where it naturally fits (for example
  `HLA_prefix_add()` / `HLA_prefix_remove()` to handle the `HLA-` / `DPB1*`
  prefixes), but do not force it.
- Keep it short and readable; clarity beats cleverness. Put the epitope logic
  in one small plain R function (outside the server) that takes an allele name
  and returns the result, so it can be demonstrated and tested on its own in
  the console.
- Write "IPD-IMGT/HLA Database" in full wherever the data source is named in the
  app or comments; never shorten it to "IMGT".
- A clean, simple layout with bslib (`page_sidebar()` or similar).

## How to check it works

Before reporting back:

1. Run the epitope function in R (`Rscript`) for at least these alleles and
   report the output: DPB1\*01:01:01:01, DPB1\*02:01:02:01, DPB1\*04:01:01:01,
   DPB1\*04:02:01:01, one allele whose 96/170 epitope is not known, and one
   null allele.
2. Start the app headlessly (for example with `shiny::runApp()` in a
   background R process, or `shiny::testServer()`) and confirm it launches
   without errors.
3. Report what you built, the test output, and anything you were unsure about.
