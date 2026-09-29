# Prompt used to build the DPB1 epitope app

*Second attempt, kept for comparison. The app it produced is [`app.R`](app.R) in this folder; the current prompt is [`../PROMPT.md`](../PROMPT.md). Links below refer to the original layout.*

This is the prompt given to an AI coding agent (Claude Code) to build the Shiny
app in [`app/`](app/). The back-and-forth that followed is recorded in
[`agent-log.md`](agent-log.md).

Notice what a good prompt contains: **what** to build, **who** it is for,
**where the data come from**, the **rules** the answer must follow, and **how to
check** that it works. The more of this you write down, the less the agent has
to guess.

This is the second version of the prompt. The first version, and what the
agent built from it, are in [`first-attempt/`](first-attempt/). Reviewing that
first app taught us what the prompt had left out; the changes are listed at
the end of this file.

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
| 55 | 55, 56, 57 | AAE, DEE, DED, EAE |
| 84 | 84, 85, 86, 87 | DEAV, GGPM, VGPM |
| 96/170 | 96 and 170 | K/I, R/T |

Write the 55 and 84 epitopes as the residues run together (for example
`DEAV`) and the 96/170 epitope as the two residues separated by a slash (for
example `K/I`).

## Data

- Get the DPB1 protein alignment from the IPD-IMGT/HLA Database with the
  HLAtools package. Load the package with `library(HLAtools)` (its functions
  need the package's built-in data, so calling them only as `HLAtools::` does
  not work), then call
  `alignmentFull(loci = "DPB1", alignType = "prot", version = "3.65.0")`.
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
  determined), `.` means there is no residue at that position (an
  insertion/deletion placeholder, or after the protein has ended), and `X`
  marks a stop codon. None of these is an amino acid: never show them as if
  they were, and treat an epitope containing any of them as "not known".

## Use immunogenetr for allele names

Handling HLA allele names is exactly what the immunogenetr package is for, and
it already gets the tricky cases right. **Wherever the code works with an
allele name, use an immunogenetr function rather than writing your own text
handling** (no hand-written regular expressions or `str_detect()` /
`startsWith()` on allele names). Read the help pages (`?HLA_truncate`,
`?HLA_validate`, `?HLA_prefix_add`, `?HLA_prefix_remove`) to find the right
function and arguments. In particular:

- Clean up what the user types (stray spaces, with or without the `HLA-` or
  `DPB1*` prefix) with immunogenetr.
- Users will often type a two- or three-field name (`04:01`, `DPB1*04:01`,
  `04:01:01`). Find the matching alleles by truncating the alignment's allele
  names to the same number of fields with `HLA_truncate()` and comparing the
  names, not by matching text. When several alleles match, show the first one
  in the alignment, and tell the user which full allele is shown (for example
  "DPB1\*04:01 (shown: DPB1\*04:01:01:01)").
- Expression suffixes matter: `DPB1*61:01N` is a different allele from
  `DPB1*61:01`, and a null allele must never be returned for a name that did
  not ask for one. Make sure typing `04:01` can never land on
  `DPB1*04:01:01:24N`, and that typing a name with no expressed match gives a
  clear "not found" message rather than a null allele.

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
   - a note if the allele is a null allele.
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
  immunogenetr.
- Keep it short and readable; clarity beats cleverness. Put the epitope logic
  in one small plain R function (outside the server) that takes an allele name
  and returns the result, so it can be demonstrated and tested on its own in
  the console.
- Write "IPD-IMGT/HLA Database" in full wherever the data source is named in the
  app or comments; never shorten it to "IMGT".
- A clean, simple layout with bslib (`page_sidebar()` or similar).

## How to check it works

Before reporting back:

1. Run the epitope function in R (`Rscript`) and report the output for:
   `DPB1*01:01:01:01`, `DPB1*02:01:02:01`, `DPB1*04:01:01:01`,
   `DPB1*04:02:01:01`, `DPB1*05:01:01:01`, `" 04:01 "` (with spaces),
   `DPB1*02:01`, `04:01:01`, `HLA-DPB1*161:01N`, `61:01`, `61:01N`,
   `401:01`, `99:99`, and one allele whose 96/170 epitope is not known.
2. Check across the whole alignment that no two- or three-field name without
   an expression suffix ever returns a null allele.
3. Start the app headlessly (for example with `shiny::runApp()` in a
   background R process, or `shiny::testServer()`) and confirm it launches
   without errors.
4. Report what you built, the test output, and anything you were unsure about.

---

## What changed from the first version

- **immunogenetr for allele names.** The first app matched short names by
  text (the typed name followed by `:`). That let `401:01` land on the null
  allele `DPB1*401:01:01N`. `HLA_truncate()` already keeps expression
  suffixes, so comparing truncated names avoids the problem without any new
  code.
- **EAE added** to the common 55 epitopes.
- **`X` (stop codon)** added to the symbols that are not amino acids; the
  first agent found it on its own.
- **`library(HLAtools)`**: the first prompt's example call did not work
  without it; the first agent found this too.
- **Two- and three-field names and stray spaces** are now in the prompt from
  the start instead of being asked for afterwards.
