# Prompt used to build the DPB1 epitope app

This is the prompt given to an AI coding agent (Claude Code) to build the Shiny
app in [`app/`](app/). The back-and-forth that followed is recorded in
[`agent-log.md`](agent-log.md).

Notice what a good prompt contains: **what** to build, **who** it is for,
**where the data come from**, the **rules** the answer must follow, and **how to
check** that it works. The more of this you write down, the less the agent has
to guess.

This is the third version of the prompt. The earlier versions, and what the
agents built from them, are in [`first-attempt/`](first-attempt/) and
[`second-attempt/`](second-attempt/). Reviewing each app taught us what the
prompt had left out; the changes are listed at the end of this file.

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
handling** (no hand-written regular expressions, `str_detect()` patterns you
build yourself, `startsWith()`, or pasting names together). Read the help
pages (`?GLstring_regex`, `?HLA_validate`, `?HLA_prefix_add`,
`?HLA_prefix_remove`, `?GLstring_drop_non_expressed`) to find the right
function and arguments. In particular:

- Clean up what the user types (stray spaces, with or without the `HLA-` or
  `DPB1*` prefix) with immunogenetr.
- Users will type names with any number of fields (`04:01`, `DPB1*04:01`,
  `04:01:01`, `DPB1*04:01:01:01`). **Do not truncate anything.** Find the
  matching alleles with `GLstring_regex()`, which builds a pattern that
  matches an allele name and any longer name that extends it (so `04:01`
  matches `04:01:01:01` but not `04:010:01`).
- Expression suffixes: `GLstring_regex()` on a name without a suffix also
  matches alleles with a suffix (`61:01` matches `DPB1*61:01N`). When several
  alleles match, show the first **expressed** allele in the alignment (one that
  is not null; Q alleles count as expressed), so `04:01` shows
  `DPB1*04:01:01:01`, not `DPB1*04:01:01:24N`. Only if every match is a null
  allele, show the first null allele. Use `GLstring_drop_non_expressed()` to
  tell whether an allele is null.
- Always tell the user which full allele is shown (for example
  "DPB1\*04:01 (shown: DPB1\*04:01:01:01)"), and give a clear "not found"
  message for names that match nothing.

## What the app does

1. The user picks a DPB1 allele from a searchable list. The list should hold
   every allele in the alignment. Typing `04:01` should find it; the user
   should not have to type the `DPB1*` prefix. (With ~3,000 alleles, use
   server-side selectize.)
2. The app shows, for the selected allele:
   - the full allele name,
   - the three epitopes (55, 84, 96/170),
   - for each epitope, whether it is one of the common sequences in the table
     above, is an uncommon sequence, or is not known.
   - **Null alleles have no epitopes**: the protein is not expressed, so there
     is nothing for antibodies to see. For a null allele, do not list any
     epitopes; say clearly that it is a null allele and not expressed.
3. Below that (for expressed alleles only), a small table of the residues at
   each individual position used (55, 56, 57, 84, 85, 86, 87, 96, 170), so the
   audience can see how the epitope strings are built.
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
   `DPB1*02:01`, `04:01:01`, `04`, `HLA-DPB1*161:01N`, `61:01`, `61:01N`,
   `401:01`, `DPB1*02:01:02:46Q`, `99:99`, and one allele whose 96/170
   epitope is not known.
2. Check across the whole alignment that a two- or three-field name without a
   suffix returns a null allele only when every allele it matches is null.
3. Start the app headlessly (for example with `shiny::runApp()` in a
   background R process, or `shiny::testServer()`) and confirm it launches
   without errors.
4. Report what you built, the test output, and anything you were unsure about.

---

## What changed from the second version

- **`GLstring_regex()` instead of truncating.** The second app found alleles
  by truncating every name in the alignment and comparing. immunogenetr has a
  function made for finding an allele from a name with any number of fields,
  so the prompt now names it.
- **Null alleles show no epitopes.** The second app listed a null allele's
  epitopes and then warned that they are not expressed, which contradicted
  itself. A null allele is now shown as null, with no epitopes.
- **`GLstring_drop_non_expressed()`** for spotting null alleles is now in the
  prompt; the second agent had to be told about it.

## What changed from the first version

- **immunogenetr for allele names.** The first app matched short names with
  its own text rule (the typed name followed by `:`). immunogenetr already
  handles HLA names, including expression suffixes, so the second prompt told
  the agent to use it (by truncating names with `HLA_truncate()`, an approach
  the third version replaced with `GLstring_regex()`).
- **EAE added** to the common 55 epitopes.
- **`X` (stop codon)** added to the symbols that are not amino acids; the
  first agent found it on its own.
- **`library(HLAtools)`**: the first prompt's example call did not work
  without it; the first agent found this too.
- **Two- and three-field names and stray spaces** are now in the prompt from
  the start instead of being asked for afterwards.
