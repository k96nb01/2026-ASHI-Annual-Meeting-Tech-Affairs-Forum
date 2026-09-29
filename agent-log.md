# Building the DPB1 epitope app with an AI coding agent

This is the record of how [`app/app.R`](app/app.R) was built. The agent was a
Claude Code sub-agent (Claude Opus), started and supervised by a second Claude
Code session working with Nick. The agent read the prompt, wrote the code, and
ran it in R to test it. What follows is every message we sent it and every
report it sent back. The reports are copied verbatim; the agent's individual
steps between messages (reading files, running R) are not shown.

This is the third attempt. Each earlier app worked, and each review taught us
something the prompt should have said:

- [`first-attempt/`](first-attempt/): matched allele names with its own text
  rules instead of using immunogenetr, the package built for HLA names; and
  the prompt's table of common epitopes left out EAE.
- [`second-attempt/`](second-attempt/): used immunogenetr, but truncated every
  name in the alignment instead of using `GLstring_regex()`, the function made
  for this; and it listed a null allele's epitopes under a warning that they
  are not expressed.

Each time, the prompt was rewritten (the changes are listed at the end of
[`PROMPT.md`](PROMPT.md)) and a fresh agent, which never saw the earlier
attempts, built the app from scratch.

Time from prompt to a tested app: about 3 minutes. The two rounds of changes
that followed (sections 3 to 7) took about 2 and 5 minutes.

---

## 1. Our message

> Repository: `2026-ASHI-Annual-Meeting-Tech-Affairs-Forum` (R 4.6.1 with
> shiny, bslib, dplyr, stringr, HLAtools 1.11.0, immunogenetr 1.5.0
> installed; use Rscript).
>
> The full task is in PROMPT.md in that repository. Read it and follow it
> exactly. Build from the prompt only; do not read the first-attempt/ or
> second-attempt/ folders.

## 2. The agent's report

> I built `app/app.R` from PROMPT.md, 270 lines. It is the only file I wrote in
> the repo, and nothing was committed. All the checks passed: the named test
> alleles, the whole-alignment check for null alleles, and a headless start of
> the app.
>
> **What it does**
> - **Version:** `ipd_imgt_hla_version <- "3.65.0"` is at the top of the file
>   and is shown in the sidebar. The alignment loads once at startup
>   (`library(HLAtools)`, `alignmentFull(...)[["prot"]][["DPB1"]]`, about 1.7 s,
>   3,170 alleles).
> - **`dpb1_epitopes(allele_typed, alignment = dpb1_alignment)`:** a plain
>   function outside the server. It returns a list: `typed`, `found`, `allele`,
>   `is_null`, `epitopes` (Epitope / Sequence / Common, Uncommon or Not known)
>   and `residues` (the nine positions). Both tables are left empty for null
>   alleles.
> - **Name cleanup:** all done with immunogenetr,
>   `HLA_validate() |> HLA_prefix_remove() |> HLA_prefix_add("DPB1*")`. It trims
>   spaces and accepts the `HLA-`/`DPB1*` prefixes.
> - **Matching:** `str_detect(allele_name, GLstring_regex(...))`, with nothing
>   truncated.
> - **Which allele is shown:** the first expressed match; only if every match is
>   null, the first null. Null is decided by
>   `is.na(GLstring_drop_non_expressed(x, suffixes = "N"))`, so Q alleles count
>   as expressed.
> - **Non-amino-acid symbols:** `*`, `.` and `X` become NA, so any epitope
>   containing one reads "not known", and the residue table says "not known"
>   rather than printing the symbol.
> - **Layout:** bslib `page_sidebar`, with the allele list filled server-side
>   and `DPB1*04:01:01:01` selected at start. For a null allele the app shows a
>   warning and no epitope or residue tables; a name that matches nothing gets a
>   clear "not found" message. Only `|>` is used, the six allowed packages, and
>   "IPD-IMGT/HLA Database" in full.
>
> **Test output**
> ```
> "DPB1*01:01:01:01"  -> DPB1*01:01:01:01 (shown: DPB1*01:01:01:01)  55=AAE  84=DEAV 96/170=K/I  (all Common)
> "DPB1*02:01:02:01"  -> (shown: DPB1*02:01:02:01)  DEE  GGPM  R/T  Common
> "DPB1*04:01:01:01"  -> (shown: DPB1*04:01:01:01)  AAE  GGPM  R/T  Common
> "DPB1*04:02:01:01"  -> (shown: DPB1*04:02:01:01)  DEE  GGPM  R/T  Common
> "DPB1*05:01:01:01"  -> (shown: DPB1*05:01:01:01)  EAE  DEAV  K/I  Common
> " 04:01 "           -> DPB1*04:01 (shown: DPB1*04:01:01:01)  AAE GGPM R/T
> "DPB1*02:01"        -> (shown: DPB1*02:01:02:01)  DEE GGPM R/T
> "04:01:01"          -> (shown: DPB1*04:01:01:01)  AAE GGPM R/T
> "04"                -> DPB1*04 (shown: DPB1*04:01:01:01)  AAE GGPM R/T
> "HLA-DPB1*161:01N"  -> (shown: DPB1*161:01N)  NULL allele, no epitopes
> "61:01"             -> (shown: DPB1*61:01N)   NULL allele (every match is null)
> "61:01N"            -> (shown: DPB1*61:01N)   NULL allele
> "401:01"            -> (shown: DPB1*401:01:01N) NULL allele (every match is null)
> "DPB1*02:01:02:46Q" -> (shown: DPB1*02:01:02:46Q) DEE GGPM R/T Common (Q counts as expressed)
> "99:99"             -> NOT FOUND (searched DPB1*99:99)
> "DPB1*01:01:04"     -> AAE DEAV, 96/170 = not known (457 expressed alleles have an unknown 96/170)
> ```
> The residue table for 04:01 reads A A E / G G P M / R T.
>
> - **Whole-alignment check:** I ran all 2,613 unique two- and three-field names
>   without a suffix (made with `HLA_truncate`). There were 0 problems: every
>   name was found, a null allele came back only when all its matches were null
>   (180 names), and each time the allele shown was the first expressed match in
>   alignment order.
> - **Headless start:** `shiny::testServer` on the app directory gave the
>   expected summary and tables for 04:01:01:01, 61:01, 99:99 and 01:01:04, with
>   no errors. `runApp` in a background R process returned HTTP 200 with the
>   page title.
>
> **Things I was unsure about**
> 1. **Null test:** "null" is taken as the N suffix only (`suffixes = "N"`). The
>    function's default of N/S/C would give the same answer here, since DPB1 in
>    3.65.0 has only N and Q suffixes.
> 2. **Other loci:** `HLA_prefix_remove()` drops any locus, so typing `A*02:01`
>    would silently be searched as `DPB1*02:01`.
> 3. **Something I added:** the allele picker has `create = TRUE`, so the user
>    can press Enter on a short name like `04:01` or `61:01` and see the
>    "(shown: …)" and "not found" behavior in the app itself. Picking from the
>    list shows only full names. Take it out if you want the list to be the only
>    way in.
> 4. **Label wording:** the summary line reads "You chose DPB1\*04:01 (shown:
>    DPB1\*04:01:01:01)". That follows the prompt's example but is not
>    word-for-word.
> 5. **`stringr` use:** it is only for `str_detect()` with the pattern from
>    `GLstring_regex()`, and for `str_c()` to join residues into epitope
>    strings. Neither handles allele names by hand.

**What to notice:** naming the right tools in the prompt made the code
shorter (270 lines, down from 340 in the second attempt) and simpler: one
`GLstring_regex()` call replaced the truncate-and-compare machinery. The
agent's own addition (item 3) is a genuinely good idea for a live demo. But
item 2 is a step backwards: the second agent rejected other genes on its own,
this one did not, because the prompt never asked. A fresh agent only knows
what you write down.

## 3. Our message

> Thanks, this is good. Keep create = TRUE; being able to press Enter on
> "04:01" or "61:01" is great for the demo. One change:
>
> Your item 2 needs fixing: a name typed with a different gene (A\*02:01,
> HLA-DQB1\*04:01) must be rejected with a clear message ("... is not an
> HLA-DPB1 allele name"), never silently searched as DPB1. A name typed with no
> gene at all ("04:01") is still read as DPB1. Do this with immunogenetr (for
> example, HLA_prefix_remove() has a keep_locus argument), not with your own
> text handling. Keep everything else as it is. Re-run the named checks plus
> "A\*02:01", "HLA-DQB1\*04:01", "DPB1\*04:01" and "04:01", confirm the app
> still launches, and report back. No need to redo the whole-alignment check
> unless the lookup itself changes.

## 4. The agent's report

> Names typed with another gene are now rejected, not searched as DPB1. All the
> named checks still pass and the app still launches. One new side effect: a
> lowercase gene (`dpb1*04:01`) is now rejected too, so decide whether that is
> acceptable for the demo (last paragraph).
>
> Everything else is unchanged, including `create = TRUE`. Nothing was
> committed.
>
> **How the check works (in `dpb1_epitopes()`)**
> It uses only immunogenetr calls and plain `!=` comparisons, with no text
> handling of my own:
> ```r
> allele_valid     <- HLA_validate(allele_typed)
> allele_fields    <- HLA_prefix_remove(allele_valid)                     # "04:01"
> allele_with_gene <- HLA_prefix_remove(allele_valid, keep_locus = TRUE)  # "A*02:01" / "DPB1*04:01" / "04:01"
> allele_clean     <- HLA_prefix_add(allele_fields, "DPB1*")
> typed_a_gene <- allele_with_gene != allele_fields
> if (typed_a_gene && allele_with_gene != allele_clean) -> wrong gene
> ```
> - The result list has a new element, `wrong_gene`. When it is TRUE, `typed`
>   holds the name with its gene.
> - The app then shows in red: "**DQB1\*04:01** is not an HLA-DPB1 allele name.
>   Type a DPB1 allele, such as 04:01 or DPB1\*04:01." No tables are shown.
> - The comments explain the check in plain language, and the function header
>   now lists `wrong_gene`.
> - The lookup itself did not change, so I did not redo the whole-alignment
>   check.
>
> **Named checks (Rscript)**
> ```
> "A*02:01"           -> REJECTED: A*02:01 is not an HLA-DPB1 allele name
> "HLA-DQB1*04:01"    -> REJECTED: DQB1*04:01 is not an HLA-DPB1 allele name
> "DPB1*04:01"        -> DPB1*04:01 (shown: DPB1*04:01:01:01)  AAE GGPM R/T, all Common
> "04:01"             -> DPB1*04:01 (shown: DPB1*04:01:01:01)  AAE GGPM R/T, all Common
> ```
> The other 16 checks came out exactly as in the first report:
> - The five full names give the same epitopes as before.
> - " 04:01 ", "04:01:01" and "04" show 04:01:01:01, and "DPB1\*02:01" shows
>   02:01:02:01.
> - 161:01N, 61:01, 61:01N and 401:01 show a null allele with no epitopes.
> - 02:01:02:46Q shows DEE GGPM R/T.
> - 99:99 is not found.
> - DPB1\*01:01:04 has 96/170 not known.
>
> **App**
> - `testServer` gave the expected output for 04:01:01:01, 61:01, 99:99,
>   01:01:04 and HLA-DQB1\*04:01; the last one renders the rejection message with
>   no result tables.
> - `runApp` in a background process returned HTTP 200 with no errors or
>   warnings in the log.
>
> **Decision for you: lowercase gene names**
> `HLA_validate()` keeps the case the user typed, so `dpb1*04:01` now comes back
> as "dpb1\*04:01 is not an HLA-DPB1 allele name". Before this change it was
> quietly read as DPB1. Lowercase was never in the required inputs, and I left
> it alone to stay within "no own text handling". If you want it accepted, the
> simplest fix is to put `stringr::str_to_upper()` on the input before
> `HLA_validate()`.

**What to notice:** we left lowercase names rejected. HLA nomenclature is
written in capitals, and in the app people pick from the list anyway.

## 5. Independent check

A separate agent that had not seen this conversation reviewed the code and
re-checked everything against the alignment:

- All sample alleles (including EAE, a Q allele, and one with an unknown
  96/170) matched its own lookup against the raw alignment.
- It re-ran the lookup for **all 2,613** two- and three-field names against an
  answer it built independently: 0 differences. A null allele came back only
  when every match was null (180 names), and no two-field group lists a null
  allele first when an expressed one exists.
- No false matches (`04:010` and `04:01:0` are not found), other genes
  rejected, empty and garbage input handled without errors.
- All allele-name handling goes through immunogenetr.
- The app starts in about 5 seconds and uses about 380 MB, well within Posit
  Cloud's limit.

It found two problems in the app, both about what the audience would see:

1. **Pressing Enter did not do what the help text said.** Typing `04:01` and
   pressing Enter picked the first full name in the list, not the typed short
   name, so the "(shown: …)" behavior could not be demonstrated. The reviewer
   worked this out by reading the source code of the list widget.
2. **Stray spaces** before punctuation on the page: "(shown:
   DPB1\*04:01:01:01 )".

It also caught a mistake in our own write-up of the first attempt (we had
described its problem wrongly), which we corrected.

## 6. Our message

> Two small fixes found by an independent review; keep everything else as it
> is:
>
> 1. Pressing Enter in the allele box: with selectize's default
>    addPrecedence = FALSE, when the list has matches the first real option is
>    active, so typing "04:01" or "61:01" and pressing Enter picks a full name
>    from the list, not the typed short name. The "(shown: ...)" behavior you
>    built create = TRUE for therefore can't be shown by pressing Enter. Add
>    addPrecedence = TRUE to the selectize options (or otherwise make Enter
>    take the typed text), and make the help text match what actually happens.
>
> 2. Stray spaces in the rendered UI: text children of p() are rendered on
>    separate lines, so the page shows "(shown: DPB1\*04:01:01:01 )" and
>    "version 3.65.0 .". Build those sentences with paste0() so the
>    punctuation sits right.
>
> Re-run testServer and a headless runApp, check the rendered HTML for the
> stray spaces, and if you can, confirm the Enter behavior (from the
> selectize.js source shipped with shiny 1.14.0 if you can't drive a browser).
> Report back.

## 7. The agent's report

> Both review fixes are in. `addPrecedence = TRUE` alone did not fix the demo
> case, so I also added a short, commented piece of JavaScript. The `paste0()`
> fix removed the stray spaces. Everything else is unchanged and nothing was
> committed.
>
> **1. Enter in the allele box**
> - **Why `addPrecedence` was not enough:** the selectize.js shipped with shiny
>   1.14.0 is v0.15.2. In `refreshOptions()` it only falls back to the
>   `addPrecedence` rule when nothing else is highlighted. In single-select
>   mode, the allele already shown is highlighted first if it appears in the
>   filtered list.
> - **What that meant in practice:** with the default DPB1\*04:01:01:01 on
>   screen, typing "04:01" and pressing Enter just re-picked
>   DPB1\*04:01:01:01. I confirmed this in a real browser run. "61:01" did work
>   with `addPrecedence` alone.
> - **The fix:** I kept `addPrecedence = TRUE` and added a small
>   `onInitialize = I("function() {...}")`. It wraps `refreshOptions()`, so
>   each time the list redraws, the "Add <typed>…" line is highlighted again.
>   The comments explain it in plain language, including what JavaScript is.
> - **Result:** Enter always looks up exactly what was typed. Clicking a list
>   item, or arrowing down and pressing Enter, still picks a full name.
> - **Help text is now:** "Type a name, for example 04:01 (the DPB1\* prefix is
>   optional). Press Enter to look up exactly what you typed, or click a full
>   allele name in the list."
>
> **Checked in a real browser.** No chromote or shinytest2 is installed, so I
> wrote a harness in the scratchpad. It runs the app's own ui/server in
> headless Chrome, types into the selectize box, sends a real keydown Enter,
> and records the result. Starting each time from the default
> DPB1\*04:01:01:01:
>
> | Typed + Enter | Selected | App shows |
> |---|---|---|
> | 04:01 | 04:01 | "You chose DPB1\*04:01 (shown: DPB1\*04:01:01:01)" |
> | 61:01 | 61:01 | "You chose DPB1\*61:01 (shown: DPB1\*61:01N)" plus the null-allele warning |
> | A\*02:01 | A\*02:01 | "A\*02:01 is not an HLA-DPB1 allele name…" |
> | 99:99 | 99:99 | "No DPB1 allele matching DPB1\*99:99 was found in the IPD-IMGT/HLA Database, version 3.65.0." |
> | 61:01, then Down, then Enter | DPB1\*61:01N | (list item picked) |
> | 04:01, then Down, then Enter | DPB1\*204:01 | (list item picked) |
>
> The last row happens because the server-side search matches any part of the
> name, so the list for "04:01" starts with 204:01.
>
> **2. Stray spaces**
> - `paste0()` now builds the "(shown: …)" line, the sidebar version note and
>   the end of the "not found" sentence. I also made the help text a single
>   string.
> - The rendered HTML from a headless `runApp` (HTTP 200, no errors or warnings
>   in the log) has "Sequences from the IPD-IMGT/HLA Database, version 3.65.0."
>   and the help text on one line.
> - The browser run shows "(shown: DPB1\*04:01:01:01)" and "version 3.65.0."
>   with no stray space.
> - Where `strong()` sits inside a sentence, the line breaks around it fall on
>   real word gaps, so the text still reads correctly.
>
> **`testServer`** gave the expected output for 04:01:01:01, 61:01, 99:99,
> 01:01:04, HLA-DQB1\*04:01 and 04:01. Only the null, not-found and wrong-gene
> inputs show no result tables, which is correct.
>
> **Decision for you:** this adds about 7 lines of JavaScript to a beginners'
> demo. The alternative is to drop Enter-to-look-up and rely on clicking list
> items, but then the "(shown: …)" behavior can't be shown in the app. I'd keep
> it.

**What to notice:** the reviewer's suggested fix was wrong. It had read the
widget's code but could not run a browser, and the one-line option it proposed
did not work. The builder did not just apply the suggestion: it built a small
test harness, drove the app in a real (headless) browser, found out why, and
fixed it properly. We kept the JavaScript. Checking each other's work, with
real tests, is what makes agents trustworthy.
