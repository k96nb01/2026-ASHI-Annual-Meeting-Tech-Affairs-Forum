# Building the DPB1 epitope app with an AI coding agent

*Second attempt, kept for comparison. The app it produced is [`app.R`](app.R) in this folder; the current prompt is [`../PROMPT.md`](../PROMPT.md). Links below refer to the original layout.*

This is the record of how [`app/app.R`](app/app.R) was built. The agent was a
Claude Code sub-agent (Claude Opus), started and supervised by a second Claude
Code session working with Nick. The agent read the prompt, wrote the code, and
ran it in R to test it. What follows is every message we sent it and every
report it sent back. The reports are copied verbatim; the agent's individual
steps between messages (reading files, running R) are not shown.

This is the second attempt. The first prompt, the first agent's log and its
app are in [`first-attempt/`](first-attempt/). That app worked, but reviewing
it showed it matched allele names by text instead of using immunogenetr, which
let a name like `401:01` land on a null allele. The prompt was rewritten
(the changes are listed at the end of [`PROMPT.md`](PROMPT.md)) and a fresh
agent, which never saw the first attempt, built the app again.

Time from prompt to a tested app: about 13 minutes, most of it spent on a
check the prompt asked for across every name in the alignment. The change in
section 3 took a few more minutes.

---

## 1. Our message

> Repository: `2026-ASHI-Annual-Meeting-Tech-Affairs-Forum` (R 4.6.1 with
> shiny, bslib, dplyr, stringr, HLAtools 1.11.0, immunogenetr 1.5.0
> installed; use Rscript).
>
> The full task is in PROMPT.md in that repository. Read it and follow it
> exactly. Build from the prompt only; do not read the first-attempt/ folder.

## 2. The agent's report

> I built `app/app.R` from PROMPT.md; it is the only new file, nothing is
> committed, and I did not read `first-attempt/`. All three checks in the
> prompt passed.
>
> **What it does**
> - **Data:** the version is pinned in `IPD_IMGT_HLA_VERSION <- "3.65.0"` and
>   shown in the sidebar. The alignment is loaded once at startup with
>   `library(HLAtools)` and `alignmentFull(...)[["prot"]][["DPB1"]]`.
> - **Epitope function:** `dpb1_epitopes(allele)` is a plain function outside
>   the server. It returns a list with `found`, `message`, `typed`, `shown`,
>   `label` (e.g. "DPB1\*04:01 (shown: DPB1\*04:01:01:01)"), `is_null`, an
>   `epitopes` table (sequence and status: common, uncommon or not known) and a
>   `residues` table for the nine positions.
> - **Allele names, all through immunogenetr:**
>   - Cleaning: `HLA_validate()` removes stray spaces and junk,
>     `HLA_prefix_remove()` strips `HLA-`/`DPB1*`, then
>     `HLA_prefix_add("DPB1*")` puts the prefix back.
>   - Field count: the smallest n for which `HLA_truncate(x, n) == x`.
>   - Matching: alignment names are truncated with
>     `HLA_truncate(..., keep_suffix = TRUE)` and compared as whole names. That
>     is why `04:01` never lands on `DPB1*04:01:01:24N`.
> - **Symbols that are not amino acids:** `*`, `.` and `X` make an epitope "not
>   known". In the residue table they appear as "not sequenced", "no residue"
>   or "stop codon", never as the symbol itself.
> - **Layout:** bslib `page_sidebar()` with a server-side selectize list of all
>   3,170 alleles, defaulting to DPB1\*04:01:01:01. A red alert appears for null
>   alleles. Base pipe only, only the six allowed packages, and "IPD-IMGT/HLA
>   Database" is always written in full.
>
> **Test output** (epitopes 55 / 84 / 96/170)
>
> | Input | Shown | 55 | 84 | 96/170 | Null |
> |---|---|---|---|---|---|
> | DPB1\*01:01:01:01 | same | AAE common | DEAV common | K/I common | no |
> | DPB1\*02:01:02:01 | same | DEE common | GGPM common | R/T common | no |
> | DPB1\*04:01:01:01 | same | AAE common | GGPM common | R/T common | no |
> | DPB1\*04:02:01:01 | same | DEE common | GGPM common | R/T common | no |
> | DPB1\*05:01:01:01 | same | EAE common | DEAV common | K/I common | no |
> | `" 04:01 "` | DPB1\*04:01:01:01 | AAE | GGPM | R/T | no |
> | DPB1\*02:01 | DPB1\*02:01:02:01 | DEE | GGPM | R/T | no |
> | 04:01:01 | DPB1\*04:01:01:01 | AAE | GGPM | R/T | no |
> | HLA-DPB1\*161:01N | DPB1\*161:01N | DED common | not known (stop codon at 85) | not known | yes |
> | 61:01 | "DPB1\*61:01 was not found in IPD-IMGT/HLA Database version 3.65.0." | | | | |
> | 61:01N | DPB1\*61:01N | DED common | not known | not known | yes |
> | 401:01 | not found | | | | |
> | 99:99 | not found | | | | |
> | DPB1\*01:01:04 (96/170 unknown) | same | AAE | DEAV | not known (96 and 170 not sequenced) | no |
> | DQB1\*04:01, abc (extra) | "... is not a valid HLA-DPB1 allele name." | | | | |
>
> **Check across the whole alignment**
> - **Unsuffixed names:** I ran all 4,450 unique two- and three-field names
>   without a suffix. None returned an allele with an N or Q suffix. 388 were
>   "not found", because their only matches are null or Q alleles.
> - **Full names:** every one of the 3,170 full names returns itself.
> - **Null flags:** 179 alleles are flagged as null.
>
> **Launch:** `shiny::testServer()` ran without errors, rendering the epitope
> table and the null alert. `shiny::runApp()` in a background R process
> answered HTTP 200 with the app page.
>
> **Things to review**
> 1. **Null detection:** there is no immunogenetr function for this. The code
>    calls `HLA_truncate(x, 4, keep_suffix = TRUE) ==
>    paste0(HLA_truncate(x, 4, keep_suffix = FALSE), "N")`. That builds a
>    string rather than parsing the name, but it borders on the "no own text
>    handling" rule. Q alleles are correctly not flagged as null.
> 2. **Gene check (not in the prompt):** a typed gene must be DPB1, so
>    `DQB1*04:01` is rejected instead of being quietly read as DPB1\*04:01. A
>    lowercase gene (`dpb1*04:01`) is also rejected, because I did not uppercase
>    the input to avoid hand-written text handling. Typing in the app only
>    searches the dropdown, so this affects console use only.
> 3. **One-field input:** `04` also works and shows the first `DPB1*04` allele;
>    the prompt did not ask for that.
> 4. **Speed:** each lookup truncates all 3,170 names, about 0.1 s per click.
>    That is fine in the app, but the full-alignment check took about 10
>    minutes.
> 5. **Null alert wording:** it says a null allele "is not expressed on the cell
>    surface, so these epitopes are not presented to antibodies." Worth checking
>    that this reads right for the talk.

**What to notice:** the rule "use immunogenetr for allele names" did most of
the work. `HLA_truncate()` already keeps expression suffixes, so comparing
truncated names makes `04:01` skip `DPB1*04:01:01:24N` and makes `401:01`
return "not found" instead of a null allele, with no special-case code. The
agent also went beyond the prompt in one sensible way (rejecting
`DQB1*04:01` rather than silently reading it as DPB1) and was honest about
the one place it still wrote its own text handling (item 1). That item is
wrong, though: immunogenetr *does* have a function for it. An agent is not an
expert in your tools; you are.

## 3. Our message

> Thanks, this is good. Rejecting other genes is right, and one-field input is
> fine; keep both. One change:
>
> Null detection: immunogenetr does have a function for this.
> GLstring_drop_non_expressed() removes non-expressed alleles, and for a single
> allele it returns NA when the allele is null (see
> ?GLstring_drop_non_expressed; its `suffixes` argument lets you choose which
> suffixes count as non-expressed). Replace the string-building in
> is_null_allele() with it, keeping Q alleles not flagged as null. Keep
> everything else as it is. Re-run the null checks (61:01N, 161:01N,
> 04:01:01:24N, a Q allele such as 1622:01Q, and 04:01:01:01) and the
> whole-alignment null count (you reported 179), confirm the app still
> launches, and report back. You don't need to redo the 10-minute
> full-alignment lookup check unless the change touches the lookup.

## 4. The agent's report

> `is_null_allele()` now uses `GLstring_drop_non_expressed()` instead of
> building a string, and the null results are the same as before: 179 alleles
> flagged, Q alleles not flagged. The app still launches without errors.
>
> The new function is `is.na(GLstring_drop_non_expressed(allele, suffixes =
> "N"))`, with its comments updated for the audience. I passed `suffixes = "N"`
> because the function's default also drops S and C alleles, and the app's note
> is specifically about null alleles. The DPB1 alignment has no S or C alleles,
> so the result is the same either way. Nothing else in `app/app.R` changed.
>
> **Null checks:**
>
> | Input | Allele shown | Null |
> |---|---|---|
> | 61:01N | DPB1\*61:01N | yes |
> | 161:01N | DPB1\*161:01N | yes |
> | 04:01:01:24N | DPB1\*04:01:01:24N | yes |
> | 1622:01Q | DPB1\*1622:01Q | no |
> | 04:01:01:01 | DPB1\*04:01:01:01 | no |
>
> **Across the whole alignment:** 201 alleles carry a suffix; 179 are flagged as
> null and the other 22 are all Q alleles.
>
> **Launch:** under `shiny::testServer()` the app flags 161:01N as null and
> 1622:01Q as not null. Started with `shiny::runApp()` in a background R
> process, it served the app page (HTTP 200) with no errors in the log.
>
> I didn't redo the 10-minute check across all two- and three-field names,
> because the lookup code is unchanged. Nothing is committed, and the test
> scripts are in the scratchpad, not the repo.

## 5. Independent check

Before the app was accepted, a separate agent that had not seen this
conversation reviewed the code and re-checked the epitopes against the
alignment. It found nothing blocking:

- All 13 sample alleles (including EAE, null, Q, and unknown-96/170
  alleles) matched an independent lookup against the raw alignment.
- Short names behave as intended: `" 04:01 "`, `04:01`, `DPB1*04:01` and
  `HLA-DPB1*04:01` all show DPB1\*04:01:01:01; `61:01`, `401:01` and `99:99`
  are not found; `61:01N` and `401:01:01N` return the null allele, flagged.
  None of 300 random unsuffixed two-field names returned a suffixed allele.
- Every allele name is handled with immunogenetr; no hand-written text
  matching is left.
- Every number in this log was recomputed and matches.
- The app starts in about 6 seconds and uses about 300 to 400 MB, well within
  Posit Cloud's limit.

It also raised one wording problem: the null-allele alert says the epitopes
"are not presented to antibodies", while the table below still labels a null
allele's 55 epitope "common". And it suggested a few ways to make the code
easier for beginners to read.

## Why there was a third attempt

Nick's review of this app: it should not truncate anything, because
immunogenetr's `GLstring_regex()` is the right tool for finding an allele from
a name with any number of fields; and a null allele should show no epitopes
at all, instead of listing them under a warning that they are not expressed.
Both went into the prompt, and a fresh agent built the app again. See
[`../PROMPT.md`](../PROMPT.md) and [`../agent-log.md`](../agent-log.md).
