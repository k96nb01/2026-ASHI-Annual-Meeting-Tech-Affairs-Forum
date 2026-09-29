# Building the DPB1 epitope app with an AI coding agent

*First attempt, kept for comparison. The app it produced is [`app.R`](app.R) in this folder; the current prompt is [`../PROMPT.md`](../PROMPT.md). Links below refer to the original layout.*

This is the record of how [`app/app.R`](app/app.R) was built. The agent was a
Claude Code sub-agent (Claude Opus), started and supervised by a second Claude
Code session working with Nick. The agent read the prompt, wrote the code, and
ran it in R to test it. What follows is every message we sent it and every
report it sent back. The reports are copied verbatim; the agent's individual
steps between messages (reading files, running R) are not shown.

Total time from prompt to a tested app: about 3 minutes for the first version
and 3 more minutes for the change.

---

## 1. Our message

> Repository: `2026-ASHI-Annual-Meeting-Tech-Affairs-Forum` (R 4.6.1 with
> shiny, bslib, dplyr, stringr, HLAtools 1.11.0, immunogenetr 1.5.0
> installed; use Rscript).
>
> The full task is in PROMPT.md in that repository. Read it and follow it
> exactly.

The whole specification is in [`PROMPT.md`](PROMPT.md). Writing it took
longer than the agent took to build the app, and that is the right way round.

## 2. The agent's report

> I built the app as a single file, `app/app.R`. The epitope function gives the
> expected results for all the required alleles, and the app starts headlessly
> without errors. I created no other files and did not commit anything.
>
> **What it does**
> - **Data:** the database version is pinned in `IPD_IMGT_HLA_VERSION <- "3.65.0"`
>   at the top of the file and shown in the sidebar. The alignment is downloaded
>   once when the app starts.
> - **Epitope logic:** it is all in one plain function, `dpb1_epitopes(allele)`.
>   It uses immunogenetr's `HLA_prefix_remove()` and
>   `HLA_prefix_add(prefix = "DPB1*")`, so `04:01:01:01`, `DPB1*04:01:01:01` and
>   `HLA-DPB1*04:01:01:01` all work. It returns a list with the allele name, a
>   null flag, the epitopes table (sequence plus common, uncommon or not known)
>   and the per-position residues.
> - **Layout:** built with `page_sidebar()`. The allele list is searched on the
>   server and holds all 3,170 alleles, with DPB1\*04:01:01:01 selected at start.
>   Null alleles get a warning box. The residue table shows `?` for any position
>   that isn't a real amino acid.
> - **Code rules:** it uses `|>` throughout, is commented in plain language,
>   only uses the six allowed packages, and names the "IPD-IMGT/HLA Database"
>   in full everywhere.
>
> **Test output** (epitopes 55 / 84 / 96/170)
>
> | Allele | 55 | 84 | 96/170 | Null |
> |---|---|---|---|---|
> | DPB1\*01:01:01:01 | AAE common | DEAV common | K/I common | no |
> | DPB1\*02:01:02:01 | DEE common | GGPM common | R/T common | no |
> | DPB1\*04:01:01:01 | AAE common | GGPM common | R/T common | no |
> | DPB1\*04:02:01:01 | DEE common | GGPM common | R/T common | no |
> | DPB1\*01:01:04 (96/170 not known: positions 96 and 170 are `*`) | AAE common | DEAV common | not known | no |
> | DPB1\*61:01N | DED common | not known | not known | yes |
> | DPB1\*161:01N (85 is `X`, then `.`) | DED common | not known | not known | yes |
>
> Across all 3,170 alleles:
>
> | Epitope | Common | Uncommon | Not known |
> |---|---|---|---|
> | 55 | 2781 | 348 | 41 |
> | 84 | 3029 | 54 | 87 |
> | 96/170 | 2531 | 20 | 619 |
>
> An unknown name such as `DPB1*99:99:99` stops with a clear message.
>
> **Launch check:** `shiny::testServer()` gave the correct epitope table for
> 04:01:01:01 and the null warning for 61:01N. `shiny::runApp()` in a
> background process started and served the page (HTTP 200) with the title
> and version 3.65.0. The log had no errors.
>
> **Things to check**
> 1. **`HLAtools` must be loaded with `library(HLAtools)`.** Calling only
>    `HLAtools::alignmentFull()`, as the prompt's example does, fails with
>    "object 'HLAgazetteer' not found". The app loads it properly and has a
>    comment explaining why.
> 2. **I treated `X` as "not known" as well as `*` and `.`.** The prompt only
>    named `*` and `.`, but the alignment also uses `X`, which marks a stop
>    codon. It appears only in null alleles, but it is not an amino acid, so
>    showing it as one would break the rule. Please confirm.
> 3. **For null alleles, "not known" is not quite the right word.** The `.`
>    there usually means the protein has already stopped, so the residue is
>    absent rather than unknown. The label follows the prompt's rule, and the
>    null-allele warning covers the difference.
> 4. **The function needs a full name that exists in the alignment.** Two-field
>    names like `04:01` give the "not found" error. In the app this can't happen
>    because the search box offers only real names.
> 5. **Two of the required alleles look identical here.** DPB1\*02:01:02:01 and
>    DPB1\*04:02:01:01 have the same residues at all nine positions
>    (DEE / GGPM / R/T). That matches the data, and may be worth pointing out in
>    the talk.

**What to notice:** the agent found two things the prompt got wrong or left
out. The example call in the prompt (`HLAtools::alignmentFull()`) does not
work on its own, and the alignment uses a third non-amino-acid symbol (`X`,
a stop codon) that the prompt never mentioned. It handled both sensibly
and, more importantly, *told us* instead of hiding them. Always read the
report, not just the code.

## 3. Our message

> Thanks, this looks good. Treating "X" as not a real amino acid is right;
> keep it. One change, please:
>
> In the talk I will call the function from the console, and people naturally
> type two-field names. Make dpb1_epitopes() also accept a two-field (or
> three-field) name such as "04:01" or "DPB1\*04:01": if the name isn't an exact
> match, use the first allele in the alignment that starts with that name
> followed by ":" (so "04:01" must not match "04:010"), and include in the
> result which full allele was used, so the app/console shows e.g.
> "DPB1\*04:01 (shown: DPB1\*04:01:01:01)". Keep the clear error for names that
> match nothing. Keep everything else as it is, re-run your checks (add
> "04:01", "DPB1\*02:01" and "99:99" to them), confirm the app still launches,
> and report back.

## 4. The agent's report

> `dpb1_epitopes()` now accepts two- and three-field names and says which full
> allele it used. All checks pass and the app still launches. Only `app/app.R`
> changed, and nothing is committed.
>
> **Behavior**
> - An exact match is tried first.
> - If there is none, it uses the first allele in the alignment that starts with
>   the typed name plus ":". The match is literal text, so "04:01" cannot pick up
>   "04:010".
> - If that also finds nothing, you get the same clear "Allele not found in the
>   IPD-IMGT/HLA Database 3.65.0 alignment: ..." error as before.
> - The result has three new fields: `allele_used` (the full allele whose
>   sequence is shown), `label`, and `is_null`, which now checks `allele_used`.
> - `label` is just the name for an exact match, or e.g.
>   "DPB1\*04:01 (shown: DPB1\*04:01:01:01)" otherwise. The app's header now
>   shows `label`.
>
> **Checks** (epitopes 55 / 84 / 96/170)
>
> | Typed | Label | Result |
> |---|---|---|
> | DPB1\*01:01:01:01 | DPB1\*01:01:01:01 | AAE / DEAV / K/I, all common (unchanged) |
> | DPB1\*02:01:02:01 | DPB1\*02:01:02:01 | DEE / GGPM / R/T, all common (unchanged) |
> | DPB1\*04:01:01:01 | DPB1\*04:01:01:01 | AAE / GGPM / R/T, all common (unchanged) |
> | DPB1\*04:02:01:01 | DPB1\*04:02:01:01 | DEE / GGPM / R/T, all common (unchanged) |
> | 04:01 | DPB1\*04:01 (shown: DPB1\*04:01:01:01) | AAE / GGPM / R/T, all common |
> | DPB1\*02:01 | DPB1\*02:01 (shown: DPB1\*02:01:02:01) | DEE / GGPM / R/T, all common |
> | 04:01:01 | DPB1\*04:01:01 (shown: DPB1\*04:01:01:01) | same as 04:01 |
> | DPB1\*04:02 | DPB1\*04:02 (shown: DPB1\*04:02:01:01) | DEE / GGPM / R/T, all common |
> | 01:01:04 | DPB1\*01:01:04 | AAE / DEAV common, 96/170 not known |
> | HLA-DPB1\*161:01N | DPB1\*161:01N | null; DED common, 84 and 96/170 not known |
> | 99:99 | error | "Allele not found ... DPB1\*99:99" |
> | 04:0 | error | "Allele not found ... DPB1\*04:0" |
>
> - **Every two-field name:** all 1,974 two-field names in the data give an
>   allele that really starts with that name plus ":", or an exact match.
> - **The ":" rule:** no real "04:01" vs "04:010" case exists in version
>   3.65.0. To test it, I put a fake "DPB1\*04:010:01" row at the top of a copy
>   of the alignment, and "04:01" still gave DPB1\*04:01:01:01.
> - **App:** `testServer` shows the header "DPB1\*04:01:01:01" and the null
>   warning for 61:01N. `runApp()` in the background served the page (HTTP 200)
>   with no errors in the log.
>
> **Decision for you:** 190 two-field names give the not-found error because the
> only allele with that name ends in a letter straight after the second field,
> such as DPB1\*61:01N, 64:01N or 120:01N. Typing "61:01" is followed by "N",
> not ":", so it doesn't match. That is exactly the rule you gave. Allowing a
> trailing expression letter (N, L, S, Q and so on) as well as ":" would fix it,
> but it would also let a two-field name land silently on a null allele. I left
> it as specified. Tell me if you want it changed, or have attendees type the
> full "61:01N".

**What to notice:** the agent tested an edge case we asked for ("04:010")
even though no such allele exists yet, by making a fake one. It also spotted a
consequence of our rule we had not thought about (null alleles), and asked
instead of guessing. We kept the rule: `61:01N` is a different allele from a
hypothetical `61:01`, and an app that quietly shows a null allele's epitopes
when you asked for an expressed one would be worse than one that says "not
found".

## 5. Independent check

A separate agent that had not seen the conversation above reviewed the code
and re-checked the epitopes against the alignment. The epitopes were all
correct, but the review showed that the note above is too optimistic: text
matching still let `401:01`, `786:01` and `1677:01` land on null alleles
(their only alleles are nulls). The app did flag them as null, but this is
exactly the kind of allele-name handling immunogenetr already does correctly,
and the prompt should have said to use it. Review also showed that EAE, a
common 55 epitope, was missing from the prompt's table.

Rather than patch this app, the prompt was rewritten and a fresh agent built
the app again. See [`../second-attempt/PROMPT.md`](../second-attempt/PROMPT.md)
and [`../second-attempt/agent-log.md`](../second-attempt/agent-log.md).
