# Getting Started with Informatics and Coding in the Clinical HLA Laboratory

Code for the Laboratory Professionals Forum (hosted by the Technologists'
Affairs Committee) at the ASHI 52nd Annual Meeting, Montréal, Wednesday
October 28, 2026.

Presenter: Nicholas Brown, PhD, F(ACHI)

No coding experience is needed. You can follow along on your laptop through
Posit Cloud, or just watch; everything shown is in this repository.

## Following along on Posit Cloud

1. Create a free account at <https://posit.cloud/>.
2. Open the Posit Cloud project link shared for the session. The packages are
   already installed there.
3. Open the files below in order and run the gray code chunks with the green
   "play" arrow.

To run the code on your own computer instead, install R and RStudio, open
`project.Rproj`, and run `00 - Install packages.R` once.

## What's here

| File | What it covers |
|---|---|
| `01 - Intro to coding with R.qmd` | Just enough R: objects, vectors, functions, the pipe, and data tables |
| `02 - immunogenetr.qmd` | The [immunogenetr](https://cran.r-project.org/package=immunogenetr) package: GL strings, mismatches, ranking donors by match grade |
| `03 - HLAtools.qmd` | The [HLAtools](https://cran.r-project.org/package=HLAtools) package: searching the IPD-IMGT/HLA Database protein alignments, using HLA-DPB1 as the example |
| `app/app.R` | A Shiny app that shows the three major antibody epitopes (55, 84, 96/170) of any HLA-DPB1 allele |
| `PROMPT.md` | The prompt used to have an AI coding agent build the app |
| `agent-log.md` | The conversation with the agent, from prompt to working app |
| `first-attempt/`, `second-attempt/` | The earlier prompts, their logs and the apps they produced, kept for comparison |

## Running the app

Open `app/app.R` and click **Run App** at the top of the editor, or run:

```r
shiny::runApp("app")
```

The app downloads the HLA-DPB1 protein alignment from the IPD-IMGT/HLA Database
(version 3.65.0) when it starts, so it needs an internet connection.

The epitope definitions follow Morris et al. (PMID
[39277973](https://pubmed.ncbi.nlm.nih.gov/39277973/)). For surrogate
single-antigen-bead alleles that share these epitopes, see the
[DPB1-surrogate-report](https://github.com/k96nb01/DPB1-surrogate-report)
repository.

## Building the app with an AI coding agent

The app was not written by hand. `PROMPT.md` is the prompt, and
`agent-log.md` records what the agent reported back, what we asked it to
change, and how its work was checked. The point is not that the AI writes
perfect code; it is that a clear description of what you want (the data, the
rules, and how to test it) gets you a working first draft in minutes, which
you then read, test, and own.

## License

MIT; see [`LICENSE`](LICENSE). You are welcome to copy and adapt this code for your own lab.
