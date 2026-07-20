# cSVD protein candidates tested for heart-failure shared-locus evidence

This repository contains the minimum code and aggregate derived results needed
to inspect and reuse the principal Mendelian-randomization (MR) and
colocalization definitions and to redraw source-independent publication displays for the
manuscript **“From Mendelian randomization signals to shared-locus
qualification: testing cerebral small-vessel disease protein candidates in
heart failure.”**

The release is intentionally narrow. It is not a mirror of provider summary
statistics and does not contain participant-level data, restricted data,
institutional data, local caches, internal review materials or every audit
artifact generated during the study.

## Repository contents

- `scripts/run_stage_g_hf_mr.py`: allele harmonization, instrument eligibility,
  Wald-ratio or multiplicative-random-effects IVW estimation, hierarchical
  Simes aggregation and protein-family multiplicity control.
- `scripts/run_shared_locus_test.R`: ABF/SuSiE routing, prespecified priors,
  PP3/PP4 criteria and explicit handling of the no-credible-set-pair state.
- `scripts/reproduce_main_figures.py`: redraws Figures 1, 2 and 4, plus the
  qualification-state and prior-sensitivity components of Figure 3, from the
  released aggregate results without accessing provider source files.
- `data/publication_results.json`: aggregate 51-protein screening states, three
  assay-specific regional evaluations, the Apo E2 prior-sensitivity grid and
  the observed 55→56→51→48→3→0 evidence flow, plus the 112-row assay-precision
  display used for Figure 4.
- `LICENSE`: BSD 3-Clause License for code.
- `DATA_LICENSE.md`: CC BY 4.0 terms for the aggregate derived result file.

## Data availability

The aggregate values needed to inspect the reported evidence states and redraw
Figures 1, 2 and 4 and the qualification/prior components of Figure 3 are openly available in `data/publication_results.json` under
CC BY 4.0. They contain no participant-level records, personal identifiers,
credentials, local paths or full provider summary-statistic archives.

Source GWAS and pQTL summary statistics remain available from their original
providers under provider-specific access and redistribution terms. The study
used public resources from the cSVD proteogenomics study, WUSTL pQTL resources,
deCODE, HERMES, BioBank Japan and reference resources described in the
manuscript. Users wishing to rerun MR or colocalization from source summary
statistics must obtain those files directly from the providers and comply with
their terms. This repository does not redistribute HERMES, deCODE, WUSTL or
BioBank Japan source archives, complete variant-level caches or LD-reference
files. UK Biobank and controlled-access individual-level resources were not
used in this release.

The exact regional exposure-pQTL and HERMES-HF association panels in main-text
Figure 3 require thousands of provider summary-statistic rows. Those rows are
not included here because this minimal release does not expand the providers'
redistribution scope. Readers may reconstruct those panels after obtaining the
official source files under the original terms.

## Installation

Figure reproduction requires Python 3.12 with `numpy`, `pandas` and
`matplotlib`. Source-level MR additionally requires `scipy` and `pyarrow`.
Colocalization requires R with `coloc` 5.2.3, `data.table`, `jsonlite` and the
dependencies used by `coloc::runsusie`.

## Method self-check

The MR script includes a source-independent self-check:

```bash
python scripts/run_stage_g_hf_mr.py --self-check
```

The colocalization script exposes its command-line contract with:

```bash
Rscript scripts/run_shared_locus_test.R --help
```

Source-level execution requires separately obtained provider inputs in the
relative paths documented by the scripts. No source data are downloaded by
these commands.

## Redraw the main figures

```bash
python scripts/reproduce_main_figures.py
```

By default, vector PDF and 300-dpi PNG files are written to `output/`. Set
`OUTPUT_DIR` to use another output directory.

## Interpretation caveats

- The released results concern one frozen, externally derived family of 51
  cSVD protein candidates; they do not show that complex models, MR signals or
  shared heart–brain biology are generally absent.
- Three HERMES protein-family FDR signals opened three assay-specific regional
  evaluations across two genomic regions. None met the prespecified
  shared-locus qualification criteria.
- Two SuSiE evaluations formed no credible-set pair. This is
  resolution-limited evidence, not proof that a shared causal signal is absent.
- The ABF evaluation was sensitive to the shared-causal prior. The higher-prior
  result was not used to rescue the primary decision.
- BioBank Japan was a post hoc cross-ancestry robustness analysis and could not
  reopen the primary route or trigger new colocalization.
- Downstream expression, spatial, proteomic and perturbation analyses were not
  opened after the genetic gate failed.

## Licenses

Code is released under the BSD 3-Clause License. The aggregate derived results
in `data/publication_results.json` are released under CC BY 4.0. Third-party
source data remain governed by their original terms and are not relicensed by
this repository.

## Citation

Citation details will be added after journal publication.

> The authors. From Mendelian randomization signals to shared-locus
> qualification: testing cerebral small-vessel disease protein candidates in
> heart failure. Manuscript under review.
