# Centromeric heterochromatin of *Capra hircus*
### A sequence-based case study with public data

How big are goat centromeres? Do they scale with chromosome size, what are they made of, and how much do the two haplotypes of one animal differ? This project tries to answer those questions from the telomere-to-telomere goat genome (Wu et al. 2024), using only public data and going from the raw download to the statistics. It started with one animal and was then extended: 20 more public assemblies were screened, and a second animal with well-resolved centromeres was added to repeat the tests.

> Personal portfolio and practice project. It is not peer reviewed and it is not a research article. Read the [limitations](#limitations) before interpreting any number. The analysis was developed with the help of an AI coding assistant (Claude Code); the decisions and the results were reviewed by the author.

![Summary figure](data/figures/summary.png)

## At a glance

| | |
|---|---|
| Data | T2T-goat2.0 genome and its two haplotypes (2.0P, 2.0M), from a male Inner Mongolia Cashmere goat (`IMCG_204112`); PacBio reads of five other breeds; 20 other public assemblies, of which a Kiko × Saanen F1 female (two haplotypes) is used as a second animal |
| Approach | Orient chromosomes, find satellite units (SRF), annotate, delimit each centromere, run the statistics, then screen other assemblies and repeat the tests |
| Notebooks | `01` one animal, `02` screening of 20 assemblies, `03` two animals |
| Stack | Python (pandas, statsmodels, SciPy), minimap2, samtools, KMC, SRF |
| Runtime | notebook 01 about 20 min (about 5 GB of downloads), notebook 02 about 2 h (about 15 GB downloaded in total, deleted as it goes), notebook 03 seconds; on an Apple-silicon laptop |
| Main finding | The centromeric satellite blocks are large, differ a lot between haplotypes, and scale only weakly with chromosome size (a slope of about 0.3 in autosomes, not seen when the X is included) |

## Questions

- H1, scaling. Does centromere length scale with chromosome size?
- H2, composition. Does the satellite composition (SatI, SatII, SatIII) explain part of the variation?
- H3, homologs. Do the two homologs differ more in the centromere than in the rest of the chromosome?
- Exploratory: does the satellite content vary between breeds?

## Key results

Every number comes from `data/results_summary.json` and the tables in `data/`. The first four items are from the first animal alone (notebook 01), the next ones from the extension (notebooks 02 and 03).

- Centromeres. The centromeric blocks of T2T-goat2.0 span 0.85 to 10.4 Mb, and the shortest is the X. SatI adds up to about 117 Mb and SatII to about 26 Mb, close to the totals in Wu et al. (118.2 and 26.5 Mb). SatIII (about 6 Mb) is under-detected.
- H1, no detectable scaling. The slope of ln(centromere) on ln(remaining chromosome length) is −0.12 (95 % CI −0.80 to 0.57), and 0.19 (−0.10 to 0.47) without the X. With 80 % power the analysis detects a slope of at least 0.39 without the X, or 0.66 with it. So a strong scaling is unlikely, but a weak or moderate one cannot be excluded. This fits the literature on the functional centromere, which is fairly uniform within a species.
- H2, weak support. No satellite family scales with chromosome size. The SatI/SatII log-ratio adds 0.73 of R² to the scaling model, but only 0.25 without the X, and it is partly circular because it describes the same block.
- H3, supported. In 18 autosomes with a complete centromere in both haplotypes, the median absolute difference is 2.7 Mb in the centromere and 0.11 Mb in the rest of the chromosome (permutation p ≈ 5×10⁻⁵). SatI drives most of it. The result holds on relative and absolute scales, without the outlier CHI21, with a strict gap criterion, and without the three blocks of haplotype M that have very little satellite (16 chromosomes, p ≈ 2×10⁻⁵).
- Breeds, exploratory. Compared with the reference animal, the SatI/SatII ratio is lower in Tibetan (−1.12, CI −1.32 to −0.92) and in Boer (−0.26, CI −0.43 to −0.09), and not different in the other three. The Tibetan library looks atypical (see below), so a technical cause is likely.

- Screening of 20 public assemblies (02). Only the two haplotypes of a Kiko × Saanen F1 female have centromeres resolved about as well as T2T-goat2.0: 29 and 27 of 30 have a block with no gap, enough satellite, and between half and twice the satellite of the same chromosome in T2T-goat2.0. The best of the others resolve 17, 15 and 13. Most of them have between 35 and 205 Mb of satellite sitting in unplaced scaffolds, so their centromeres are fragmented rather than absent.
- H1 with two animals (03). With the X, no haplotype shows scaling (slopes between −0.06 and 0.19). In autosomes the slope is about 0.3 in both animals (0.30 for IMCG, 0.30 for Kiko × Saanen, 0.30 pooled, 95 % CI 0.10 to 0.49), and the two animals do not differ from each other. This is inside the interval of the first animal alone (0.19, −0.10 to 0.47), which did not have the power to see it. So the reading of H1 is a weak scaling of about 0.3 in autosomes, flattened by the X.
- H3 with two animals (03). It replicates in the second animal with a smaller effect: 1.17 Mb in the centromere against 0.18 Mb in the rest (permutation p = 0.03 on the absolute scale, below 0.0001 on the relative one). Over the 44 pairs of both animals it is clear (1.28 against 0.14 Mb, p = 0.0001). It depends on three chromosomes with large expansions of SatI, which are discussed below.
- Within an animal against between animals (03). In 15 autosomes usable in all four haplotypes, the median difference in centromere length is larger between the two haplotypes of IMCG (0.38 in |ln ratio|) than between KikoP and SaanenM (0.17) or between animals (0.21). With 15 chromosomes this is a hint and not a result (p = 0.16).

## How the analysis went

The original design was to delimit the centromeres by combining CENP-A ChIP-seq, satellite density and sequence complexity, then compare chromosomes, haplotypes and animals. What follows is what happened when the design met the public data. Several steps changed the plan.

### 1. Getting the data
The three assemblies come from NCBI with MD5 checks, and the chromosomes are renamed to `CHI1…CHI29, CHIX, CHIY`. In the maternal haplotype the X is called `30`, so it is matched by length. The PacBio reads of the other breeds come from ENA. They are streamed and cut at 0.5 Gb per run, because the full files are 50 to 120 GB each. The first version ran on Google Colab. It ran out of memory in SRF and timed out on a base-level alignment, so everything moved to a laptop, with a different annotation strategy (see step 3).

### 2. Orienting the chromosomes
The goat chromosomes are acrocentric, so each centromere should be near one end. I first aligned the SatI reference sequence to both ends to see which one holds it. That reference is short and diverged, and for most chromosomes it hits nothing at either end, so the test alone says very little. A real check came later, after SRF: with the full set of satellite units, the centromere of every chromosome in the three assemblies (89 in total) sits at the left end of the deposited sequence, as the genome description says. The check is part of the notebook (section 2.4).

### 3. Finding and classifying the satellites
SRF counts frequent 151-mers and assembles them into repeat units. It found 21 units. They are grouped by sequence into SatI (816-bp monomer), SatII (702 bp) and SatIII (a 22-bp motif). A first version classified units by length only, and it counted units of 132, 770 and 2115 bp as SatII even though they are mostly repeats of the SatIII motif. Classifying by sequence brought the SatII total from about 33 Mb to about 26 Mb, in line with the paper. The units are then mapped back to each chromosome with minimap2, without base-level alignment, and each 100-bp bin is assigned to its best hit.

### 4. What is a centromere here?
It is the satellite block with the most sequence in the proximal window, after merging pieces closer than 500 kb. This is an operational definition. It marks the satellite domain and not the CENP-A domain, so it can be longer than published values. CENP-A ChIP-seq was tried on a subsample of 2 million read pairs and showed no enrichment in any satellite family. Doing it per chromosome would need the full runs (about 29 GB) and unique mapping, because short reads inside nearly identical arrays cannot be assigned to a chromosome. It was dropped for that cost, not because it is impossible ([details](docs/NOTE_chipseq_removed.Rmd)).

![Composition of each centromeric block](data/figures/composition_by_chromosome.png)

Each block is mostly SatI, with SatII in every chromosome and SatIII in 24 of 30. Between about 15 and 55 % of each block is not satellite, since it includes the stretches between merged pieces.

### 5. The odd case of the X chromosome
The X has the shortest centromere, and its only satellite block lies 7.2 Mb from the end, while in the autosomes it starts at about 30 kb. I checked whether this was an error ([details](docs/X_CHROMOSOME_CHECK.md)). The first 7 Mb of the X match the Y (the pseudoautosomal region described by Wu et al.), they are gene-rich, SRF run on the X alone finds no new satellite, and the X of haplotype M is collinear with it. The goat X is acrocentric with a short arm of about 5 % of its length. It is not an error, but it is an influential point in H1 and H2.

![Satellite along the first 15 Mb of the X](data/figures/x_chromosome.png)

### 6. Reworking the statistics
The first round of statistics had four problems that were fixed before writing anything up:
- H1 asked whether the slope differs from 1, which is almost trivial when the slope is near 0. The main contrast is now slope = 0, with intervals, influence diagnostics and a power analysis, so "no scaling" can be told apart from "no way of seeing it".
- H2 was partly circular, because the proportions describe the same block whose length is being explained. The main analysis is now the allometry of each family, and composition enters as a single SatI/SatII log-ratio.
- H3 compared relative differences, which favours the centromere simply because it is a short segment. It now also compares absolute differences in Mb.
- The breed comparison had no uncertainty. It now uses a bootstrap over reads, with quality control of each run.

### 7. Haplotype M: a problem found late
When the repository was reviewed, three centromeres of haplotype M (CHI7, CHI9 and CHI29) turned out to have very little satellite, between 10 and 21 % of the block, against at least 44 % in every other case.

![Satellite density of each block](data/figures/satellite_density_by_block.png)

In CHI7 and CHI29 the chromosome of M is shorter than in P and in T2T by roughly the missing amount, so sequence seems to be absent from the assembly. In CHI9 the sequence is there, but about 4 Mb of it does not match the SRF units, which were built from the collapsed assembly only. Because a density threshold is not a biological definition of a centromere, these blocks are not removed. They get a flag, `low_satellite`, meaning that their length cannot be measured reliably here, and H3 and the haplotype replicate of H1 are reported with and without them. H3 holds either way.

![The two haplotypes, centromere by centromere](data/figures/haplotypes_by_chromosome.png)

### 8. Breeds
The reads are mapped to the SatI and SatII units, and each breed is compared with the reference animal.

![Satellite ratio by breed](data/figures/breeds.png)

The Tibetan run stands out for two reasons. About 12.5 % of its reads carry SatI or SatII, against 5 to 7 % in the others, and its read N50 is longer. Yet its SatI and SatII add up to only about 2.5 % of the bases, so a good part of its satellite signal may be a family this analysis does not measure. For that reason the Tibetan difference is not read as a change in the SatI/SatII ratio.

### 9. Looking for more genomes
There are 70 *Capra hircus* assemblies in NCBI. Twenty were chosen for the screening: chromosome level (one scaffold level) and mostly a high contig N50, plus the RefSeq reference ARS1.2. The T2T-goat versions of the original animal were left out. The assembly of the Inner Mongolia Agricultural University has a different submitter and is assumed to be another animal, which was not verified. Every one was processed with the same method as the first animal, one at a time, deleting the big files after each. The result is the table of `data/screening/screen_summary.csv`.

![Screening of 20 assemblies](data/figures/screen_overview.png)

A high contig N50 does not mean that the centromeres are resolved. Satellite arrays are often collapsed or cut in long-read assemblies. The yardstick is T2T-goat2.0, so a block counts as resolved only if it is complete and its satellite is between half and twice that of the same chromosome there. A real expansion beyond twice would be missed.

### 10. A second animal
The Kiko × Saanen F1 female has two haplotypes of two breeds, which gives a second pair of homologs. Four haplotypes of two animals are not enough to estimate variation between individuals, so that part is descriptive. The tests of H1 and H3 were repeated, using only centromeres that are complete and not flagged `low_satellite`. These criteria do not use T2T-goat2.0 as a yardstick, so they do not bias the comparison between animals.

![Two animals](data/figures/two_animals.png)

### 11. Four centromeres reviewed by hand
Four centromeres of the second animal have more than twice the satellite of the same chromosome in T2T-goat2.0 or have a gap: Kiko CHI14 and Saanen CHI5, CHI7 and CHI15. In CHI14, CHI7 and CHI15 the chromosome is longer by about the same amount as the centromere, so the extra sequence is in the assembly and is not a detection error. It may be a real expansion of the SatI array or a duplication made by the assembler, and long reads from the same animal would be needed to tell them apart. In Saanen CHI15 about 3 Mb of the block is spacer between two small satellite islands that were merged. Saanen CHI5 reaches the end of the 15-Mb window, so its block is probably cut, and it has a gap; it is left out. These cases stay in the analysis, and they matter for reading H3: without CHI7, CHI14 and CHI15 the replication in the second animal drops to p = 0.17 on the absolute scale, while the pooled result (p = 0.0005) and the slope of H1 (0.28) hardly move. Removing them only because their differences are large would be selecting against the hypothesis, so they are kept and this is stated as a caveat.

## Limitations

- Two animals. The first analysis used one individual, and the extension adds a second one. The two haplotypes of Kiko × Saanen come from two breeds, so their difference mixes individual and breed. Two animals are not enough to estimate variation between individuals. The differences between haplotypes may also include assembly artefacts.
- The screening uses T2T-goat2.0 as its yardstick, and its SatI-III units come from that animal only. A satellite that the T2T animal does not have would look like missing sequence in another assembly. The SRF check on each assembly found only a few units that do not resemble the T2T ones, and they were not characterised.
- Long reads of the second animal were not used, so a large array cannot be told apart from an assembly duplication.
- The block of the centromere is searched in a 15-Mb window. A block that reaches the edge, like Saanen CHI5, is cut.
- Operational centromere. It is the satellite domain, not the functional CENP-A domain, and includes some non-satellite sequence. It could not be validated with CENP-A data. H3 uses the block length, and the check with satellite only is pending.
- Haplotype M. Three blocks look incomplete, and the cause of CHI9 is not resolved.
- SatIII, a 22-bp motif, is poorly detected by alignment.
- Breeds. PacBio CLR subreads at about 0.17× coverage, genome-wide fractions only, and SatIII is not measured.
- Power. With n = 29 or 30 chromosomes, a weak scaling cannot be excluded.
- The X is a highly influential point, so results are reported with and without it.

## Reproducibility

Three levels, from cheapest to most complete:

1. Read the executed notebooks in `notebooks/`.
2. Re-run the statistics in seconds from the tables committed in `data/`, with no downloads: section 3 of notebook 01 and the whole of notebook 03. The breeds step is skipped if `data/reads/` is absent.
3. Re-run everything from the raw public data. Notebook 01 downloads about 5 GB, takes about 20 min and leaves about 13 GB of intermediate files in `work/`. Notebook 02 needs 01 first, takes about 2 h and downloads about 15 GB in total, deleting each assembly when it finishes. Notebook 03 needs the tables of 01 and 02.

```bash
# macOS (Apple silicon). External tools via Homebrew:
brew install minimap2 samtools

python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m ipykernel install --user --name capra --display-name "Capra hircus"

jupyter lab notebooks/01_one_animal.ipynb      # run from the repository root or from notebooks/
# or, headless:  python scripts/run_notebook.py 0 notebooks/01_one_animal.ipynb
# notebook 02 can be tested with one assembly:  CAPRA_ONLY=KikoP python scripts/run_notebook.py 0 notebooks/02_screening_20_assemblies.ipynb
```

KMC, SRF and k8 are downloaded and built into `tools/` automatically (macOS arm64 only). The notebook skips downloads whose MD5 already matches.

## Repository structure

```
├── notebooks/
│   ├── 01_one_animal.ipynb                  download, satellites, centromeres and statistics for IMCG_204112
│   ├── 02_screening_20_assemblies.ipynb     the same method on 20 public assemblies
│   └── 03_two_animals.ipynb                 H1, H3 and within/between comparison with a second animal
├── data/                         tables produced by the pipeline (raw downloads and reads are git-ignored)
│   ├── screening/                summary and per-assembly tables of notebook 02
│   ├── two_animals/              tables of notebook 03
│   └── figures/                  summary and per-topic figures
├── docs/
│   ├── X_CHROMOSOME_CHECK.md     evidence on the position of the X centromere
│   ├── NOTE_chipseq_removed.Rmd  why CENP-A ChIP-seq was dropped
│   └── REFERENCES.md
├── scripts/run_notebook.py       headless, cell-by-cell runner with a log
└── requirements.txt  CITATION.cff  LICENSE  LICENSE-CC-BY-4.0.txt
```

The state of the repository before the extension is kept as the tag `v1.0-one-animal`.

## License

- Code: MIT ([LICENSE](LICENSE)).
- Text, figures and documentation: CC BY 4.0 ([LICENSE-CC-BY-4.0.txt](LICENSE-CC-BY-4.0.txt)).
- The sequence data are public (NCBI, ENA) and are not redistributed here. The tables in `data/` are derived results.

## Key references

Wu et al. 2024 (*Nat Commun* 15:10041) for the genome; Zhang et al. 2023 (*Genome Res*) for SRF; Plačková et al. 2021 and Zhang and Dawe 2012 on centromere size and genome size; Logsdon et al. 2024 and 2025 (*Nature*) on variation of centromeric satellite arrays; Greenacre et al. 2022 (*Stat Sci*) for pairwise log-ratios. The full list is in [docs/REFERENCES.md](docs/REFERENCES.md).
