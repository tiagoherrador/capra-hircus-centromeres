# Verification of the centromere of the *Capra hircus* X chromosome

**Question.** In T2T-goat2.0 the only satellite block of the X is ~7.2 Mb from the end, while in the autosomes it starts at ~30 kb. Is this a detection, orientation or assembly error?

**Conclusion.** No. The result is consistent with the cytogenetics of the X and with what Wu et al. (2024) describe. Details are in `data/x_chromosome_check.csv` and in section 2.7 of the notebook.

## 1. What the literature documents

| Aspect | Evidence |
|---|---|
| Morphology of the X | **Acrocentric** in most descriptions: Hansen 2009 (Q-band karyotype: autosomes acrocentric or telocentric, X acrocentric, Y submetacentric); Kowal-Mierzwa et al. 2026 (goat: X acrocentric; cattle: submetacentric); Wu et al. 2024 (all autosomes and the X are acrocentric). |
| Exception | Putri et al. 2018 (Peranakan Etawa) classify the X as subtelocentric with a centromeric index of 19–33 % (light microscopy; classification inconsistent between sexes). It is an outlier relative to the rest of the literature. |
| Arms | The X is one of the largest chromosomes. Wu et al. describe a **pseudoautosomal region (PAR) of ~7 Mb on the short arms of the X and the Y**. |
| Centromere sequence | Wu et al.: the X centromere is the shortest (0.6 Mb) and **SatII co-occurs with CENP-A peaks on CHI22, CHI28 and the X**. |
| Caution | Kowal-Mierzwa et al. 2026 detected, by FISH, a ~3-Mb inversion between the bovine and caprine X and **two assembly errors in the goat genome** (the abstract does not say in which assembly). |

## 2. What was found in the sequence (T2T-goat2.0)

| Check | Result |
|---|---|
| Position of the only satellite block of the X | 7.233–8.085 Mb (853 kb long; 638 kb of satellite: 592 kb SatII and 47 kb SatI). **5.0 %** of the length from the end. |
| Centromeric index (short arm / total) | ≈ 5 %; arm ratio q/p ≈ 19 → **acrocentric**. It does not match the 19–33 % of Putri et al. |
| Other satellite blocks on the X | None ≥ 30 kb in 144.7 Mb. |
| SRF applied to the X alone | Finds a 2801-bp unit (≈ 4 × 702 bp: a SatII higher-order repeat) and microsatellites; **there is no X-specific satellite that we are missing**. |
| Short arm = PAR | The first 7.04 Mb of the X align with 7.02 Mb of the Y (+ strand, identity ≈ 89 %), as Wu et al. describe. |
| Genes of the short arm | 62 genes in the first 7 Mb (8.9 genes/Mb) versus 8.5 genes/Mb over the whole X: ordinary euchromatin. |
| Orientation | The X of haplotype M is collinear with the X of T2T-goat2.0 over 144.68 of 144.69 Mb, same strand. |

## 3. Interpretation

The short arm of the X (≈ 7.2 Mb) contains the pseudoautosomal region, gene-rich and homologous to the Y. The centromere (almost pure SatII, ≈ 0.64 Mb of satellite) comes right after it, and the rest is the long arm. This is the structure of an acrocentric chromosome with a short arm of ≈ 5 % of its length. In the autosomes, in contrast, the satellite block starts at the end because they have no euchromatic short arm.

There is no evidence of a detection or orientation error. What cannot be ruled out with these data is a local assembly error of the class described by Kowal-Mierzwa et al. 2026 (it was not checked in which assembly).

## 4. Consequence for the analysis

The value for the X is **kept**. It remains a highly influential point in H1 and H2 (Cook's distance 1.19) because it combines a very short centromere with a long chromosome; that is biology, not an error, and results are reported with and without the X.

## References
- Hansen KB. Q-band karyotype of the goat (*Capra hircus*) and the relation between goat and bovine Q-bands. *Hereditas* (2009).
- Kowal-Mierzwa W, et al. Identification of cryptic evolutionary rearrangements between cattle (Bovinae) and goat (Caprinae) X chromosomes. *Anim Genet* 57(2):e70094 (2026). https://doi.org/10.1002/age.70094
- Putri A, et al. Chromosome characteristic of Peranakan Etawa (PE) goat (*Capra hircus* Linn.) as Indonesian local breed. *IOP Conf Ser Earth Environ Sci* (2018).
- Schibler L, Di Meo GP, Cribiu EP, Iannuzzi L. Molecular cytogenetics and comparative mapping in goats (*Capra hircus*, 2n = 60). *Cytogenet Genome Res* 126:77–85 (2009). https://doi.org/10.1159/000245908
- Wu H, et al. Telomere-to-telomere genome assembly of a male goat reveals variants associated with cashmere traits. *Nat Commun* 15:10041 (2024). https://doi.org/10.1038/s41467-024-54188-z
