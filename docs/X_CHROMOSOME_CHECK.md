# Checking the centromere of the *Capra hircus* X chromosome

In T2T-goat2.0 the only satellite block of the X sits about 7.2 Mb from the end, while in the autosomes the satellite starts at about 30 kb. I wanted to know whether that was a detection, orientation or assembly error. It is not. What I found fits the cytogenetics of the X and what Wu et al. (2024) describe. The numbers are in `data/x_chromosome_check.csv` and in section 2.7 of notebook 01.

## What the literature says

Most descriptions agree that the goat X is acrocentric. Hansen (2009) describes a Q-band karyotype with acrocentric or telocentric autosomes, an acrocentric X and a submetacentric Y. Kowal-Mierzwa et al. (2026) find the goat X acrocentric, while the cattle X is submetacentric, and Wu et al. (2024) report that all the autosomes and the X are acrocentric. The exception is Putri et al. (2018), who classify the X of the Peranakan Etawa breed as subtelocentric, with a centromeric index of 19 to 33 % measured under the light microscope. Their classification is inconsistent between sexes, so I treat it as an outlier.

The X is one of the largest chromosomes. Wu et al. describe a pseudoautosomal region (PAR) of about 7 Mb on the short arms of the X and the Y, and say that the X has the shortest centromere (0.6 Mb), with SatII coinciding with CENP-A peaks on CHI22, CHI28 and the X. One caution: Kowal-Mierzwa et al. found by FISH a 3-Mb inversion between the bovine and the goat X, and two assembly errors in the goat genome. Their abstract does not say in which assembly.

## What I found in the sequence

The only satellite block of the X goes from 7.233 to 8.085 Mb. It is 853 kb long and holds 638 kb of satellite (592 kb of SatII and 47 kb of SatI), so it starts at 5.0 % of the chromosome length. That gives a centromeric index of about 5 % and an arm ratio (q/p) of about 19, which is an acrocentric chromosome and does not match the 19 to 33 % of Putri et al. There is no other satellite block of 30 kb or more in the 144.7 Mb of the chromosome. When I run SRF on the X alone, it finds a 2801-bp unit (about 4 × 702 bp, a higher-order repeat of SatII) and some microsatellites, so there is no X-specific satellite that I am missing.

The first 7.04 Mb of the X align with 7.02 Mb of the Y, on the same strand and with an identity of about 0.89, as Wu et al. describe for the PAR. That short arm is ordinary euchromatin: it has 62 genes in 7 Mb (8.9 genes per Mb) against 8.5 genes per Mb over the whole X. Finally, the X of haplotype M lines up with the X of T2T-goat2.0 over 144.68 of 144.69 Mb, on the same strand.

## How I read it

The short arm of the X, about 7.2 Mb, holds the pseudoautosomal region, which is gene-rich and homologous to the Y. The centromere comes right after it, made of almost pure SatII with about 0.64 Mb of satellite, and the rest is the long arm. That is the structure of an acrocentric chromosome with a short arm of about 5 % of its length. In the autosomes the satellite block starts at the end simply because they have no euchromatic short arm.

I found no sign of a detection or orientation error. What I cannot rule out with these data is a local assembly error like the ones Kowal-Mierzwa et al. describe, since I did not check in which assembly they found theirs.

## What it means for the analysis

I keep the value of the X. It stays a highly influential point in H1 and H2 (Cook's distance 1.19) because it combines a very short centromere with a long chromosome. That is biology and not an error, so I report the results with and without the X.

## References
- Hansen KB. Q-band karyotype of the goat (*Capra hircus*) and the relation between goat and bovine Q-bands. *Hereditas* (2009).
- Kowal-Mierzwa W, et al. Identification of cryptic evolutionary rearrangements between cattle (Bovinae) and goat (Caprinae) X chromosomes. *Anim Genet* 57(2):e70094 (2026). https://doi.org/10.1002/age.70094
- Putri A, et al. Chromosome characteristic of Peranakan Etawa (PE) goat (*Capra hircus* Linn.) as Indonesian local breed. *IOP Conf Ser Earth Environ Sci* (2018).
- Schibler L, Di Meo GP, Cribiu EP, Iannuzzi L. Molecular cytogenetics and comparative mapping in goats (*Capra hircus*, 2n = 60). *Cytogenet Genome Res* 126:77–85 (2009). https://doi.org/10.1159/000245908
- Wu H, et al. Telomere-to-telomere genome assembly of a male goat reveals variants associated with cashmere traits. *Nat Commun* 15:10041 (2024). https://doi.org/10.1038/s41467-024-54188-z
