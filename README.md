## Sex-Specific Differences in MicroRNA Expression During Human Fetal Lung Development

**Key Links:** \
Manuscript: [Frontiers in Genetics](https://www.frontiersin.org/articles/10.3389/fgene.2022.762834) (open access) \
GEO: [GSE200153](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200153) (raw data/.fastqs)

**Summary:** We profiled whole prenatal human lung samples using microRNA-sequencing, demonstrating differences in autosomal microRNA (i) average expression levels and (ii) expression patterns across gestational age between males and females--even in very early lung development. The sex differences identified include microRNAs thought to regulate prenatal lung development and post-natal disease susceptibility.

## Repo Files

This repo includes analysis scripts and key processed data files.

In addition to the count matrices provided above, see above GEO link for acces to the raw miRNA-seq data.

```
File Structure
├─ PublicData
│   └─ miRNA_counts.tsv            # filtered counts, 898 miRNAs x 298 samples
│   └─ sample_metadata.tsv         # 298 samples x 6 metadata variables
└─ Scripts
│   └─ 01 through 03               # miRNA-seq .fastq -> counts
│   └─ 11 through 18               # general dataset QC/EDA, count modeling
│   └─ 21 through 26               # sex dimorphism analyses for this manuscript
```

## Citation

Lin NW, Liu C, Yang IV, Maier LA, DeMeo DL, Wood C, Ye S, Cruse MH, Smith VL, Vyhlidal CA, Kechris K. Sex-Specific Differences in MicroRNA Expression During Human Fetal Lung Development. Frontiers in Genetics. 2022;13.

PMID: [35480332](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9037032/) \
doi: [10.3389/fgene.2022.762834](https://doi.org/10.3389/fgene.2022.762834)
