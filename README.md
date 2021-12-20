This repo stores analysis code and processed data accompanying "Sex-Specific Differences in MicroRNA Expression during Human Fetal Lung Development" (research article; under review).

We profiled whole prenatal human lung samples using microRNA-sequencing and found differences in microRNA (i) average expression levels and (ii) expression patterns across gestational age between males and females, even in very early lung development.

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

In addition to the processed data provided here, we are currently also working on making all raw data files accessible via GEO/SRA.

This repo will be updated with a finalized abstract, links to the published manuscript, & .fastq accession information as it becomes available.

