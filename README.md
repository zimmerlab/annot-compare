# Annot-Compare

![Logo](assets/logo.png)

![GitHub release](https://img.shields.io/github/v/release/zimmerlab/annot-compare)
[![Release](https://github.com/zimmerlab/annot-compare/actions/workflows/release.yml/badge.svg)](https://github.com/zimmerlab/annot-compare/actions/workflows/release.yml)
![Java](https://img.shields.io/badge/Java-21-blue)
![License](https://img.shields.io/github/license/zimmerlab/annot-compare)

Genome annotation changes between releases and across sources, those changes carry real consequences for downstream analyses. When a lab reruns an RNA-seq or variant analysis under a new Ensembl release, read counts shift because exon boundaries moved, variant consequences change because CDS boundaries shifted, and transcripts disappear or are reassigned. Most of the time, this goes undetected: there is no systematic way to know which results were affected, or how badly.

Annot-Compare makes that risk visible. Given two GTF annotation files, two Ensembl releases, or two sources such as Ensembl and GENCODE: it determines which genes and transcripts correspond to each other across the two annotations and reports every difference at the gene, transcript, and feature level (exons, CDS, UTRs, introns, start/stop codons). Each difference is tagged with an impact level so results can be filtered to what actually matters for a given analysis.

**Typical uses:**

- Determining which results need to be rerun after an annotation update, rather than rerunning everything or nothing
- Checking whether a gene or transcript of interest changed in a way that affects a prior finding before publication
- Assessing whether two studies using different annotation versions or sources are directly comparable

**What makes this non-trivial:** annotation identifiers are not stable across sources, and are not reliably stable even between consecutive Ensembl releases. A diff based on shared IDs misses the majority of real correspondences. Annot-Compare matches genes and transcripts by CDS structure and coordinate overlap, identity by function rather than by name, and falls back to minimap2 sequence alignment for transcripts that cannot be matched by position alone. This makes it applicable across sources, not only between consecutive releases of the same annotation.

---

## Requirements

- Java 21 or later
- Two GTF files to compare
- A genome FASTA file for each annotation, indexed with `samtools faidx`

---

## Quick Start

```bash
java -jar annot-compare.jar newMapping \
    --target-gtf <path-to-target-gtf> \
    --target-fasta <path-to-target-fasta> \
    --target-fai <path-to-target-fai> \
    --query-gtf <path-to-query-gtf> \
    --query-fasta <path-to-query-fasta> \
    --query-fai <path-to-query-fai> \
    --output <path-to-output-file>
```

The **target** is the reference annotation (e.g., the version currently in use). The **query** is the annotation being compared against it (e.g., a new release).

For full documentation of inputs and output format, see the [User Guide](USER_GUIDE.md).

---

## Example

A ready-to-run example comparing Ensembl releases 113 and 114 on a subset of the human genome is provided in the [`example/`](example) directory. The genome FASTA must be downloaded separately:

```
https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
```

Unzip it into `example/input/`, then run `make` from inside the `example/` directory.
