# Annot-Compare User Guide

## Overview

Annot-Compare compares two GTF genome annotation files and reports differences at every level of the annotation hierarchy: genes, transcripts, and transcript features (exons, CDS, UTRs, introns, start/stop codons). Each difference is written as a row in a TSV output file tagged with an impact level, so results can be filtered and prioritized for the research question at hand.

Typical use cases:
- Tracking what changed between two annotation releases before updating a pipeline
- Deciding which analyses need to be rerun because the underlying annotation changed
- Checking whether genes or transcripts of interest changed in ways that affect prior findings
- Comparing annotations from different sources (e.g., Ensembl vs. RefSeq)

---

## Requirements

- Java 21 or later
- Two GTF files to compare (target and query)
- A genome FASTA file indexed with `samtools faidx` (produces a `.fai` file)

---

## Getting the JAR

Build from source with Maven:

```bash
mvn package -DskipTests
```

The JAR is produced at `target/annot-compare-*.jar`.

---

## Quick Start

```bash
java -jar annot-compare.jar analysis \
    --target-gtf target.gtf \
    --query-gtf query.gtf \
    --fasta genome.fa \
    --fidx genome.fa.fai \
    --config config.json \
    --o results.tsv
```

The **target** is the reference annotation (e.g., the version currently in use). The **query** is the annotation being compared against it (e.g., a new release). Both must cover the same genome assembly, with contigs listed in the same order.

---

## Input Files

### GTF Files

Standard GTF format. The tool reads gene and transcript entries and their `gene_id`, `transcript_id`, `gene_biotype`, and `transcript_biotype` attributes. Both files must use the same contig names and be ordered consistently.

### FASTA and Index

Index the genome FASTA before running:

```bash
samtools faidx genome.fa
```

This creates `genome.fa.fai`. Pass both files via `--fasta` and `--fidx`. The index must correspond to the same FASTA file.

---

## Configuration

All comparison behaviour is controlled by a JSON file passed via `--config`. See `example/config.json` for a ready-to-use template.

### Top-Level Structure

```json
{
  "features": { ... },
  "transcript_features": { ... },
  "gene_biotypes": { ... }
}
```

---

### `features` — Comparator Settings

Controls individual comparators that check specific properties. Each comparator supports some combination of `enabled`, `threshold`, and `impact`.

| Key | What it checks | `threshold` supported |
|---|---|---|
| `transcript_length` | Transcript total length | yes |
| `transcript_start` | Transcript start coordinate | yes |
| `transcript_stop` | Transcript stop coordinate | yes |
| `transcript_strand` | Strand change | no |
| `transcript_biotype` | Biotype change | no |
| `start` | Feature start coordinate | yes |
| `stop` | Feature stop coordinate | yes |
| `length` | Feature length | yes |
| `sequence` | DNA sequence change | no |
| `same_protein` | Protein sequence change (CDS only) | no |

- `enabled` (boolean): whether this comparator runs at all.
- `threshold` (integer, bp): differences within this many base pairs are not reported as changes.
- `impact` (string): override the default impact level for this comparator. Valid values: `HIGH`, `MODERATE`, `LOW`.

Example:

```json
"features": {
  "transcript_start": { "enabled": true, "threshold": 0, "impact": "HIGH" },
  "sequence": { "enabled": true },
  "same_protein": { "enabled": true }
}
```

---

### `transcript_features` — Feature Types to Analyze

Selects which GTF feature types are included in the feature-level comparison.

| Key | Description |
|---|---|
| `transcript` | Transcript-level boundaries |
| `exon` | Exon regions |
| `intron` | Intron regions (derived from exon pairs) |
| `cds` | Coding sequence regions |
| `5UTR` | Five-prime UTR |
| `3UTR` | Three-prime UTR |
| `utr` | UTR (unspecified direction) |
| `start_codon` | Start codon position |
| `stop_codon` | Stop codon position |

Example:

```json
"transcript_features": {
  "exon": { "enabled": true },
  "cds": { "enabled": true },
  "intron": { "enabled": true },
  "5UTR": { "enabled": false },
  "3UTR": { "enabled": false }
}
```

---

### `gene_biotypes` — Biotype Filter

Restricts comparison to transcript pairs where at least one side has a listed biotype.

```json
"gene_biotypes": {
  "allowed": ["protein_coding"]
}
```

An empty list (`[]`) disables the filter and includes all biotypes.

The biotype values must match what appears in the `transcript_biotype` attribute of your GTF files exactly. An unrecognized value silently filters out all transcripts.

---

## CLI Reference

| Flag | Required | Description |
|---|---|---|
| `--target-gtf` | Yes | Path to the target (reference) GTF file |
| `--query-gtf` | Yes | Path to the query GTF file |
| `--fasta` | Yes | Path to the genome FASTA file |
| `--fidx` | Yes | Path to the FASTA index file (`.fai`) |
| `--config` | Yes | Path to the JSON configuration file |
| `--o` | Yes | Path for the output TSV file |
| `--map-with-strand` | No | Whether strand is considered when mapping transcripts (`true`/`false`, default: `true`) |

---

## Output

### Main Output File (`--o`)

A tab-separated file with a header row followed by one row per detected difference.

| Column | Description |
|---|---|
| `impact` | Severity: `HIGH`, `MODERATE`, `LOW`, or `NONE` |
| `contig` | Chromosome / contig name |
| `targetGeneId` | Gene ID in the target annotation |
| `queryGeneId` | Gene ID in the query annotation |
| `targetBioType` | Gene biotype in the target |
| `queryBiotype` | Gene biotype in the query |
| `featureType` | Level of the difference: `gene`, `transcript`, `exon`, `cds`, `intron`, `5UTR`, `3UTR`, `start_codon`, `stop_codon` |
| `difference` | Type of difference (see tables below) |
| `targetTranscriptId` | Transcript ID in the target (empty for gene-level rows) |
| `queryTranscriptId` | Transcript ID in the query (empty for gene-level rows) |
| `targetTranscriptBiotype` | Transcript biotype in the target |
| `queryTranscriptBiotype` | Transcript biotype in the query |
| `targetFeatureStart` | Feature start coordinate in the target |
| `queryFeatureStart` | Feature start coordinate in the query |
| `targetFeatureStop` | Feature stop coordinate in the target |
| `queryFeatureStop` | Feature stop coordinate in the query |
| `targetStrand` | Strand in the target (`+` or `-`) |
| `queryStrand` | Strand in the query (`+` or `-`) |

---

### `difference` Values

**Gene level** (`featureType` = `gene`)

| Value | Meaning |
|---|---|
| `missingInQueryFile` | Gene present in target but absent in query |
| `missingInTargetFile` | Gene present in query but absent in target |
| `start` | Gene start coordinate changed |
| `stop` | Gene stop coordinate changed |
| `strand` | Strand changed |
| `seq` | Sequence changed |
| `length` | Gene length changed |
| `contig` | Gene is on a different contig |

**Transcript level** (`featureType` = `transcript`)

| Value | Meaning |
|---|---|
| `none` | Transcripts are identical (informational row) |
| `start` | Transcript start coordinate changed |
| `stop` | Transcript stop coordinate changed |
| `seq` | Transcript sequence changed |
| `length` | Transcript length changed |
| `strand` | Strand changed |
| `biotype` | Biotype changed |
| `missingInTarget` | Transcript present in query but absent in the matched target gene |
| `missingInFileQuery` | Transcript present in target but absent in the matched query gene |

**Feature level** (`featureType` = `exon`, `cds`, `intron`, etc.)

| Value | Meaning |
|---|---|
| `featureTypeMissingInTargetTranscript` | This feature type is entirely absent in the target transcript |
| `featureTypeMissingInQueryTranscript` | This feature type is entirely absent in the query transcript |
| `missingFeatureEntryFileInTarget` | A specific feature instance is absent in the target |
| `missingFeatureEntryFileInQuery` | A specific feature instance is absent in the query |
| `start` | Feature start coordinate changed |
| `stop` | Feature stop coordinate changed |
| `length` | Feature length changed |
| `seq` | Feature DNA sequence changed |
| `protein` | CDS protein sequence changed |

---

### Supplementary Output Files

Three additional files are written alongside the main output.

| File | Content |
|---|---|
| `<output>.minimap2` | Same format as the main output, but for transcript pairs found via minimap2 sequence alignment rather than coordinate overlap |
| `<output>.unmapped_queries` | Transcripts from the query that could not be matched to any target transcript (tab-separated: contig, transcript_id) |
| `<output>.unmapped_targets` | Transcripts from the target that could not be matched to any query transcript |

Transcripts in the unmapped files represent structural additions (unmapped queries) or removals (unmapped targets) with no sequence-similar counterpart in the other annotation.

---

## Example

A complete example is provided in the `example/` directory. It compares Ensembl releases 113 and 114 on a subset of the human genome.

To run it:

1. Download the full human genome FASTA from the Ensembl FTP (link in the README) and place the unzipped file in `example/input/`.
2. Run `make` from inside the `example/` directory.

The configuration in `example/config.json` can be customized before running.

---

## Troubleshooting

**The program exits with a contig mismatch error.**
Both GTF files must list contigs in the same order. Re-sort the files to match if they differ.

**The output has very few or no rows despite expecting changes.**
Check that the biotype values in `gene_biotypes.allowed` exactly match the `transcript_biotype` attribute in your GTF files. A mismatch silently filters out all transcripts. You can temporarily set `allowed` to `[]` to disable the filter and confirm the tool is running.

**Parsing fails with a FASTA index error.**
The `.fai` file must be generated by `samtools faidx` from the exact FASTA file passed to `--fasta`. Regenerate it if there is any mismatch.

**Many transcripts appear in the unmapped files.**
This is expected when comparing annotations across substantially different sources or versions. The `.minimap2` output file may recover some of these pairs via sequence alignment.
