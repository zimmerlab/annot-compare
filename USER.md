# Annot-Compare User Guide

## Overview

Annot-Compare compares two GTF genome annotation files and reports differences at every level of the annotation hierarchy: genes, transcripts, and transcript features (exons, CDS, UTRs, introns, start/stop codons). It first establishes which genes and transcripts in one annotation correspond to which in the other, then characterises every structural difference between matched pairs.

The tool provides four modes. `newMapping` and `newTranscriptMapping` produce a mapping file recording which genes and transcripts correspond across the two annotations. `newMappingVal` validates a mapping by checking sequence identity. `filterGeneMappings` filters a mapping file to retain only pairs matched by specific criteria.

---

## Requirements

- Java 21 or later
- Two GTF files to compare
- A genome FASTA file for each annotation (uncompressed or BGZF-compressed/bgzipped)
- A genome FASTA index (`.fai`) for each FASTA file, generated with `samtools faidx`

---

## Getting the JAR

Build from source with Maven:

```bash
mvn package -DskipTests
```

The JAR is produced at `target/annot-compare-*.jar`.

---

## Modes

### `newMapping`

Produces a gene and transcript mapping between two annotations. For each contig, all gene pairs that share at least one structurally similar transcript are identified and written to a mapping file together with the origins that support the match (shared gene ID, shared transcript ID, shared gene name, coordinate overlap, or genomic distance).

```bash
java -jar annot-compare.jar newMapping \
    --target-gtf <path> \
    --target-fasta <path> \
    --query-gtf <path> \
    --query-fasta <path> \
    --output <path> \
    [--target-fai <path>] \
    [--query-fai <path>] \
    [--allowed-types <comma-separated-feature-types>] \
    [--useHomology <true|false>]
```

| Flag | Required | Description |
|---|---|---|
| `--target-gtf` | Yes | Path to the target GTF file |
| `--target-fasta` | Yes | Path to the target genome FASTA (can be `.fa` or `.fa.gz` if bgzipped) |
| `--target-fai` | No | Path to the target FASTA index (defaults to `<target-fasta>.fai`) |
| `--query-gtf` | Yes | Path to the query GTF file |
| `--query-fasta` | Yes | Path to the query genome FASTA (can be `.fa` or `.fa.gz` if bgzipped) |
| `--query-fai` | No | Path to the query FASTA index (defaults to `<query-fasta>.fai`) |
| `--output` | Yes | Path for the output mapping file |
| `--allowed-types` | No | Comma-separated list of feature types to use for structure-based matching (e.g. `exon,CDS`). Defaults to all types. |
| `--useHomology` | No | Whether to use sequence homology in addition to structure for matching (`true`/`false`, default: `true`) |

**Note on FASTA and FAI:**
The program supports both plain-text FASTAs and BGZF-compressed FASTAs (typically produced by `bgzip`). In both cases, a `.fai` index is required. If the `--target-fai` or `--query-fai` flags are omitted, the program will look for a file with the `.fai` extension appended to the FASTA path (e.g., `genome.fa` -> `genome.fa.fai` or `genome.fa.gz` -> `genome.fa.gz.fai`).

Two output files are written:

- `<output>` - the mapping file (see Output Format below)
- `<output>.unmapped` - genes that could not be matched to any gene in the other annotation (columns: `contig`, `geneId`, `name`, `origin`)

---

### `newTranscriptMapping`

Produces a transcript-level mapping using a pre-computed gene-level mapping file. Genes are not re-matched; instead, transcripts are matched within each already-established gene pair. Use this when you have a gene mapping from a previous `newMapping` run and want to refine or re-run transcript matching independently.

```bash
java -jar annot-compare.jar newTranscriptMapping \
    --target-gtf <path> \
    --target-fasta <path> \
    --target-fai <path> \
    --query-gtf <path> \
    --query-fasta <path> \
    --query-fai <path> \
    --gene-mapping <path> \
    --output <path> \
    [--allowed-types <comma-separated-feature-types>] \
    [--useHomology <true|false>]
```

| Flag | Required | Description |
|---|---|---|
| `--target-gtf` | Yes | Path to the target GTF file |
| `--target-fasta` | Yes | Path to the target genome FASTA |
| `--target-fai` | No | Path to the target FASTA index (defaults to `<target-fasta>.fai`) |
| `--query-gtf` | Yes | Path to the query GTF file |
| `--query-fasta` | Yes | Path to the query genome FASTA |
| `--query-fai` | No | Path to the query FASTA index (defaults to `<query-fasta>.fai`) |
| `--gene-mapping` | Yes | Path to an existing gene-level mapping file (output of `newMapping`) |
| `--output` | Yes | Path for the output mapping file |
| `--allowed-types` | No | Feature types to use for structure matching (default: all) |
| `--useHomology` | No | Whether to use sequence homology (`true`/`false`, default: `true`) |

Output files are the same format as `newMapping`.

---

### `newMappingVal`

Validates a mapping file by checking whether the genomic sequences of each mapped transcript pair are identical. Produces a TSV listing every pair with a boolean `sequenceIdentical` column. Useful for assessing mapping quality before downstream analysis.

```bash
java -jar annot-compare.jar newMappingVal \
    --target-gtf <path> \
    --target-fasta <path> \
    --target-fai <path> \
    --query-gtf <path> \
    --query-fasta <path> \
    --query-fai <path> \
    --mapping <path> \
    --output <path> \
    [--allowed-types <comma-separated-feature-types>]
```

| Flag | Required | Description |
|---|---|---|
| `--target-gtf` | Yes | Path to the target GTF file |
| `--target-fasta` | Yes | Path to the target genome FASTA |
| `--target-fai` | No | Path to the target FASTA index (defaults to `<target-fasta>.fai`) |
| `--query-gtf` | Yes | Path to the query GTF file |
| `--query-fasta` | Yes | Path to the query genome FASTA |
| `--query-fai` | No | Path to the query FASTA index (defaults to `<query-fasta>.fai`) |
| `--mapping` | Yes | Path to the mapping file to validate |
| `--output` | Yes | Path for the output validation file |
| `--allowed-types` | No | Feature types to include (default: all) |

Output columns: `targetGeneId`, `queryGeneId`, `targetTranscriptId`, `queryTranscriptId`, `sequenceIdentical` (`true` or `false`).

---

### `filterGeneMappings`

Filters a mapping file to retain only pairs supported by specified mapping origins. Use this to restrict a mapping to high-confidence pairs (for example, only pairs where genes share an ID or a name) before running downstream analysis.

```bash
java -jar annot-compare.jar filterGeneMappings \
    --mapping-file <path> \
    --output-file <path> \
    --allowed-mappings <comma-separated-origins> \
    [--require-all <true|false>]
```

| Flag | Required | Description |
|---|---|---|
| `--mapping-file` | Yes | Path to the mapping file to filter |
| `--output-file` | Yes | Path for the filtered output file |
| `--allowed-mappings` | Yes | Comma-separated list of origins to allow (see table below) |
| `--require-all` | No | If `true`, a pair must have all listed origins to be retained. If `false`, any one is sufficient. Default: `true`. |

Allowed origin values for `--allowed-mappings`:

| Value | Meaning |
|---|---|
| `gene-id` | The two genes share the same gene identifier |
| `gene-name` | The two genes share the same `gene_name` attribute |
| `transcript-id` | The matched transcript pair shares the same transcript identifier |
| `overlapping` | The two genes overlap on the same contig (genomic distance = 0) |
| `distance` | The two genes were matched by genomic proximity |

---

## Output Format

The mapping modes (`newMapping` and `newTranscriptMapping`) produce three output files:

### 1. Main Mapping File (`<output>`)

A tab-separated file recording every transcript pair that supported a gene-level match.

| Column | Description |
|---|---|
| `contig` | Chromosome or contig name |
| `queryId` | Gene ID in the query annotation |
| `targetId` | Gene ID in the target annotation |
| `queryTranscriptId` | Transcript ID in the query that supported the match |
| `targetTranscriptId` | Transcript ID in the target that supported the match |
| `mapping_origins` | Comma-separated list of origins (see below) followed by `,DISTANCE:<bp>` |

**Mapping Origins:**
- `GENE_ID_MAPPING`: Genes share the same identifier.
- `TRANSCRIPT_ID_MAPPING`: At least one transcript pair shares the same transcript identifier.
- `NAME_MAPPING`: Genes share the same `gene_name`.
- `STRUCTURE_BASED_MAPPING`: At least one transcript pair with identical structure cigars (based on feature types and lengths).
- `OVERLAPPING`: Genes overlap on the same contig (distance = 0).
- `DISTANCE`: Genes were matched by proximity (distance > 0).

The `DISTANCE:<bp>` suffix at the end of the `mapping_origins` column indicates the genomic distance between the two genes in base pairs.

### 2. Unmapped Genes File (`<output>.unmapped`)

Records genes that could not be matched to any gene in the other annotation.

| Column | Description |
|---|---|
| `contig` | Chromosome or contig name |
| `geneId` | Identifier of the unmapped gene |
| `name` | Name of the unmapped gene (if available) |
| `origin` | Either `QUERY` (gene exists only in query) or `TARGET` (gene exists only in target) |

### 3. Execution Time (`<output>.overall_time`)

A single-line file containing the total execution time of the mapping process in milliseconds.

---

## Troubleshooting

**The program exits with a contig mismatch error.**
Both GTF files must list contigs in the same order. Re-sort the files to match if they differ.

**The mapping file is empty or has very few entries.**
Check that at least one annotation contains `protein_coding` genes. The structure-based matching requires at least one protein-coding gene in each pair. If both annotations use non-standard biotype values, no pairs will be found.

**Parsing fails with a FASTA index error.**
Each `.fai` file must be generated by `samtools faidx` from the exact FASTA file passed alongside it. Regenerate the index if there is any mismatch.
