# Annot-Compare Developer Guide

## Motivation and Vision

### The Problem

Genome annotation is not static. Projects like Ensembl release multiple annotation versions per year, and independent sources (Ensembl, RefSeq, GENCODE, and others) produce annotations for the same reference genome using different methods and evidence sets. This creates a pervasive problem in genomics research: **comparison inconsistency**.

When two studies use different annotation versions, it is not clear whether their results are directly comparable. For expression analysis, a change in exon boundaries shifts read counts even on identical sequencing data. For variant interpretation, a changed CDS boundary alters the protein consequence. More concretely:

- Is a study performed under annotation version A still valid under version B?
- If a pipeline is rerun with a new annotation, which results actually need to be recomputed and which are unchanged?
- Have genes or transcripts of interest been modified in a way that invalidates a published finding?
- Over many releases, how stable is the annotation of a given gene - is its structure settled or still frequently revised?

Existing tools either compare annotations at a coarse statistical level or rely on stable IDs across versions, an assumption that cannot be made across sources or even across consecutive Ensembl releases.

### The Vision

Annot-Compare approaches the problem from first principles: to compare annotations you must first determine which gene in one annotation corresponds to which gene in the other. This correspondence cannot be assumed from identifiers alone.

**Gene identity by function.** A gene is ultimately defined by what it encodes. If two genes across annotations represent the same biological entity, they should produce at least one similar transcript, and for coding genes, a similar protein. The primary matching criterion is therefore CDS structure similarity: two genes are candidates for comparison if they share at least one transcript with a similar coding sequence.

**Coordinate overlap as a practical filter.** Genomic position is a strong additional signal: the same gene in two annotation versions will almost always occupy the same chromosomal region. Coordinate overlap (via interval trees and bipartite matching) efficiently generates and ranks candidate transcript pairs without exhaustive pairwise comparison. For very short transcripts where overlap alone is less discriminative, this filter prevents a combinatorial explosion in candidates.

**Sequence alignment as a fallback.** A minority of transcripts will not overlap by coordinate, for example, when a gene is relocated in a new assembly or annotation boundaries shift dramatically. For these, minimap2 sequence alignment is used as a fallback to find matches on the basis of sequence similarity alone.

**Hierarchical comparison.** The ultimate goal is transcript-level and feature-level differences: what changed in the exon structure, the CDS, the UTRs. But transcript comparison requires first knowing which transcripts to pair, and that in turn requires knowing which genes correspond to each other. The comparison is therefore structured hierarchically — gene matching → transcript pairing → feature-level diffing, and results are produced at every level so they can be consumed at whichever resolution the research question demands.

---

## Processing pipeline ideas

- different levels of comparison depth are possible like gene / transcript or missing / changed / what changed exactly
- change can be defined differently for various use cases -> no single definition of "change" useful
- genes are defined as an entity that has a protein / transcript as product (NOT by name, id, etc. as those are set by various sources and subject to arbitrary change)
- annotations are always compared pairwise (source -> target)
- every gene (can be filtered, i.e. protein_coding) from source is mapped to the "same" gene entity in target
- this mapping is performed by comparing the transcripts from gene A from source to gene B from target
- if any transcript is considered to be equal, the gene A is mapped to gene B
- transcript equality is defined by a custom feature vector
- the feature vector is a string built by using the feature length in bases and the feature type as char (where the feature to be used is configurable) such as 15C200C12C (15 bases CDS, intron, 200 bases CDS, intron 12 bases CDS)
- for short transcripts (with one feature) this can be ambigous and therefore the actual sequence can also be used for equality (or similarity)
- the outcome of the gene mapping can be: not mapped (in source or target), mapped one2one, mapped one2many, mapped many2one
- one2many or many2one can be if genes are split or merged (which rarely happens but it does)
- these gene mappings are only the first step as such mappings can be wrong
- then for every mapped gene, additional features can be compared such as gene id, gene name, etc
- for annotations with the same underlying genome assembly, the position can also be used
- when the underlying genome assembly changes, the feature vectors should not change but if they do, sequence similarity can be used with minimap2 for this mapping
- after genes are mapped, a similar mapping step is performed for the transcripts of each gene mapping pair to find out which transcripts are missing/added/changed/unchanged

---

## Data on LFE File System

### Ensembl data
Data was downloaded from the public Ensembl FTP server:

- `/mnt/raidbio2/extdata/projekte/annot-compare/ensembl_fastas/`: contains the fasta files (and indices) for the human reference genomes used by Ensembl
- `/mnt/raidbio2/extdata/projekte/annot-compare/ensembl_sorted_gtfs`: contains the human Ensembl genome annotations

(Fabian): some gtf files required manual intervention because of missing data like missing gene entries or inconsistent chr annotations -> not easily downloadable from ensembl ftp (should be kept)

### Project data

`/mnt/raidbio2/extproj/projekte/annotation/annot-compare/`: contains scripts and results from previous work (BA, hiwi, etc.) which is not required in the current project state

### Results for Ensembl comparison

Ensembl annotations are compared consecutively from 46 up to 115

`/mnt/raidbio2/extproj/projekte/annotation/annot-compare/BA_mathis/new_mapping/output/`: the results folder

These comparisons are done on gene and transcript level as well as with sequence similarity (falsely called homology _hom currently) using minimap2

---

## Building

### Internal dependency on zimmerlab/gtf-utils

Private repo requires authentication

#### aquire github personal access token (PAT)

`https://github.com/settings/tokens/` with `read:packages` access

#### set PAT in maven config

`settings.servers.server.id` in `settings.xml` must match `repositories.repository.id` in `pom.xml`

```xml
# ~/.m2/settings.xml
<settings>
	<servers>
		<server>
			<id>github-zimmerlab</id>
			<username>USERNAME</username>
			<password>PAT</password>
		</server>
	</servers>
</settings>
```

```bash
mvn package -DskipTests
```

Output JAR: `target/annot-compare-*.jar`

To run the bundled example:

```bash
cd example && make
```

---

## Run Modes

Modes are selected by passing the mode name as the first argument to the JAR. Spring Boot activates the corresponding `@Profile` and the matching `CommandLineRunner` handles the rest. All modes remain runnable regardless of which ones are shown in the help output.

### User-facing modes

**`newMapping`**
Maps genes and transcripts between two annotations using structure-based similarity (transcript feature cigars), with optional sequence homology. Produces a mapping TSV and an unmapped gene list.
Flags: `--target-gtf`, `--target-fasta`, `--target-fai`, `--query-gtf`, `--query-fasta`, `--query-fai`, `--output`, `--allowed-types` (optional), `--useHomology` (optional, default `true`).

**`newTranscriptMapping`**
Maps transcripts within gene pairs defined by an existing gene mapping file. Same inputs as `newMapping` plus `--gene-mapping`.

**`newMappingVal`**
Validates a mapping file by comparing the genomic sequences of each transcript pair. Reports which pairs are sequence-identical and which are not.
Flags: same as `newMapping` plus `--mapping` (the file to validate).

**`filterGeneMappings`**
Filters a mapping file to rows whose `mapping_origins` column contains the specified origin types.
Flags: `--mapping-file`, `--output-file`, `--allowed-mappings` (comma-separated: `gene-id`, `gene-name`, `transcript-id`, `overlapping`, `distance`), `--require-all` (default `true`).

### Research and utility modes

**`analysis`**
The original annotation comparison pipeline. Maps transcripts by coordinate overlap and bipartite matching, with a minimap2 sequence-alignment fallback, then runs all configured comparators on matched pairs. Produces a hierarchical TSV of differences tagged by impact level.
Flags: `--target-gtf`, `--query-gtf`, `--fasta`, `--fidx`, `--config`, `--o`, `--map-with-strand` (optional), `--mapping-only-cds` (optional).

**`cliqueAnalysis`**
Analyzes the overlap cluster structure of genes across the two annotations without running a full comparison. Reports how many gene pairs are 1:1, 1:many, or many:many matches. Writes a per-entry TSV and a summary `clique_report.csv`.
Flags: `--target-gtf`, `--query-gtf`, `--output`.

**`liftedDifferences`**
Compares a query annotation against a target using a lifted version of the query (e.g. produced by liftoff). Accepts separate FASTAs for each assembly. Isolates annotation-driven differences from assembly-driven differences.
Flags: `--target-gtf`, `--query-gtf`, `--lifted-query-gtf`, `--target-fasta`, `--target-fidx`, `--query-fasta`, `--query-fidx`, `--output`.

**`gtfStats`**
Computes descriptive statistics (gene/transcript/feature counts, biotype distributions, length distributions) for one or more GTF files in a directory. Writes one set of TSVs per input file.
Flags: `--inDir`, `--outDir`, `--recursive` (optional).

**`addMetaFeatures`**
Rewrites a GTF, deriving missing transcript-level start/stop coordinates from child feature records. Use this to normalize GTF files where the transcript row is absent or has empty coordinate fields.
Flags: `--gtf`, `--o`.

**`benchmark`**
Aggregates Java Flight Recorder (JFR) profiling files from performance benchmark runs into summary TSVs.
Flags: `--dir` (directory of subdirectories each containing one `.jfr` file), `--o`, `--pre` (optional filename prefix).

### Internal / development modes

**`firstAnalysis`**
Predecessor to `analysis`. Loads entire GTF files into memory and matches genes only by shared identifier. Output format differs from `analysis`. Superseded and not maintained.

**`seqExtractor`**
Extracts sequences for two hardcoded gene IDs from two GTFs and writes them as a FASTA. A one-off development artifact.

**`test`**
Development scratch runner. Hardcoded to chromosome 3, no meaningful output. Not for general use.

---

## Package Map

```
com.github.zimmerlab.gtfcompare
├── Application.java                  — Spring Boot entry point; routes to runners by mode
├── AnnotComparator.java              — Comparison engine; drives per-pair comparison and result accumulation
│
├── compare/                          — Comparator interfaces and implementations
│   ├── ComparisonFeature             — Interface for feature-level comparators (exon, UTR, intron, …)
│   ├── CDSComparisonFeature          — Interface for CDS-specific comparators
│   ├── TranscriptComparisonFeature   — Interface for transcript-level comparators
│   ├── GeneComparisonFeature         — Interface for gene-level comparators
│   ├── ComparisonContext             — Input bundle passed to every comparator
│   └── comparators/
│       ├── transcriptfeatures/       — Transcript-level: start, stop, length, strand, biotype
│       ├── transcript/               — Feature-level: start, stop, length, sequence
│       └── gene/                     — Gene-level: start, stop, length, strand, contig
│
├── mapping/                          — Transcript and exon matching (for analysis mode)
│   ├── OverlappingTranscripts        — Interval tree overlap + bipartite matching for transcript pairing
│   ├── ExonMapping                   — Needleman-Wunsch gap alignment for exon pairing within a transcript pair
│   ├── Minimap2Bundler               — Extracts the bundled minimap2 binary to a temp directory at runtime
│   └── Minimap2Validator             — Runs minimap2 on unmapped transcripts; returns additional transcript pairs
│
├── newmapping/                       — Structure-based gene and transcript mapping
│   ├── Mapping                       — Main mapping entry point
│   ├── Similarity                    — Structure similarity and sequence homology logic
│   ├── model/                        — Data models for mapping results and origins
│   └── outpututil/                   — Writers for mapping and unmapped output
│
├── newmappingval/                    — Mapping validation logic
│
├── model/                            — Data models (primarily for analysis mode)
│   ├── TranscriptPair                — A matched target/query transcript pair
│   ├── FeaturePair                   — A matched target/query feature instance (exon, CDS, …)
│   ├── GenePair                      — A matched target/query gene
│   ├── Impact                        — Enum: HIGH, MODERATE, LOW, NONE
│   ├── config/                       — Configuration model (ComparisonConfig, ConfigJSON, FeatureConfig)
│   └── comparison/                   — Result types at each hierarchy level
│
├── parser/
│   ├── FidxParser                    — Parses a FASTA index (.fai) for random-access sequence extraction
│   └── EnsemblMappingParser          — Parses legacy Ensembl mapping files
│
├── runner/                           — CLI entry points, one per mode (Spring @Profile)
│   ├── NewMappingRunner              — Entry point for newMapping mode
│   ├── AnalysisRunner                — Entry point for detailed analysis mode
│   ├── NewMappingValidationRunner    — Entry point for mapping validation
│   ├── FilterMappingsRunner          — Entry point for filtering mapping files
│   └── ...                           — Other utility and research runners
│
└── utils/
    ├── GenomeSequenceExtractor       — Random-access FASTA sequence extraction via the .fai index
    ├── ResultWriter                  — Flattens the result hierarchy into sorted TSV rows
    ├── OutputLine                    — A single output row with its sort key
    └── Constants                     — Feature type name constants and comparator name constants
```

---

## Core Processing Pipeline

The main analysis mode (`AnalysisRunner`) processes the two GTF files **one contig at a time** to bound memory usage. For each contig the following steps run in sequence:

1. **Parse** — `GtfFile.parseNextContig()` loads all genes and transcripts for the current contig from both GTF files.

2. **Coordinate-overlap mapping** — `OverlappingTranscripts.map()` builds interval trees over both transcript sets, identifies overlapping pairs, clusters them into loci (gene-level groups), and runs bipartite matching within each locus to produce an optimal set of `TranscriptPair` objects. Transcripts with no match on either side are collected as unmapped.

3. **Minimap2 fallback** — `Minimap2Validator.validateWithMinimap2()` aligns the unmapped query transcripts against the unmapped target transcripts using the bundled minimap2 binary. Matches produce a second set of `TranscriptPair` objects.

4. **Main comparison** — `AnnotComparator.compare()` iterates all transcript pairs from step 2. For each pair:
   - Biotype filtering is applied per the config allowlist.
   - `ExonMapping.pairExonsByGapAlignment()` aligns exons between the two transcripts using Needleman-Wunsch gap alignment to produce `FeaturePair` lists.
   - All enabled `TranscriptComparisonFeature` comparators run via `ServiceLoader`.
   - All enabled `ComparisonFeature` / `CDSComparisonFeature` comparators run on each feature pair.
   - Results accumulate into the `ComparisonResult` hierarchy.

5. **Minimap2 comparison** — The same `AnnotComparator` runs a second time on the minimap2 pairs, writing to a separate `.minimap2` output file.

6. Output — `ResultWriter.writeComparisonResult()` flattens the result hierarchy into `OutputLine` rows, sorts them by type and content, and appends them to the output TSV.

Transcripts still unmatched after both passes are written to the `.unmapped_queries` and `.unmapped_targets` files.

---

## New Mapping Pipeline

The `newMapping` mode implements a structure-based mapping approach that is assembly-independent.

1. **Transcript Cigars** — For each transcript, a "cigar" string is generated based on the sequence of its features (e.g., `100E200I100E` for 100bp exon, 200bp intron, 100bp exon).
2. **Initial Matching** — Genes are matched if they share at least one transcript with an identical cigar string (structure-based).
3. **Identifier Fallback** — Shared gene IDs, transcript IDs, and gene names are also used to establish correspondences.
4. **Sequence Homology** — If enabled, transcripts that don't match by structure or ID are compared using BioJava's Smith-Waterman alignment to find high-similarity matches.
5. **Output** — Results are written to a mapping TSV where each row represents a transcript pair, tagged with the origins of the match.

---

## Plugin Architecture (ServiceLoader)

Comparators are discovered at runtime via Java's `ServiceLoader`. Four interfaces exist:

| Interface | Scope |
|---|---|
| `ComparisonFeature` | Any feature type (exon, UTR, intron, …) |
| `CDSComparisonFeature` | CDS features only |
| `TranscriptComparisonFeature` | Transcript-level properties |
| `GeneComparisonFeature` | Gene-level properties |

Each interface is registered in a corresponding file under:

```
src/main/resources/META-INF/services/
```

Each file lists the fully qualified class name of every implementation, one per line. `ServiceLoader` discovers them at startup without any further wiring.

### Adding a New Comparator

1. Implement the appropriate interface. `getName()` returns the config key that enables this comparator. `compare(ComparisonContext ctx)` returns `true` when a difference is detected.
2. Add the fully qualified class name to the corresponding `META-INF/services/` file.
3. If the comparator needs configuration (threshold, impact override), add the key to `ComparisonConfig` and handle it in `ConfigJSON`.
4. Handle the comparator's name in the relevant `addTo*Comparison()` switch in `AnnotComparator` to propagate the result into the output model.

---

## Key Algorithms

| Algorithm | Class |
|---|---|
| Interval tree construction and overlap queries | `mapping/OverlappingTranscripts.java` |
| Bipartite matching within a locus | `mapping/OverlappingTranscripts.java` (uses jgrapht) |
| Needleman-Wunsch exon gap alignment | `mapping/ExonMapping.java` |
| Minimap2 invocation and result parsing | `mapping/Minimap2Validator.java`, `mapping/Minimap2Bundler.java` |
| Random-access FASTA sequence extraction | `utils/GenomeSequenceExtractor.java`, `parser/FidxParser.java` |

---

## Configuration System

`ConfigJSON` is parsed from the user's JSON file by Jackson and converted to a `ComparisonConfig` via `ComparisonConfig.getComparisonConfig(ConfigJSON)`. Comparators access configuration exclusively through the `ComparisonContext` they receive, which carries the `ComparisonConfig`. Key accessors:

- `isEnabled(String name)` — whether a named comparator or feature type is active
- `getAllowedGeneBiotypes()` — the biotype allowlist (empty set = all allowed)
- `getImpactLevels()` — map from comparator/feature name to `Impact` enum value

---

## Data Model

The result hierarchy mirrors the comparison hierarchy:

```
ComparisonResult                     (one per gene pair)
└── TranscriptComparisonResult[]     (one per matched transcript pair)
    └── FeatureComparisonResult[]    (one per feature type: exon, CDS, …)
        └── RegionComparisonResult[] (one per matched feature instance)
```

`ResultWriter` flattens this tree into `OutputLine` rows. Each `OutputLine` carries an integer sort key that determines its position in the output file.

---

## Dependencies

| Library | Version | Purpose |
|---|---|---|
| `gtf-utils` | 0.7 | GTF file parsing and data model (`GtfFile`, `TranscriptFeature`, etc.) |
| `htsjdk` | 4.3.0 | Interval tree implementation used for transcript overlap detection |
| `jgrapht` | 1.5.1 | Bipartite matching within transcript loci |
| `biojava-structure` | 7.2.3 | Sequence alignment and similarity calculation |
| `jackson` | 2.15.4 | JSON config file parsing |
| `commons-cli` | 1.9.0 | CLI argument parsing in runners |
| `log4j2` | 2.20.0 | Logging |
| Spring Boot | 3.4.3 | CLI mode routing via Spring `@Profile` annotations |
| minimap2 | bundled | Sequence alignment fallback for unmatched transcripts (Linux x86_64 binary in `src/main/resources/minimap2/`) |
