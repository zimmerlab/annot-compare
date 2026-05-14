# Annot-Compare Developer Guide

## Motivation and Vision

### The Problem

Genome annotation is not static. Projects like Ensembl release multiple annotation versions per year, and independent sources (Ensembl, RefSeq, GENCODE, and others) produce annotations for the same reference genome using different methods and evidence sets. This creates a pervasive problem in genomics research: **comparison inconsistency**.

When two studies use different annotation versions, it is not clear whether their results are directly comparable. For expression analysis, a change in exon boundaries shifts read counts even on identical sequencing data. For variant interpretation, a changed CDS boundary alters the protein consequence. More concretely:

- Is a study performed under annotation version A still valid under version B?
- If a pipeline is rerun with a new annotation, which results actually need to be recomputed and which are unchanged?
- Have genes or transcripts of interest been modified in a way that invalidates a published finding?
- Over many releases, how stable is the annotation of a given gene — is its structure settled or still frequently revised?

Existing tools either compare annotations at a coarse statistical level or rely on stable IDs across versions, an assumption that cannot be made across sources or even across consecutive Ensembl releases.

### The Vision

Annot-Compare approaches the problem from first principles: to compare annotations you must first determine which gene in one annotation corresponds to which gene in the other. This correspondence cannot be assumed from identifiers alone.

**Gene identity by function.** A gene is ultimately defined by what it encodes. If two genes across annotations represent the same biological entity, they should produce at least one similar transcript — and for coding genes, a similar protein. The primary matching criterion is therefore CDS structure similarity: two genes are candidates for comparison if they share at least one transcript with a similar coding sequence.

**Coordinate overlap as a practical filter.** Genomic position is a strong additional signal: the same gene in two annotation versions will almost always occupy the same chromosomal region. Coordinate overlap (via interval trees and bipartite matching) efficiently generates and ranks candidate transcript pairs without exhaustive pairwise comparison. For very short transcripts where overlap alone is less discriminative, this filter prevents a combinatorial explosion in candidates.

**Sequence alignment as a fallback.** A minority of transcripts will not overlap by coordinate — for example, when a gene is relocated in a new assembly or annotation boundaries shift dramatically. For these, minimap2 sequence alignment is used as a fallback to find matches on the basis of sequence similarity alone.

**Hierarchical comparison.** The ultimate goal is transcript-level and feature-level differences: what changed in the exon structure, the CDS, the UTRs. But transcript comparison requires first knowing which transcripts to pair, and that in turn requires knowing which genes correspond to each other. The comparison is therefore structured hierarchically — gene matching → transcript pairing → feature-level diffing — and results are produced at every level so they can be consumed at whichever resolution the research question demands.

---

## Building

```bash
mvn package -DskipTests
```

Output JAR: `target/annot-compare-*.jar`

To run the bundled example:

```bash
cd example && make
```

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
├── mapping/                          — Transcript and exon matching
│   ├── OverlappingTranscripts        — Interval tree overlap + bipartite matching for transcript pairing
│   ├── ExonMapping                   — Needleman-Wunsch gap alignment for exon pairing within a transcript pair
│   ├── Minimap2Bundler               — Extracts the bundled minimap2 binary to a temp directory at runtime
│   └── Minimap2Validator             — Runs minimap2 on unmapped transcripts; returns additional transcript pairs
│
├── model/                            — Data models
│   ├── TranscriptPair                — A matched target/query transcript pair
│   ├── FeaturePair                   — A matched target/query feature instance (exon, CDS, …)
│   ├── GenePair                      — A matched target/query gene
│   ├── Impact                        — Enum: HIGH, MODERATE, LOW, NONE
│   ├── config/                       — Configuration model (ComparisonConfig, ConfigJSON, FeatureConfig)
│   └── comparison/                   — Result types at each hierarchy level
│       ├── ComparisonResult          — Gene-level result; contains transcript results
│       ├── TranscriptComparisonResult — Transcript-level result; contains feature results
│       ├── FeatureComparisonResult   — Per feature-type result; contains region results
│       └── RegionComparisonResult    — Result for a single paired feature instance
│
├── parser/
│   └── FidxParser                    — Parses a FASTA index (.fai) for random-access sequence extraction
│
├── runner/                           — CLI entry points, one per mode (Spring @Profile)
│   ├── AnalysisRunner                — Main analysis mode
│   ├── FirstAnalysisRunner           — First-pass analysis variant
│   ├── SequenceExtractorRunner       — Extracts sequences from a FASTA
│   ├── AddMetaFeaturesToGTFRunner    — Adds derived features to a GTF
│   ├── GtfStatsRunner                — Reports GTF statistics
│   └── BenchmarkRunner              — Performance benchmarking
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

6. **Output** — `ResultWriter.writeComparisonResult()` flattens the result hierarchy into `OutputLine` rows, sorts them by type and content, and appends them to the output TSV.

Transcripts still unmatched after both passes are written to the `.unmapped_queries` and `.unmapped_targets` files.

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
| `jackson` | 2.15.4 | JSON config file parsing |
| `commons-cli` | 1.9.0 | CLI argument parsing in runners |
| `log4j2` | 2.20.0 | Logging |
| Spring Boot | 3.4.3 | CLI mode routing via Spring `@Profile` annotations |
| minimap2 | bundled | Sequence alignment fallback for unmatched transcripts (Linux x86_64 binary in `src/main/resources/minimap2/`) |
