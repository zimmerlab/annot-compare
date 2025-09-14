## Example Usage

To run an analysis with **Annot-Compare**, use:

```bash
java -jar annot-compare.jar analysis \
    --target-gtf <path-to-target-gtf> \
    --query-gtf <path-to-query-gtf> \
    --fasta <path-to-fasta> \
    --faidx <path-to-fai> \
    --config <path-to-config> \
    --o <path-to-output-file> \
    [--map-with-strand <boolean>]
