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
```
### Ready-to-Use Example

A complete, almost ready-to-use example (including input files and configuration) is provided in the [`example/`](example) folder.  
You can run it directly to test the program. The FASTA file needs to be downloaded separately from this [link](https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz).
Put the FASTA unzipped into the input folder to be able to run the program with the makefile. The config can be customized as desired.
