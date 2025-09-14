JAR=annot-compare.jar
JAVA=java

all: run

run:
	$(JAVA) -jar $(JAR) analysis --target-gtf input/ensembl_113_example.gtf --query-gtf input/ensembl_114_example.gtf --fidx input/Homo_sapiens.GRCh38.114.dna.toplevel.fa.fai --config config.json --fasta PATH/TO/FASTA --o output/113_114_test.tsv
