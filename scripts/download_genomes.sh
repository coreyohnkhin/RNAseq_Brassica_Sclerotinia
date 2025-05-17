mkdir data/reference_genomes
cd data/reference_genomes
# unmasked genomes (no difference between softmasked and unmasked for in Bowtie2,
# hard masked deletes repeats)
# sclerotinia sclerotiorum
wget http://ftp.ensemblgenomes.org/pub/fungi/release-61/fasta/sclerotinia_sclerotiorum/dna/Sclerotinia_sclerotiorum.ASM14694v1.dna.toplevel.fa.gz
# b oleracea
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/brassica_oleracea/dna/Brassica_oleracea.BOL.dna.toplevel.fa.gz
cd ../..
