mkdir -p data/annotations
cd data/annotations
# B. oleracea annotations (GFF3)
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gff3/brassica_oleracea/Brassica_oleracea.BOL.61.gff3.gz
# S. sclerotiorum annotatinons (GFF3)
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-61/gff3/sclerotinia_sclerotiorum/Sclerotinia_sclerotiorum.ASM14694v1.61.gff3.gz
gunzip *gz
cd ../..
