wget https://ftp.ensembl.org/pub/release-113/fasta/anas_platyrhynchos_platyrhynchos/dna/Anas_platyrhynchos_platyrhynchos.CAU_duck1.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/anas_platyrhynchos_platyrhynchos/Anas_platyrhynchos_platyrhynchos.CAU_duck1.0.113.gtf.gz

cellranger mkgtf \
    Anas_platyrhynchos_platyrhynchos.CAU_duck1.0.113.gtf \
    Anas_platyrhynchos_platyrhynchos.CAU_duck1.0.113.filtered.gtf \
    --attribute=gene_biotype:protein_coding

cellranger mkref \
    --genome=CAU_1_0 \
    --genes=Anas_platyrhynchos_platyrhynchos.CAU_duck1.0.113.filtered.gtf  \
    --fasta=Anas_platyrhynchos_platyrhynchos.CAU_duck1.0.dna_sm.toplevel.fa
