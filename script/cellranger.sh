
VDJ_DIR=/data2/wenkanl2/Tomas//data/2411_DuckHuman/Project_WU_VDJ/
GEX_DIR=/data2/wenkanl2/Tomas//data/2411_DuckHuman/Project_WU_GEX/
REF_DIR=./ref/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0
REF_DIRG=./ref/CAU_1_0


echo '
cellranger vdj --id=duck_G1_VDJ --denovo  \
    --inner-enrichment-primers data/Primer.txt \
    --fastqs=$VDJ_DIR \
    --sample=duck_G1_VDJ \
    --localcores=60 \
    --localmem=256 --chain IG

cellranger vdj --id=duck_Group3_VDJ --denovo  \
    --inner-enrichment-primers data/Primer.txt \
    --fastqs=$VDJ_DIR \
    --sample=duck_Group3_VDJ \ --localcores=60 \
    --localmem=256 --chain IG

cellranger vdj --id=miniHA3_G1_VDJ --denovo  \
    --inner-enrichment-primers data/Primer.txt \
    --fastqs=$VDJ_DIR \
    --sample=miniHA3_G1_VDJ \
    --localcores=60 \
    --localmem=256 --chain IG

cellranger vdj --id=miniHA3_G3_VDJ --denovo  \
    --inner-enrichment-primers data/Primer.txt \
    --fastqs=$VDJ_DIR \
    --sample=miniHA3_G3_VDJ \
    --localcores=60 \
    --localmem=256 --chain IG

'

cellranger count --id=duck_G1_GEX \
    --create-bam=false \
    --fastqs=$GEX_DIR \
    --localcores=60 \
    --localmem=256 \
    --sample=duck_G1_GEX\
    --transcriptome=$REF_DIRG

cellranger count --id=duck_Group3_GEX \
    --create-bam=false \
    --fastqs=$GEX_DIR \
    --localcores=60 \
    --localmem=256 \
    --sample=duck_Group3_GEX\
    --transcriptome=$REF_DIRG

cellranger count --id=miniHA3_G1_GEX \
    --create-bam=false \
    --fastqs=$GEX_DIR \
    --localcores=60 \
    --localmem=256 \
    --sample=miniHA3_G1_GEX\
    --transcriptome=$REF_DIRG

cellranger count --id=miniHA3_G3_GEX \
    --create-bam=false \
    --fastqs=$GEX_DIR \
    --localcores=60 \
    --localmem=256 \
    --sample=miniHA3_G3_GEX\
    --transcriptome=$REF_DIRG
