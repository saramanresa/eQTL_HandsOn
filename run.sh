#!/bin/sh

## Author: Sara Manresa
## Date: 14/10/2018
## eQTL Hands-On


## ========= Task 1 =========

echo "****"
echo "Task 1 -> wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf{.gz,.gz.tbi} --directory-prefix input/unprocessed/1000g "
echo "****"

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf{.gz,.gz.tbi} --directory-prefix input/unprocessed/1000g

echo "****"
echo "Task 1 -> cut -f1 input/unprocessed/geuvadis/geuvadis.metadata.txt | sed '1d' | sort | uniq > tmp/geuvadis.samples.txt "
echo "****"


## ========= Task 2 =========
# Get GEUVADIS samples from the metadata
cut -f1 input/unprocessed/geuvadis/geuvadis.metadata.txt | sed '1d' | sort | uniq > tmp/geuvadis.samples.txt 

echo "****"
echo "Task 2 -> (IT MAY TAKE A WHILE) bcftools view -v snps,indels -m 2 -M 2 -q 0.05:minor -S tmp/geuvadis.samples.txt -Ob input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools norm -d all -Oz -o tmp/genotypes.chr22.vcf.gz "
echo "****" 

# Subset the VCF (common samples, biallelic SNPs and indels, MAF >= 0.05, no duplicates)
bcftools view -v snps,indels -m 2 -M 2 -q 0.05:minor -S tmp/geuvadis.samples.txt -Ob input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools norm -d all -Oz -o tmp/genotypes.chr22.vcf.gz

# Q1: What do the bcftools options employed mean?
# Command: bcftools options view -v snps,indels -m 2 -M 2 -q 0.05:minor -S mean: 
# view: VCF/BCF conversion, view, subset and filter VCF/BCF files
# -v: select/exclude comma-separated list of variant types: snps,indels,mnps,ref,bnd,other : ****  in this case snps and indels
# -m: minimum/maximum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)
# -q: minimum/maximum frequency for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
# -S: file of samples to include (or exclude with "^" prefix) **** In this case: tmp/geuvadis.samples.txt
# -Ob: output type: b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v] ****  In this case: BCF
# pipe: with the result of previous command, pass to next one (bcftools norm)
# norm: Left-align and normalize indels; check if REF alleles match the reference; split multiallelic sites into multiple rows; recover multiallelics from multiple rows.
# -d all remove duplicate snps|indels|both|all|none *** in this case: all
# -Oz: output type: b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v] ****  In this case: compressed VCF
# -o:  write output to a file

echo "****"
echo "Task 2 -> PATH=$PATH:bin/"
echo "****"

PATH=$PATH:bin/

# Subset the VCF so that there are at least 10 individuals per genotype group and compress it (for indexing we require 'bgzip' compression)

echo "****"
echo "Task 2 -> filter.genotype.py -t 10 -g <(zcat tmp/genotypes.chr22.vcf.gz) | bgzip > input/processed/genotypes.chr22.vcf.gz"
echo "****"
filter.genotype.py -t 10 -g <(zcat tmp/genotypes.chr22.vcf.gz) | bgzip > input/processed/genotypes.chr22.vcf.gz

# Index the VCF
echo "****"
echo "Task 2 -> tabix -p vcf input/processed/genotypes.chr22.vcf.gz"
echo "****"
tabix -p vcf input/processed/genotypes.chr22.vcf.gz
 
# Q2:How many variants do you get in input/processed/genotypes.chr22.vcf.gz?
# Run bcftools 

echo "****"
echo "Task 2 -> bcftools stats input/processed/genotypes.chr22.vcf.gz > file.stats"
echo "****"
bcftools stats input/processed/genotypes.chr22.vcf.gz > file.stats
#También se puede hacer abriendo archivo con un less -S  y luego con un zcat   pipe -v "#" pipe less -S
# Y después otro zcat pipe -v "#" pipe wc -l ?
#
# Results in file.stats (after)
# SN      0       number of samples:      445
# SN      0       number of records:      74656
# SN      0       number of no-ALTs:      0
# SN      0       number of SNPs: 65324
# SN      0       number of MNPs: 0
# SN      0       number of indels:       9332
# SN      0       number of others:       0
# SN      0       number of multiallelic sites:   0
# SN      0       number of multiallelic SNP sites:       0

#RESPUESTA: records:      74656
# With indels + SNPs = records --> 9332 + 65324
    
    
# Q3: How many samples do you have before and after subsetting? 
# Hint: The genotype file contains information from many samples, 
# but only a subset of them has gene expression data in the GEUVADIS project. 
# Note that the GEUVADIS metadata contains duplicated sample IDs, as gene expression 
# in GEUVADIS is calculated from paired-end RNA-Seq data. You can use bcftools to print the sample IDs and count them

# We know that tmp/genotypes.chr22.vcf.gz has genotype information BEFORE subsetting

echo "****"
echo "Task 2 -> bcftools stats tmp/genotypes.chr22.vcf.gz  > filebefore.stats"
echo "****"
bcftools stats tmp/genotypes.chr22.vcf.gz  > filebefore.stats
#before es con el archivo ALL, vamos a la carpeta input, que es donde esta, y realizo:
echo "****"
echo "Task 2 -> bcftools stats input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > filebefore.stats"
echo "****"
bcftools stats input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > filebefore.stats
# ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 
# Before subsetting:

#SN    [2]id   [3]key  [4]value
#SN      0       number of samples:      2504
#SN      0       number of records:      1103547
#SN      0       number of no-ALTs:      0
#SN      0       number of SNPs: 1060388
#SN      0       number of MNPs: 0
#SN      0       number of indels:       43230
#SN      0       number of others:       801
#SN      0       number of multiallelic sites:   6348
#SN      0       number of multiallelic SNP sites:       4063

# Before subsetting we have 2504 and 1103547
# After subsetting we have 445 and 74656           

## ========= Task 3 =========
#(Q3 en statistics)
#Q1: Which version of GENCODE is GEUVADIS using? 
#usa la version 12v

#Q2: To which genome assembly does this annotation correspond? 
# The genome assembly is GRCh37 

#Q3: How many protein coding genes are annotated in the last version (v29)?
# in the last version v29 with genome assembly version GRCh38
# Entonces vamos en el GENCODE de v29 y allí accedemos a STATISTICS
# EN STATISTICS observamos que el nº de 19940 protein coding genes 

#Q4: Which command do you use to do this?
# export PATH=$PATH:$PWD/bin


# Set the variable 'release' to the version of GENCODE used in GEUVADIS (e.g. `release=99`) and download the corresponding GENCODE annotation
echo "****"
echo "Task 3 -> release=12"
echo "****"
release=12

echo "****"
echo "Task 3 -> wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$release/gencode.v$release.annotation.gtf.gz"
echo "****"
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$release/gencode.v$release.annotation.gtf.gz

echo "****"
echo "Task 3 -> gencode.v$release.annotation.gtf.gz input/unprocessed/gencode/gencode.annotation.gtf.gz"
echo "****"
mv gencode.v$release.annotation.gtf.gz input/unprocessed/gencode/gencode.annotation.gtf.gz

# Obtain a BED file from the GTF, selecting just the 'protein coding' and 'lincRNA' genes
echo "****"
echo "Task 3 -> zcat input/unprocessed/gencode/gencode.annotation.gtf.gz | grep \"gene_type \"protein_coding\"\|gene_type \"lincRNA\"\" | gtf2bed.sh > tmp/gencode.annotation.bed"
echo "****"
zcat input/unprocessed/gencode/gencode.annotation.gtf.gz | grep "gene_type \"protein_coding\"\|gene_type \"lincRNA\"" | gtf2bed.sh > tmp/gencode.annotation.bed

#Q5: But how to get the TSS ((transcription start site) positions and the gene lengths from it?
# Primero miramos la ultima columna del fichero BED, que indica si es positivo o negativo la cadena 
# la posición inicial del tss es la primera posición de la cadena y para ello miramos la segunda columna que indica la posición inicial (start position) 
# Por ejemplo, tenemos: 
# chr1    36080   36081   ENSG00000237613.2       .       -
# chr1    69090   69091   ENSG00000186092.4       .       + 
# chr1    91104   91105   ENSG00000239945.1       .       -
# La primera posición del principio de una cadena en este ejemplo es 69090
# La longitud es la resta de la posición final del final de la cadena y el inicio de la cadena, aqui seria 91105 - 69090 

#Q6: to which BED coordinates would correspond the GTF coordinates chr1 10 20?
# To BED positions [9,20) because BED starts with 0 and gtf not. 


# Compute gene lengths 
echo "****"
echo "Task 3 -> awk 'BEGIN{OFS=\"\t\"}{print $1,$2,$3,$4,$3-$2,$6}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed"
echo "****"
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$3-$2,$6}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed

# Compute TSS positions. Note that for genes in the '+' strand, the TSS is the start position, and for genes in the '-' strand it is the end position!
echo "****"
echo "Task 3 -> awk 'BEGIN{OFS=\"\t\"}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}else{print $1,$3-1,$3,$4,$5,$6}}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed"
echo "****"
awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}else{print $1,$3-1,$3,$4,$5,$6}}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed

# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
echo "****"
echo "Task 3 -> sed -i \"s/^chr//\" tmp/gencode.annotation.bed"
echo "****"
sed -i "s/^chr//" tmp/gencode.annotation.bed


#Q7: Why do we need to use tmpfile below?
# usamos el tempfile para que el siguiente comando siempre se base en el anterior
# trampa: en este ejemplo el4 tempfile machaca con el mv el fichero tmp/gencode.annotation.bed, por lo que no haria falta

## ========= Task 4 =========

# Join the bed file with the expression file
# (Both files should be row-ordered by gene ID. Column order and header are lost in the output file)
echo "****"
echo "Task 4 -> join -1 4 -2 1 -t $'\t' <(sort -k4,4 tmp/gencode.annotation.bed) <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | sort -k1,1) > tmp/joint.tsv"
echo "****"
join -1 4 -2 1 -t $'\t' <(sort -k4,4 tmp/gencode.annotation.bed) <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | sort -k1,1) > tmp/joint.tsv

# Subset chr22 (same as the VCF file)
echo "****"
echo "Task 4 -> tmp/joint.tsv > tmp/joint.chr22.tsv"
echo "****"
awk '$2==22' tmp/joint.tsv > tmp/joint.chr22.tsv

# Recover the column order, sort rows by chr and start position (WARNING: this command may not work within the docker container for WSL users)
echo "****"
echo "Task 4 -> paste <(awk 'BEGIN{OFS=\"\t\"}{print $2,$3,$4,$1,$5,$6}' tmp/joint.chr22.tsv) <(cut -f1-6 --complement tmp/joint.chr22.tsv) | sort -k1,1V -k2,2n > tmp/joint.chr22.bed"
echo "****"
paste <(awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6}' tmp/joint.chr22.tsv) <(cut -f1-6 --complement tmp/joint.chr22.tsv) | sort -k1,1V -k2,2n > tmp/joint.chr22.bed

# Recover the header
echo "****"
echo "Task 4 -> cat <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | head -1 | sed \"s/TargetID/#chr\tstart\tend\tgene\tlength\tstrand/\") tmp/joint.chr22.bed > tmp/genes.chr22.rpkm.bed"
echo "****"
cat <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | head -1 | sed "s/TargetID/#chr\tstart\tend\tgene\tlength\tstrand/") tmp/joint.chr22.bed > tmp/genes.chr22.rpkm.bed

#  Q1: Of all genes considered, which have lower expression levels, protein-coding or lincRNA? 
# la menor expresion son las que valen 0 ya que no hay expresion, no hay línea celular 
# lincRNA tiene menos expression ya que la proteina no esta codificada, por lo tanto el protein-coding tendra mas expression.

# Q2: Why do we need gene expression to be normal? 
# para evitar alteraciones en la expresion del gen (si no es normalizada no sabremos si esta alterada o no)

# Q3: How would you check that quantile normalization worked? 
# miramos que sean valores normales introduciendo los máximos y minimos que corresponden a las alteraciones de decrecimiento de linea celular. 

# Q4: and that gene expression of a gene follows a normal distribution?
# miramos que la expresión se encuentra dentro del intervalo de distribución normalizada, es decir, que no se encuentra en alteración

# Run gene expression normalization (quantile normalization + gene expression to normal distribution)
# Filter out genes with less than 0.1 RPKM in 50% of the samples
echo "****"
echo "Task 4 -> normalize.R -i tmp/genes.chr22.rpkm.bed -o tmp/genes.chr22.norm.bed"
echo "****"
normalize.R -i tmp/genes.chr22.rpkm.bed -o tmp/genes.chr22.norm.bed

# Compress and index the final gene expression file
echo "****"
echo "Task 4 -> bgzip tmp/genes.chr22.norm.bed"
echo "****"
bgzip tmp/genes.chr22.norm.bed

echo "****"
echo "Task 4 -> tabix -p bed tmp/genes.chr22.norm.bed.gz"
echo "****"
tabix -p bed tmp/genes.chr22.norm.bed.gz

echo "****"
echo "Task 4 -> mv tmp/genes.chr22.norm.bed.gz* input/processed"
echo "****"
mv tmp/genes.chr22.norm.bed.gz* input/processed

## ========= Task 5 =========

# Before:
echo "****"
echo "Task 5 -> check.norm.R -i tmp/genes.chr22.rpkm.bed -o result/plots/check.not.norm.pdf"
echo "****"
# utilizamos el check.norm para generar un pdf igual que el de after y compararlos. Como input usamos sin embargo el fichero sin normalizar. 
check.norm.R -i tmp/genes.chr22.rpkm.bed -o result/plots/check.not.norm.pdf

# After:
echo "****"
echo "Task 5 -> check.norm.R -i input/processed/genes.chr22.norm.bed.gz -o result/plots/check.norm.pdf"
echo "****"
check.norm.R -i input/processed/genes.chr22.norm.bed.gz -o result/plots/check.norm.pdf

# Q1: What can you see?
# vemos que en la expresion del gen normalizada podemos ver las medianas y predeterminar cuales estan alteradas o fuera de rango. Sin embargo en la 
# expresion del gen sin normalizar no podemos determinar la mediana ni saber cuales estan alteradas o no. Tambien vemos el que el q-q plot en 
# caso normalizado tiene una distribucion lineal y en el sin normalizar hay una distribucion logaritmica. 

## ========= Task 6 =========

# we have observed all files (with -cd we decompress just for next command using gzip -cd)
echo "****"
echo "Task 6 -> gzip -cd input/processed/genes.chr22.norm.bed.gz | head -1 | awk 'BEGIN{FS=OFS=\"\t\"}{for (i=1;i<=NF;i++) {print i, $i}}' "
echo "****"
gzip -cd input/processed/genes.chr22.norm.bed.gz | head -1 | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'
echo "****"
echo "Task 6 -> gzip -cd input/processed/genotypes.chr22.vcf.gz | head -1 | awk 'BEGIN{FS=OFS=\"\t\"}{for (i=1;i<=NF;i++) {print i, $i}}' "
echo "****"
gzip -cd input/processed/genotypes.chr22.vcf.gz | head -1 | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'
echo "****"
echo "Task 6 -> head -1 input/unprocessed/geuvadis/geuvadis.metadata.txt | awk 'BEGIN{FS=OFS=\"\t\"}{for (i=1;i<=NF;i++) {print i, $i}}' "
echo "****"
head -1 input/unprocessed/geuvadis/geuvadis.metadata.txt | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'
echo "****"
echo "Task 6 -> head -1 input/unprocessed/1000g/1000g.phase3_metadata.txt | awk 'BEGIN{FS=OFS=\"\t\"}{for (i=1;i<=NF;i++) {print i, $i}}' "
echo "****"
head -1 input/unprocessed/1000g/1000g.phase3_metadata.txt | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'


#Q1: Which ones would you select? Hint: Maybe you find useful this small oneliner:
# we use input/unprocessed/1000g/1000g.phase3_metadata.txt and genotypes.chr22.vcf (decompressed)
# after executing head,... we obtain: 
#1	sample
#2	pop
#3	super_pop
#4	gender


## ========= Task 7 =========

# Expression PCA

echo "****"
echo "Task 7 -> QTLtools pca --bed input/processed/genes.chr22.norm.bed.gz --scale --center --out result/expression "
echo "****"
QTLtools pca --bed input/processed/genes.chr22.norm.bed.gz --scale --center --out result/expression 

# Genotypes PCA
echo "****"
echo "Task 7 -> QTLtools pca --vcf input/processed/genotypes.chr22.vcf.gz --scale --center --maf 0.05 --distance 50000 --out result/genotypes "
echo "****"
QTLtools pca --vcf input/processed/genotypes.chr22.vcf.gz --scale --center --maf 0.05 --distance 50000 --out result/genotypes

# Q1: What do the parameters employed mean? 
# pca: Calculate principal components for a BED/VCF/BCF file
# --bed: Phenotypes in BED format
# --vcf: Genotypes in VCF/BCF/BED format
# --scale: Scale the quantifications or genotypes to unit variance before PCA.
# --center: Center the quantifications or genotypes before PCA.
# --out: Output file prefix
# --maf: Exclude sites with MAF less than this.
# --distance: Only include sites separated with this many bp (bitmap pixels?)

 
# Q2: Which information do the output files contain?
# In .pca file, the header line gives the sample IDs, the first column the ID of the PCs, starting from the first one and each successive line the coordinates of the samples on the PCs.
# In .pca_stats file, The 3 first lines give you for each successive PC:
 #   1. The standard deviation
 #  2. The variance explained
 #  3. The cumulative variance explained

# Note the input should coincide with the output of QTLtools pca on each case
echo "****"
echo "Task 7 -> pcaPlot.R -i result/expression -o result/plots/expression.pca.pdf "
echo "****"
pcaPlot.R -i result/expression -o result/plots/expression.pca.pdf

echo "****"
echo "Task 7 -> pcaPlot.R -i result/genotypes -o result/plots/genotypes.pca.pdf "
echo "****"
pcaPlot.R -i result/genotypes -o result/plots/genotypes.pca.pdf

# Q3: What can you observe in the plots?
# observamos que en el expression.pca.pdf hay una dispersion muy amplia, es decir, hay poca relacion entre PC1 y PC2 
# sin embargo, en genotypes vemos en el diagrama de dispersion hay mas relacion entre las variables PC1 y PC2, sobretodo cuando PC1 es menor que 0 y PC2 cualquier valor. 
# o un PC1 alrededor de 20 y PC2 esta entre 5 y -5. En concreto, en genotypes tenemos una alta correlacion. 

echo "****"
echo "Task 7 -> pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color super_pop --out result/plots/genotypes.pca.super_pop.pdf "
echo "****"
pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color super_pop --out result/plots/genotypes.pca.super_pop.pdf

echo "****"
echo "Task 7 -> pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color pop --out result/plots/genotypes.pca.pop.pdf "
echo "****"
pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color pop --out result/plots/genotypes.pca.pop.pdf

echo "****"
echo "Task 7 -> pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color gender --out result/plots/genotypes.pca.gender.pdf "
echo "****"
pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color gender --out result/plots/genotypes.pca.gender.pdf

# Q4: With this information, which covariates seem more relevant to explain the variability in the data?
# pop, super_pop and gender
# in pop we can mesure the proteins: CEU, FIN, GBR, TSI, YRI

# Generate a common metadata with info about the population, gender and laboratory.
echo "****"
echo "Task 7 -> join -j 1 -t $'\t' <(sort -k1,1 input/unprocessed/1000g/1000g.phase3_metadata.txt) <(cut -f1,20 input/unprocessed/geuvadis/geuvadis.metadata.txt | sort -k1,1 | uniq) > tmp/metadata.txt "
echo "****"
join -j 1 -t $'\t' <(sort -k1,1 input/unprocessed/1000g/1000g.phase3_metadata.txt) <(cut -f1,20 input/unprocessed/geuvadis/geuvadis.metadata.txt | sort -k1,1 | uniq) > tmp/metadata.txt

# Set names for the new metadata
echo "****"
echo "Task 7 -> sed -i '1s/^/sampleID\tpop\tsuper_pop\tgender\tlab\n/' tmp/metadata.txt "
echo "****"
sed -i '1s/^/sampleID\tpop\tsuper_pop\tgender\tlab\n/' tmp/metadata.txt

# Build a linear model and plot the contribution of each factor in the metadata to the total variance
echo "****"
echo "Task 7 -> var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m tmp/metadata.txt --formula \"~ (1|gender) + (1|pop) + (1|lab)\" -o result/plots/vp.pdf "
echo "****"
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m tmp/metadata.txt --formula "~ (1|gender) + (1|pop) + (1|lab)" -o result/plots/vp.pdf

# Q5: Which are the factors that explain more variance?
# The factor that explains more variance is the "Residuals" as we can see in vp.pdf

## ========= Task 8 =========

# Compute 10 PEER factors
echo "****"
echo "Task 8 -> peer.R -i input/processed/genes.chr22.norm.bed.gz -p 10 -o tmp/peer.tsv "
echo "****"
peer.R -i input/processed/genes.chr22.norm.bed.gz -p 10 -o tmp/peer.tsv

# Check how much variance do the first 5 PEER explain in comparison with the known factors
echo "****"
echo "Task 8 -> var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m <(paste tmp/peer.tsv tmp/metadata.txt) -f \"~ (1|pop) + (1|lab) + PEER1 + PEER2 + PEER3 + PEER4 + PEER5\" -o result/plots/vp.peer.pdf "
echo "****"
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m <(paste tmp/peer.tsv tmp/metadata.txt) -f "~ (1|pop) + (1|lab) + PEER1 + PEER2 + PEER3 + PEER4 + PEER5" -o result/plots/vp.peer.pdf

# 'Rscript -e' is just a trick to run an R script without opening an interactive R session in the console. ;)
echo "****"
echo "Task 8 -> join -j 1 -t $'\t' tmp/metadata.txt tmp/peer.tsv  | Rscript -e 'write.table(t(read.table(file(\"stdin\", open = "r", blocking = T), h = F)), file = \"input/processed/covariates.tsv\", quote = F, sep = \"\t\", col.names = F, row.names = F)' "
echo "****"
join -j 1 -t $'\t' tmp/metadata.txt tmp/peer.tsv  | Rscript -e 'write.table(t(read.table(file("stdin", open = "r", blocking = T), h = F)), file = "input/processed/covariates.tsv", quote = F, sep = "\t", col.names = F, row.names = F)'

# Compress it
echo "****"
echo "Task 8 -> gzip input/processed/covariates.tsv "
echo "****"
gzip -f input/processed/covariates.tsv

# Q1: How much variance do they explain? On average is it more or less than the explained by the known factors?
# we have 5 hidden covariates: on average is more than the explained known factors but we lost "gender" factor. 


## ========= Task 9 =========


echo "****"
echo "Task 9 -> QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --nominal 1 --out result/nominals.txt "
echo "****"
# we put --nominal 1 instead of 0.01 because exercise requests from 0,1
QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --nominal 1 --out result/nominals.txt

echo "****"
echo "Task 9 -> pvdist.R -i result/nominals.txt --col 12 -o result/plots/pvdist.pdf "
echo "****"
pvdist.R -i result/nominals.txt --col 12 -o result/plots/pvdist.pdf

# Q1: Are there pairs genotype-phenotype with exactly the same p-value and effect size (β)? How is this possible?
# Si, ya que como podemos en el histograma, por ejemplo, con p-value = 0.2 hay aproximadamente 60000 effect-size, mismo effect-size con aproximadamente 0.4 o 0.6
# una posible causa de que haya tanta similitud de p-values i effect size es debido a que partimos de una normalizacion de los genes 

# Q2: What do you observe?
# en el grafico logaritmico observamos que a partir de 2 p-values expected tenemos un incremento logaritmico respecto a los observados
# en el histograma vemos que todos los p-value tienen una frecuencia aproximada de entre 5000 y 6000

# Q3: Which SNPs did you select? What do you observe?
# we select two lines with same p-value = 0.700757, snps are rs9265 and rs165655
# ENSG00000237476.1 22 20957092 20957092 + 2926 -999461 rs9265 22 19957631 19957631 0.700757 0.0156431 0
# ENSG00000237476.1 22 20957092 20957092 + 2926 -999329 rs165655 22 19957763 19957763 0.700757 0.0156431 0
echo "****"
echo "Task 9 -> plink --ld rs9265 rs165655 --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/ld2 "
echo "****"
plink --ld rs9265 rs165655 --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/ld2

# Result: 
# Note: No phenotypes present.
#
#--ld rs9265 rs165655:
#
#   R-sq = 1              D' = 1
#
#   Haplotype     Frequency    Expectation under LE
#   ---------     ---------    --------------------
#          CG      0.389888                0.152012
#          AG     -0                       0.237875
#          CA     -0                       0.237875
#          AA      0.610112                0.372237
#
#   In phase alleles are CG/AA

# Observamos que no tenemos presencia de fenotipos y que tenemos frecuencias de haplotypes baja, es decir, que tenemos una recombinacion genetica de baja frecuencia


## ========= Task 10 =========

# generate with 1 thread
#echo "****"
#echo "Task 10 -> QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --out result/permutations.txt "
#echo "****"
#QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --out result/permutations.txt

# if we want to generate it with 4 threads
echo "****"
echo "Task 10 -> Generating permutations with 4 threads in result/permutations.txt "
echo "****"
for j in $(seq 1 16); do
  echo "cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --chunk $j 16 --out result/permutations_$j.txt"
done | xargs -P4 -n14 QTLtools
cat result/permutations_*.txt > result/permutations.txt; rm result/permutations_*.txt 

# we preapre commands for R script 
echo "****"
echo "Task 10 -> preparing commands for R script using printf... and saving into R.permutation "
echo "****"
printf "p <- read.table(\"result/permutations.txt\")                                                      # Read input file\n
pdf(\"result/plots/pv-correlation.pdf\",  paper = 'a4r', width = 9, height = 6)                   # Open PDF device\n
plot(p[, 18], p[, 19], xlab = \"pv (perm)\", ylab = \"pv (beta)\")                                  # Plot p-values\n
abline(0, 1, col = \"red\")                                                                       # Add red line 1=1\n
plot(-log10(p[, 18]), -log10(p[, 19]), xlab = \"-log10 pv (perm)\", ylab = \"-log10 pv (beta)\")    # Repeat in -log10 space to check the behaviour of the small p-values.\n
abline(0, 1, col = \"red\")\n
dev.off()                                                                                       # Close device\n
quit(\"no\")\n" > R.permutation

# We run R to generate the correlation pdf
echo "****"
echo "Task 10 -> R --no-save < R.permutation "
echo "****"
R --no-save < R.permutation

## ========= Task 11 =========

# we run:

echo "****"
echo "Task 11 -> mtc.R -n result/nominals.txt -p result/permutations.txt -m bonferroni --alpha=0.05 --out tmp/bonferroni.txt "
echo "****"
mtc.R -n result/nominals.txt -p result/permutations.txt -m bonferroni --alpha=0.05 --out tmp/bonferroni.txt

echo "****"
echo "Task 11 -> mtc.R -n result/nominals.txt -p result/permutations.txt -m fdr --alpha=0.05 --out tmp/fdr.txt "
echo "****"
mtc.R -n result/nominals.txt -p result/permutations.txt -m fdr --alpha=0.05 --out tmp/fdr.txt

echo "****"
echo "Task 11 -> mtc.R -n result/nominals.txt -p result/permutations.txt -m perm-fdr --alpha=0.05 --out result/eqtls.tsv "
echo "****"
mtc.R -n result/nominals.txt -p result/permutations.txt -m perm-fdr --alpha=0.05 --out result/eqtls.tsv

# Q1: How many significant eQTLs do we find in each case in comparison with the nominal pass?
# we run wc to count lines in each file 
echo "****"
echo "Task 11 -> wc -l result/nominals.txt "
echo "****"
wc -l result/nominals.txt
echo "****"
echo "Task 11 -> wc -l tmp/bonferroni.txt "
echo "****"
wc -l tmp/bonferroni.txt
echo "****"
echo "Task 11 -> wc -l tmp/fdr.txt "
echo "****"
wc -l tmp/fdr.txt
echo "****"
echo "Task 11 -> wc -l result/eqtls.tsv "
echo "****"
wc -l result/eqtls.tsv
# results: 
# 1406669 result/nominals.txt
# 6008 tmp/bonferroni.txt
# 17268 tmp/fdr.txt 
# 12032 result/eqtls.tsv
# In comparison with nominal pass we have much less. 

## ========= Task 12 =========
echo "****"
echo "Task 12 -> eQTLviewer.R -i <(head -n 10 result/eqtls.tsv) -g input/processed/genotypes.chr22.vcf.gz -e input/processed/genes.chr22.norm.bed.gz -o result/plots/eQTLs_head.pdf --verbose "
echo "****"
eQTLviewer.R -i <(head -n 10 result/eqtls.tsv) -g input/processed/genotypes.chr22.vcf.gz -e input/processed/genes.chr22.norm.bed.gz -o result/plots/eQTLs_head.pdf --verbose

## ========= Task 13 =========

echo "****"
echo "Task 13 -> rsync -av rsync://ftp.ensembl.org/ensembl/pub/grch37/release-86/regulation/homo_sapiens/AnnotatedFeatures.gff.gz input/unprocessed/ensembl "
echo "****"
# Download from ftp server
rsync -av rsync://ftp.ensembl.org/ensembl/pub/grch37/release-86/regulation/homo_sapiens/AnnotatedFeatures.gff.gz input/unprocessed/ensembl

# Get chr, start, end and feature name in BED format
echo "****"
echo "Task 13 -> (IT MAY TAKE A WHILE) zcat input/unprocessed/ensembl/AnnotatedFeatures.gff.gz | awk 'BEGIN{FS=OFS=\"\t\"}{print $1, $4-1, $5, $9}' | sed -r 's/Name=([^;]+);.*/\1/' | grep -v '^GL' | sort -V > tmp/ERB.bed "
echo "****"
zcat input/unprocessed/ensembl/AnnotatedFeatures.gff.gz | awk 'BEGIN{FS=OFS="\t"}{print $1, $4-1, $5, $9}' | sed -r 's/Name=([^;]+);.*/\1/' | grep -v '^GL' | sort -V > tmp/ERB.bed

# Merge overlapping features of the same type 
# e.g. chr1 100 200 feat1            chr1 100 300 feat
#      chr1 150 300 feat1     =>     chr1 100 250 feat2
#      chr1 100 250 feat2
echo "****"
echo "Task 13 -> (IT MAY TAKE A WHILE) Merging into input/processed/ERB.collapsed.bed "
echo "****"
for feat in $(cut -f4 tmp/ERB.bed | sort | uniq); do 
  bedtools merge -i <(grep -Fw $feat tmp/ERB.bed) -c 4 -o distinct
done > input/processed/ERB.collapsed.bed

# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
echo "****"
echo "Task 13 -> (IT MAY TAKE A WHILE) sed -i "s/^chr//" input/processed/ERB.collapsed.bed "
echo "****"
sed -i "s/^chr//" input/processed/ERB.collapsed.bed

echo "****"
echo "Task 13 -> (IT MAY TAKE A WHILE) Enriching into result/enrichments.txt "
echo "****"
for feat in $(cut -f4 input/processed/ERB.collapsed.bed | sort | uniq); do 
  QTLtools fenrich --qtl <(sed '1d' result/eqtls.tsv | awk '{print $9, $10-1, $10, $8, $1, "."}') --tss tmp/gencode.annotation.bed  --bed <(grep -Fw $feat input/processed/ERB.collapsed.bed) --out tmp/enrich.txt > /dev/null; echo "$(cat tmp/enrich.txt) $feat" 
done | grep -Fwv inf | grep -Fwv nan > result/enrichments.txt

echo "****"
echo "Task 13 -> plot.enrich.R -i result/enrichments.txt -o result/plots/enrich.pdf "
echo "****"
plot.enrich.R -i result/enrichments.txt -o result/plots/enrich.pdf

# Q1: Which are the top enriched features? Which kind of factors are they? 
# enriched feature H3K36me3, H3K4me1, PolII, DNase1, H3K79me2, H3K4me1 are the most enriched features because they have a higher p-value 
# they are proteins

# Q2: What does an odds ratio lower than one mean?
# la proteina con odds ration < 1 tiene menos probabilidades de ejercer su funcion caracteristica a pesar del enriquecimiento (la proteina en este caso H3K27me3, lini1, BCLAF1 y ZBTB7A sera no funcional)
## ========= Task 14 =========

echo "****"
echo "Task 14 -> sed '1d' result/eqtls.tsv | cut -f8 | sort | uniq > tmp/eqtls_snps.tsv "
echo "****"
sed '1d' result/eqtls.tsv | cut -f8 | sort | uniq > tmp/eqtls_snps.tsv

# !!!! WARNING: El VEP no funcionava con todos los eqtls_snps.tsv de modo que solo hemos procesado los 100 primeros. Entre los 100 primeros no hay ninguno con HIGH impact.  

# Q1: Which kind of consequences have they, according to the VEP? In which proportion?

# Consequences: 
# intron_variant: 33%
# downstream_gene_variant: 25%
# upstream_gene_variant: 13%
# non_coding_transcript_variant: 10%
# NMD_transcript_variant: 5%
# 3_prime_UTR_variant: 5%
# regulatory_region_variant: 2%
# non_coding_transcript_exon_variant: 2%
# synonymous_variant: 2%
# Others

# Q2: How many eQTLs are high impact variants? Which consequences are related to those high impact variants?
# There is no HIGH impact, but there are three variants with MODERATE impact, the consequences of those variants are "missense_variant" which is the 38% of the coding consequences


# Q3: Out of all high impact variants, how many of them are falling in acceptor splice sites of protein coding genes?
# we filter using "Intron is defined" to detect variants with splice and we get 76 variants. 


## ========= Task 15 =========

# Q1: In which biological processes are your eGenes enriched? Which molecular functions and components correspond to those processes?
# Processes: 
# - (GO:0032496) response to lipopolysaccharide	
# - (GO:0002237) response to molecule of bacterial origin	
# - (GO:0042445) hormone metabolic process	
# Functions: 
# - (GO:0005088) Ras guanyl-nucleotide exchange factor activity	
# Component: 
# - (GO:0005783) endoplasmic reticulum

# Generate a list of sGenes
echo "****"
echo "Task 15 -> cut -f1 result/eqtls.tsv | sed '1d' | sed 's/\..\+//' | sort | uniq > tmp/egenes.txt "
echo "****"
cut -f1 result/eqtls.tsv | sed '1d' | sed 's/\..\+//' | sort | uniq > tmp/egenes.txt

# We will use as background all the genes (PC and lincRNA) in chr22
echo "****"
echo "Task 15 -> awk '{if($1==22) print $4}' tmp/gencode.annotation.bed | sed 's/\..\+//' | sort | uniq > tmp/bg.txt "
echo "****"
awk '{if($1==22) print $4}' tmp/gencode.annotation.bed | sed 's/\..\+//' | sort | uniq > tmp/bg.txt

## ========= Task 16 =========

# Generate input files for QTLtools rtc
echo "****"
echo "Task 16 -> grep -Fwf <(cut -f1 result/eqtls.tsv ) result/permutations.txt > tmp/rtc_input "
echo "****"
grep -Fwf <(cut -f1 result/eqtls.tsv ) result/permutations.txt > tmp/rtc_input

echo "****"
echo "Task 16 -> cut -f4,7 input/unprocessed/gwas/gwas.catalog.hg19.bed > tmp/gwas_trait "
echo "****"
cut -f4,7 input/unprocessed/gwas/gwas.catalog.hg19.bed > tmp/gwas_trait

# Download the file 'hotspots_b37_hg19.bed' from QTLtools website
echo "****"
echo "Task 16 -> wget http://jungle.unige.ch/QTLtools_examples/hotspots_b37_hg19.bed --directory-prefix tmp "
echo "****"
wget http://jungle.unige.ch/QTLtools_examples/hotspots_b37_hg19.bed --directory-prefix tmp

# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
echo "****"
echo "Task 16 -> sed -i 's/^chr//' tmp/hotspots_b37_hg19.bed "
echo "****"
sed -i 's/^chr//' tmp/hotspots_b37_hg19.bed

# Run RTC
echo "****"
echo "Task 16 -> QTLtools rtc --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --hotspot tmp/hotspots_b37_hg19.bed --gwas-cis tmp/gwas_trait tmp/rtc_input --out result/rtc.txt "
echo "****"
QTLtools rtc --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --hotspot tmp/hotspots_b37_hg19.bed --gwas-cis tmp/gwas_trait tmp/rtc_input --out result/rtc.txt

# Q1: How many pairs of variants have a RTC value above 0.9? 
# we have 37 pairs of variants with RT higher than 0.9 

# Q2: For each pair, we have a GWAS hit and an eQTL. Find the trait/disease to which the GWAS variant is associated, and the gene to which the eQTL is associated. Does this gene have any relationship to the trait/disease? Explore the biological databases that you know to gather more information. 
# Yes, for example, the GWAS covariant rs5747327 associated to eQTL covariant rs10537282 are used for the treatment of Granulocyte count and this covariant and myeloid white cell count. 
# Also, you can check the file tmp/gwas_trait and check that rs5747327 is treating Granulocyte. 


# Q3: Which consequences, according to the variant effect predictor, do these variants have?
# Para obtener solo los covariantes usados en el rtc. ejecutamos: 
# cat result/rtc.txt |  awk '{ print $1 } { print $2 }' | tail -n+3 > tmp/rtc_covariants.txt
# COn eso cojemos las dos primeras columnas de rtc.txt y las juntamos para poder passarlo al VEP.
# Consequences: 
# intron_variant: 39%
# upstream_gene_variant: 20%
# downstream_gene_variant: 17%
# non_coding_transcript_variant: 8%
# NMD_transcript_variant: 6%
# regulatory_region_variant: 3%
# non_coding_transcript_exon_variant: 3%
# synonymous_variant: 2%
# missense_variant: 1%
# Others

## ========= Task 17 =========


# Generate the ID/Z-scores input file. Select your favourite gene (e.g. gene=ENS00000000000.0).
# Set k (number of variants) to 50
echo "****"
echo "Task 17 -> gene=ENSG00000237476.1 "
echo "****"
gene=ENSG00000237476.1

echo "****"
echo "Task 17 -> compZscore.R --gene $gene --nominal result/nominals.txt -k 50 --output tmp/$gene.rs_z "
echo "****"
compZscore.R --gene $gene --nominal result/nominals.txt -k 50 --output tmp/$gene.rs_z

# Generate the LD matrix 
echo "****"
echo "Task 17 -> plink --r square --snps $(cut -f1 tmp/$gene.rs_z) --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/$gene "
echo "****"
plink --r square --snps $(cut -f1 tmp/$gene.rs_z) --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/$gene

# run caviar
echo "****"
echo "Task 17 -> CAVIAR -z tmp/$gene.rs_z -l tmp/$gene.ld -o result/$gene "
echo "****"
CAVIAR -z tmp/$gene.rs_z -l tmp/$gene.ld -o result/$gene

# Q1: How many variants are there in the credible (ρ=0.95) set? For each of these variants, which is the probability to be causal?
# Our result/ENSG00000237476.1_set indicates that the following variants are credible (p=0.95 by default) 
# rs46455
# rs178250
# rs362232
# We check the file result/ENSG00000237476.1_post to check the casual probabilities: 
# rs46455 has 0.148131 probability to be casual in the set
# rs178250 has 0.335514 probability to be casual in the set
# rs362232 has 0.483646 probability to be casual in the set

# Q2: Which are the p-values and effect sizes of these variants? How are they in comparison to the p-values and effect sizes of other variants tested for the same gene? 
# Ejecutamos un grep sobre el nominals.txt de nuestro gen y nuestras variantes para buscar sus p-values 
# grep 'rs46455' result/nominals.txt | grep 'ENSG00000237476'
# Obtenemos ENSG00000237476.1 22 21311380 21311380 + 2773 -320 rs46455 22 21311060 21311060 7.50339e-05 0.222139 1
# p-value = 7.50339e-05
# grep 'rs178250' result/nominals.txt | grep 'ENSG00000237476'
# Obtenemos ENSG00000237476.1 22 21311380 21311380 + 2773 -1045 rs178250 22 21310335 21310335 0.000103796 0.217946 0
# p-value = 0.000103796
# grep 'rs362232' result/nominals.txt | grep 'ENSG00000237476'
# Obtenemos ENSG00000237476.1 22 21311380 21311380 + 2773 -85966 rs362232 22 21225413 21225414 0.000180008 -0.260878 0
# p-value = 0.000180008
# Hemos ejecutado el siguiente comando para ver los p-values de las otras variantes de nuestro gen para compararlo: 
# grep 'ENSG00000237476' result/nominals.txt | grep -v 'rs46455' | grep -v 'rs178250' | grep -v 'rs362232' | awk '{print $12}'
# Podemos ver que p-value es mayor en la mayoria o casi todas de las otras variantes respecto a las nuestras y tambien tienen un effect-size mucho mayor. 


# Q4: Which consequences, according to the variant effect predictor, do these variants have?
# Consequences:
# downstream_gene_variant: 64%
# intron_variant: 18%
# upstream_gene_variant: 18%

## ========= Task 18 =========
echo "****"
echo "Task 18 -> cat <(echo \"MarkerName P.value\") <(grep $gene result/nominals.txt | cut -d " " -f8,12) > tmp/metal.$gene "
echo "****"

# Define the gene corresponding to the co-localized or fine-mapped variants of interest 
cat <(echo "MarkerName P.value") <(grep $gene result/nominals.txt | cut -d " " -f8,12) > tmp/metal.$gene

