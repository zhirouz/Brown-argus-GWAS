# Brown-argus-GWAS

#【GWAS on BA – e-lab-book and scripts】

#main working folder:
/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics

#【pre-GWAS file processing】

# 1) bcftools
#-The vcf file contains 251 individuals which is fewer than in the phenotype file. 
#-need to extract the individuals with phenotypes from the vcf file. You can write a list of the names of samples in the vcf file with: bcftools query -l AA251.FINAL.MAF0.01.missing0.5perpop.vcf > names, where $VCF is your vcf file name. 
#-Compare this list of names with the list of phenotypes and makes sure you have a matching list (in the same order as well). 
 
#installing and compiling bcftools:
wget https://github.com/samtools/bcftools/releases/download/1.15/bcftools-1.15.tar.bz2
tar xvf bcftools-1.15.tar.bz2
mkdir bcftools_install
cd bcftools-1.15
./configure --prefix=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/bcftools_install
make
make install
 
export PATH=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/bcftools_install/bin:$PATH
bcftools=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/bcftools_install/bin/bcftools
 
#to get list of sample names from the vcf file:
bcftools query -l AA251.FINAL.MAF0.01.missing0.5perpop.vcf > names
 
# 2) vcftools
#-After comparing sample named from vcf and sample names from phenotypic data table, I obtained 239 samples with both genetic and phenotypic data
#-You can extract a list of individuals from a vcf file using --keep from vcftools. Have a look at their manual. http://vcftools.sourceforge.net/man_latest.html
#-You can also rename the individuals using bcftools reheader e.g. if the samples have annoying name ends
 
#installing and compiling vcftools:
git clone https://github.com/vcftools/vcftools.git
cd vcftools
./autogen.sh
./configure --prefix=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/bcftools_install
make
make install
 
export PATH=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/vcftools_install/bin:$PATH
vcftools=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/vcftools_install/bin/vcftools
 
#to generate a new vcf file with 239 filtered individuals with both vcf and phenotypic information
vcftools  --vcf AA251.FINAL.MAF0.01.missing0.5perpop.vcf  --out AA239.MAF0.01.missing0.5.22Jun --keep samples_to_keep --recode
 
#check samples in the new file:
bcftools query -l AA239.MAF0.01.missing0.5.22Jun.recode.vcf > names2
 
# 3) plink
#-Needed to prepare PLINK map and ped files, and PLINK binary files
 
#installing and compiling PLINK:
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20220603.zip
plink=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/plink2

#convert vcf file to PLINK binary files: https://www.cog-genomics.org/plink/2.0/input#pheno
$plink --vcf AA239.MAF0.01.missing0.5.22Jun.recode.vcf --allow-extra-chr --make-bed --out AA239.MAF0.01.missing0.5.binary
#output files: AA239.MAF0.01.missing0.5.binary .bed  & .bim & .fam & .cov
#to view the binary .bed file
xxd -b AA239.MAF0.01.missing0.5.binary.bed

#to make both sample IDs and family IDs as sample ID
$plink --vcf AA208.MAF0.01.missing0.5.08Jul.recode.vcf --allow-extra-chr --make-bed --double-id --out AA208.MAF0.01. 15Jul.binary

#make relatedness matrix: https://www.cog-genomics.org/plink/2.0/distance
#format A: matrix
$plink --bfile AA239.MAF0.01.missing0.5.binary --allow-extra-chr --make-rel square
#output files: GRM written to plink2.rel and IDs to plink2.rel.id

# 4) bimbam file generated using QCTools
#-BLSMM in GEMMA does not take in a separate covariate file. To incorporate and use covariates, we need BIMBAM mean genotype file format for input files instead (see GEMMA manual http://www.xzlab.org/software/GEMMAmanual.pdf)

#QCTools may be able to convert vcf to BIMBAM mean genotype files 
#https://github.com/genetics-statistics/GEMMA/issues/103
#https://enkre.net/cgi-bin/code/qctool/doc/tip/doc/html/qctool/documentation/download.html

#running QCtools
qctool=/share/apps/genomics/qctoolv2/bin/qctool_v2.0.5
$qctool -g AA208.MAF0.01.missing0.5.08Jul.recode.vcf -ofiletype bimbam_dosage -og AA208.bimbam.txt

#add rows of covariates to the bimbam mean genotype file
#first, in the excel sheet for sample data and info, code levels of all categorical covariates with numbers to make them numerical; and transform the covariate table into landscape form with covariates in rows and sample ID in columns. Following this, add 2nd and 3rd columns in this transformed table: 0 and 1 which indicates the two levels of covariates; these two columns match the two columns in the bimbam genotype file where minor and major alleles are.
#copy rows of covariates into a file in the server and then into the bimbam mean genotype file
cp AA239.bimbam.txt AA239.bimbam_with_covariates.txt
nano covariate_table
#copy content of the covariates table to the end of the bimbam file
cat covariate_table >> AA239.bimbam_with_covariates.txt
#check the first two columns and the last few rows of the new bimbam file to make sure the format is correct
awk '{print $1, $2}' AA239.bimbam_with_covariates.txt > test
tail test

#check the row number of covariates in the mean genotype file
#count total number of rows
wc -l AA239.bimbam_with_covariates.txt
#get the last few rows and see what they are
sed '61211,61214!d' AA239.bimbam_with_covariates.txt > test
awk '{print $1, $2}' test > test2
cat test2

# 5) Check relatedness for exceptional relations
#using the KING-robust kinship estimator in plink
$plink  --bfile AA239.MAF0.01.missing0.5.binary --allow-extra-chr --make-king square
#output files are plink2.king and plink2.king.id
#check kinship in R

# 6) Filter individuals again based on total number of eggs laid >=10.
#-After further filtering, I get 208 individuals with both vcf and phenotypic information available and with number of eggs >=10
vcftools  --vcf AA251.FINAL.MAF0.01.missing0.5perpop.vcf  --out AA208.MAF0.01.missing0.5.08Jul --keep samples_to_keep_AA208 --recode
$plink --vcf AA208.MAF0.01.missing0.5.08Jul.recode.vcf --allow-extra-chr --make-bed --out AA208.08Jul.binary
 

#【GWAS】

# A. GEMMA BLSMM:
#entering interactive node:
#Note from GEMMA manual: “Notice that a large memory is needed to fit BSLMM (e.g. may need 20 GB for a data set with 4000 individuals and 400,000 SNPs)”
qrsh -l tmem=10G,h_vmem=10G

#executing gemma: 
export PATH=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/software/gemma-0.98.5:$PATH
gemma=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/software/gemma-0.98.5/gemma
$gemma -h

#running BLSMM: $gemma -bfile [prefix] -bslmm [num] -o [prefix]

# [1st run]
#data: all data pooled together
#for now, I use: a standard linear BLSMM; no covariate file; automatic relatedness matrix calculation by GEMMA; and default settings
$gemma -bfile AA239.MAF0.01.missing0.5.binary -bslmm 1 -o AA239.test.26Jun

GEMMA 0.98.5 (2021-08-25) by Xiang Zhou, Pjotr Prins and team (C) 2012-2021
Reading Files ...
##number of total individuals = 239
##number of analyzed individuals = 239
##number of covariates = 1
##number of phenotypes = 1
##number of total SNPs/var        =    61210
##number of analyzed SNPs         =    31451
Start Eigen-Decomposition...
pve estimate =0.267816
se(pve) =0.22585
Calculating UtX...
initial value of h = 0.267816
initial value of rho = 1
initial value of pi = 0.000317955
initial value of |gamma| = 10
================================================== 100%    0.28
**** INFO: Done.

#P.S. the options Kristian used:

$gemma -bfile AA239.MAF0.01.missing0.5.binary -bslmm 1 -s 25,000,000 -w 2,500,000 -rpace 100 -wpace 10,000 -o AA239.MAF0.01.missing0.5.BLSMM1.25Jun

GEMMA 0.98.5 (2021-08-25) by Xiang Zhou, Pjotr Prins and team (C) 2012-2021
Reading Files ...
##number of total individuals = 239
##number of analyzed individuals = 239
##number of covariates = 1
##number of phenotypes = 1
##number of total SNPs/var        =    61210
##number of analyzed SNPs         =    31451
Start Eigen-Decomposition...
pve estimate =0.267816
se(pve) =0.22585
Calculating UtX...
initial value of h = 0.267816
initial value of rho = 1
initial value of pi = 0.000317955
initial value of |gamma| = 10
================================================== 100%    0.15
**** INFO: Done.

#output files are stored in a directory named “output”

#Basic and BLSMM options:
-w [num]: number of burn-in iterations that will be discarded (default 100,000)
-s [num]: number of sampling iterations that will be saved (default 1,000,000)
-rpace [num]: specify recording pace, record one state in every [num] steps (default 10)
-wpace [num]: specify writing pace, write values down in every [num] recorded steps (default 1000

#[QUESTIONS/OPTIONS FOR BLSMM to be considered:]
#Which BLSMM model to use?
#Number of iterations for MCMC
#Other options in MCMC chain – follow default?
#How to account for LD?


[2nd run]

#this time using bimbam input files
#with covariates; for covariates, from manual: “Finally, the mean genotype file can accommodate values other than SNP genotypes. One can use the “-notsnp” option to disable the minor allele frequency cutoff and to use any numerical values as covariates.”
##note that in the vcf file, our SNPs were already filtered with MAF 0.01 =)
#data: all data pooled together
#standard linear BLSMM
$gemma -g AA239.bimbam.txt -notsnp -p AA239.phenotype.bimbam.txt -bslmm 1 -o AA239.test.covariates.30Jun


# Understanding The output files from BLSMM:
#-log.txt: running parameters, computation time, PVE estimate and its standard error in the null linear mixed model (not the BSLMM)
#-hyp.txt: MCMC posterior samples for the hyper-parameters (h, PVE, ρ, PGE, π and |γ|), for every 10th iteration
#(interpretation of hyper-parameters is from the GEMMA manual and the GEMMA google group: https://groups.google.com/g/gemma-discussion/c/Ow3EUP_qJuU/m/X6rpKMGgAQAJ)
#----h: a parameter of the model
#----PVE: the proportion of variance in phenotypes explained by the sparse effects (Xβ) and random effects terms (u) together 🡪 total phenotypic variation explained by all SNPs (i.e. genetic variance) 🡪 When PVE is low, it indicates either: (1) heritability is relatively low, or (2) the heritability is not low, but the proportion of genetic variance captured by your available SNP genotypes is small.
#----rho: rho is the approximation to proportion of genetic variance explained by variants with major effects, and when rho is close to 0, that means that we are more in a LMM case (highly polygenic basis of the phenotype) and when rho is close to 1, it is more a BVSR case, which means that there are a few major effects loci. E.g. when rho = 0.6, it might mean that it is not well defined for either LMM or BSVR (??) probably because PVE is quite low (not much phenotypic variation explained by the set of SNPs tested and therefore, not clear which genetic architecture is behind this phenotype).
#----PGE: the proportion of genetic variance explained by the sparse effects terms (Xβ) 🡪 the proportion of PVE explained by those SNPs with significance effects on phenotypic variation
#----pi: the proportion of variants with non-zero “large” effects
#----n_gamma: the number of major-effect SNPs that contribute to PGE
#-param.txt: the posterior mean estimates for the effect size parameters (α, β|γ == 1 and γ)
#----chr: chromosome
#----rs: SNP ID
#----ps: SNP position?
#----n_miss:
#----alpha: captures the small effects that all SNPs have
#----beta: captures the additional effects of some large effect SNPs
#----gamma:
#*From manual: Notice that the beta column contains the posterior mean estimate for βi |γi == 1 rather than βi beta. Therefore, the effect size estimate for the additional effect is βiγi, and in the special case K = XXT /p, the total effect size estimate is αi + βiγi
#-bv.txt: a column of breeding value estimates uˆ. Individuals with missing phenotypes will have values of “NA”
#-gamma.txt: contains the posterior samples for the gamma, for every 10th iteration. Each row lists the SNPs included in the model in that iteration, and those SNPs are represented by their row numbers (+1) in the prefix.param.txt file.

#PIP calculation:
#-proportion of MCMC iterations that an SNP is identified as having a non-zero beta
#-for model parameters and PIP calculation, refer to (Guan & Stephens, 2011 DOI: 10.1214/11-AOAS455) 


# [A major problem in using BSLMM:]
#the phenotype data used here, i.e. phenotypic ratio (of number of eggs laid on two different host plants) is not normally distributed, and is inflated by 0 and 1.
#log transformation doesn’t work, and quantile transformation only works when transformation is done on the residuals from a linear regression of phenotypic ratio on the covariates – however, this is hard to interpret.

#solutions: 
#1) try other methods: permutation method (permGWAS); treating phenotype as categorical traits (POLMM)
#2) further process phenotype data: transformation; downsampling the inflated GM-extreme; subsampling
#3) the hybrid two-step method: analyse binary zero vs non-zero traits first, then analyse the numerical non-zero values (details see below)
#4) Separating GM-preference and RR-preference, to make both phenotype datasets right-skewed and hence earier for transformation


# B. Try treating nonnormal, zero-inflated phenotypic data as categorical data – a GWAS tool for categorical data, POLMM (a R package)

#POLMM paper: 
https://www.sciencedirect.com/science/article/pii/S0002929721001038#:~:text=In%20genome%2Dwide%20association%20studies,inappropriately%20to%20analyze%20categorical%20phenotypes.
#POLMM Github page:
https://github.com/WenjianBI/POLMM

#installing POLMM in R
library(devtools)  # author version: 2.3.0
install_github("WenjianBi/POLMM")
library(POLMM)
?POLMM  # manual of POLMM()

#first, build a regression model: POLMM_Null_Model()
#[error message:] Error in fitPOLMMcpp(t_flagSparseGRM = flagSparseGRM, t_flagGMatRatio = flagGMatRatio,  : 
  Index out of bounds: [index='scaffold_MT'].

#to fixt this: filter my vcf file with only chromosomes and not the scaffold MT
#using vcftools to filter for chromosomes http://vcftools.sourceforge.net/man_latest.html
export PATH=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/vcftools_install/bin:$PATH
vcftools=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/vcftools_install/bin/vcftools
vcftools  --vcf AA208.MAF0.01.missing0.5.08Jul.recode.vcf  --out AA208.MAF0.01.chr.only --not-chr scaffold_MT --recode

#make new plink files
plink=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/plink2
$plink --vcf AA208.MAF0.01.chr.only.recode.vcf --allow-extra-chr --make-bed --double-id --out AA208.MAF0.01.15Jul.binary

#Several ways to categroise phenotype:
#1) into bins of 0.1 intervals
#2) into bins of 0.2 intervals
#3) into GM (<0.3), RR (>0.7), and INTermediate (0.4< <0.6)


# C. Try a permutation-based test for LMM - permGWAS

#Github page for permGWAS:
https://github.com/grimmlab/permGWAS

#install
git clone https://github.com/grimmlab/permGWAS.git

#running commend in the repository
cd permGWAS
export PATH=/share/apps/python-3.8.5-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/python-3.8.5-shared/lib:$LD_LIBRARY_PATH
python3 permGWAS.py -x /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/AA208.08Jul.binary.bed -y /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/AA208.phenotype.txt --cov /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/AA208.covariates.csv --out_dir /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/permGWAS_results_18Jul

#preparing input files
#to add an additional column / combine two files by putting the columns side by side; the -d” ” means that the columns are separated by single space/tab
paste -d" " FID AA208.phenotype.txt > AA208.pheno.new.txt

#use normal plink map and ped files as genotype input
#generate map and ped files from vcf
$plink --vcf AA208.MAF0.01.missing0.5.08Jul.recode.vcf --recode ped --allow-extra-chr --out AA208.09Aug
#tidy up format
#-For what actual map and ped files should look like and what the columns mean: https://easygwas.ethz.ch/faq/view/15/
#-For trying to change chromosome IDs from SUPER_X to just number X:
https://stackoverflow.com/questions/12400217/replace-a-field-with-values-specified-in-another-file
https://stackoverflow.com/questions/46157292/replace-in-one-file-with-value-from-another-file-not-working-properly/46157609#46157609

#search for an element/text and replace it with another text
sed -i 's/SUPER_9/9/g' AA208.map
for i in {1..22}; do sed -i 's/SUPER_$i/$i/g' AA208.map; done

#use new map and ped files:
/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/new2.map
/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/new2.ped

#change column separator in the covariate.csv file
awk -v OFS=',' '{ $1 = $1; print }' AA208.covariates.csv > new.cov.csv

#replace tab with single space as the column separator
expand -t 1 new2.ped > new3.ped

#remove the underscores in sample names – by finding _ and replace it with nothing
sed -i 's/_//g' AA208.phenotype.txt
sed -i 's/_//g' new4.ped (and remember to copy map file into a new4.map file in order to keep the prefix same)
sed -i 's/_//g' new.cov.csv

#further changing the sample names to integer form
sed -i 's/BAR/01/g' AA208.phenotype.txt
sed -i 's/BCH/02/g' AA208.phenotype.txt
sed -i 's/BRO/03/g' AA208.phenotype.txt
sed -i 's/FOR/04/g' AA208.phenotype.txt
sed -i 's/HOD/05/g' AA208.phenotype.txt
sed -i 's/LYD/06/g' AA208.phenotype.txt
sed -i 's/MOF/07/g' AA208.phenotype.txt
sed -i 's/SWD/08/g' AA208.phenotype.txt
sed -i 's/WIS/09/g' AA208.phenotype.txt

sed -i 's/BAR/01/g' new4.ped
sed -i 's/BCH/02/g' new4.ped
sed -i 's/BRO/03/g' new4.ped
sed -i 's/FOR/04/g' new4.ped
sed -i 's/HOD/05/g' new4.ped
sed -i 's/LYD/06/g' new4.ped
sed -i 's/MOF/07/g' new4.ped
sed -i 's/SWD/08/g' new4.ped
sed -i 's/WIS/09/g' new4.ped

sed -i 's/BAR/01/g' new.cov.csv
sed -i 's/BCH/02/g' new.cov.csv
sed -i 's/BRO/03/g' new.cov.csv
sed -i 's/FOR/04/g' new.cov.csv
sed -i 's/HOD/05/g' new.cov.csv
sed -i 's/LYD/06/g' new.cov.csv
sed -i 's/MOF/07/g' new.cov.csv
sed -i 's/SWD/08/g' new.cov.csv
sed -i 's/WIS/09/g' new.cov.csv

#The formatted input files:
#----genotype plink: new4.map + new4.ped
#----covariate: new.cov.csv
#----phenotype: AA208.phenotype.txt
export PATH=/share/apps/python-3.8.5-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/python-3.8.5-shared/lib:$LD_LIBRARY_PATH
python3 permGWAS.py -x /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/new4.ped -y /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/AA208.phenotype.txt --cov /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/new.cov.csv --out_dir /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/permGWAS_results_18Jul

#[Error message: ModuleNotFoundError: No module named ‘xxx’]
pip install xxx

#[Error message: Exception: Phenotype phenotype_value is not in phenotype file /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/AA208.phenotype.txt]
#🡪use .txt phenotype file, with title of the phenotypic value column as phenotype_value, and with columns being separated by single space

#[Error message: Samples of genotype and phenotype do not match.]
#🡪solution: format the input files as described above to make them exactly the same as how the example files are formatted

#[Error message:]
Checked if all specified files exist. Start loading data.
Traceback (most recent call last):
  File "permGWAS.py", line 88, in <module>
    X, y, K, covs, positions, chrom, X_index = prep.load_and_prepare_data(args)
  File "/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/permGWAS/preprocessing/prepare_data.py", line 75, in load_and_prepare_data
    X = load_files.load_genotype_matrix(arguments.x, sample_index=sample_index[1])
  File "/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/permGWAS/preprocessing/load_files.py", line 129, in load_genotype_matrix
    snps.append(iupac_map[tmp[j] + tmp[j + 1]])
KeyError: '00'

#🡪emailed the developers, and developer reply: “Currently, permGWAS only supports fully imputed genotype matrices. You can impute your data for example using the software Beagle or Impute2.”
#🡪Our data contains missing genotypes for some individuals/SNPs
 
# [Imputation]
https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#home
https://faculty.washington.edu/browning/beagle/beagle.html#introduction

# [Imputation by Beagle] (adapted from Eve Taylor-Cox’s script for imputation)
 
#software needed
##bcftools
export PATH=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/bcftools_install/bin:$PATH
bcftools=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/bcftools_install/bin/bcftools
##bgzip and tabix
wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2
tar -xvjf htslib-1.15.1.tar.bz2
cd  htslib-1.15.1
./configure --prefix=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/Imputation/htslib_install
make
make install
export PATH=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/Imputation/htslib_install/bin:$PATH
bgzip=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/Imputation/htslib_install/bin/bgzip
tabix=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/Imputation/htslib_install/bin/tabix
##java
export PATH=/share/apps/java/bin:$PATH
java=/share/apps/java/bin/java
 
#Remove duplicate using bcftools
bcftools norm -d all AA208.MAF0.01.chr.only.recode.vcf -Ov -o AA208.MAF0.01.chr.only.rmDups.vcf
bgzip -c AA208.MAF0.01.chr.only.rmDups.vcf > AA208.MAF0.01.chr.only.rmDups.vcf.gz
tabix -p vcf AA208.MAF0.01.chr.only.rmDups.vcf.gz

#Change missing genotype pattern - as explained here: https://www.biostars.org/p/433624/
echo " Change missing genotype pattern"
zcat AA208.MAF0.01.chr.only.rmDups.vcf.gz | perl -pe "s/\s\.:/\t.\/.:/g" | bgzip -c > AA208.MAF0.01.chr.only.new.rmDups.vcf.gz
printf "\n"

#Make crude .map file for Beagle
echo "Make crude map fpr vcf file"
zcat AA208.MAF0.01.chr.only.new.rmDups.vcf.gz | grep -v "#" | awk -v OFMT='%f' -v OFS='\t' '{a=$2/1000000;b=a*4.2;print $1,".", b,$2;}' > crude_4.2rate.bgl.map
printf "\n\n"

#Any scaffolds that only have one SNP need to be excluded and placed into a text file
echo "Get scaffolds with only one occurence"
awk 'NR==FNR { a[$1]++ } NR!=FNR && a[$1]==1' crude_4.2rate.bgl.map crude_4.2rate.bgl.map > uniq_scaffs
#Remove the any scaffolds excluded from the chrom map
echo "Remove these scaffolds from crude map"
grep -xvf uniq_scaffs crude_4.2rate.bgl.map > new_crude_4.2rate.bgl.map
#edit eclude to be in correct format for excludemarkers
echo "uniq scaff to CHROM:POS"
awk '{print $1":"$4}' uniq_scaffs > exclude.txt
printf "\n\n"

#Run with map
echo "Run Beagle with crude map"
java -jar /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/beagle.22Jul22.46e.jar gt=AA208.MAF0.01.chr.only.new.rmDups.vcf.gz map=crude_4.2rate.bgl.map gp=TRUE out=beagle_out nthreads=64
echo "DONE"
printf "\n\n"

#Run without map
echo "Run beagle without crude map"
java -jar /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/beagle.22Jul22.46e.jar gt=AA208.MAF0.01.chr.only.new.rmDups.vcf.gz gp=TRUE out=beagle_out_nomap nthreads=64
echo "DONE"
printf "\n\n"

#[then proceed with making input files for permGWAS]
gunzip beagle_out.vcf.gz
gunzip  beagle_out_nomap.vcf.gz
plink=/SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/plink2
$plink --vcf beagle_out.vcf --recode ped --allow-extra-chr --out beagle

##tidy up format - map file
##replace chromosome names with pure numbers
sed -i 's/SUPER_//g' beagle.map
##re-format SNP ID to chr_pos by reformatting the columns:
awk '{print $1" "$1"_"$4" "$3" "$4}' beagle.map > beagle.new.map

##ped file
##replace whole columns with another column
awk '{$1 = $2; print}' beagle.ped > beagle.1.ped
awk '{$3 = $6; print}' beagle.1.ped > beagle.2.ped
awk '{$6 = $4; print}' beagle.2.ped > beagle.3.ped
awk '{$4 = $3; print}' beagle.3.ped > beagle.4.ped
awk '{$5 = $3; print}' beagle.4.ped > beagle.new.ped
##remove the underscores in sample names – by finding _ and replace it with nothing
sed -i 's/_//g' beagle.new.ped
#further changing the sample names to integer form
sed -i 's/BAR/01/g' beagle.new.ped
sed -i 's/BCH/02/g' beagle.new.ped
sed -i 's/BRO/03/g' beagle.new.ped
sed -i 's/FOR/04/g' beagle.new.ped
sed -i 's/HOD/05/g' beagle.new.ped
sed -i 's/LYD/06/g' beagle.new.ped
sed -i 's/MOF/07/g' beagle.new.ped
sed -i 's/SWD/08/g' beagle.new.ped
sed -i 's/WIS/09/g' beagle.new.ped

#The formatted input files:
#----genotype plink: beagle.mew.ped + beagle.new.map
#----covariate: new.cov.csv
#----phenotype: AA208.phenotype.txt
export PATH=/share/apps/python-3.8.5-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/python-3.8.5-shared/lib:$LD_LIBRARY_PATH
python3 permGWAS.py -x /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/Imputation/beagle.new.ped -y /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/AA208.phenotype.txt --cov /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/new.cov.csv --out_dir /SAN/ugi/LepGenomics/1_Aricia.agestis_PopGenomics/permGWAS_results

# D. Try a data transformation tool – WarpedLMM

https://github.com/PMBio/warpedLMM
#installing:
export PATH=/share/apps/python-3.8.5-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/python-3.8.5-shared/lib:$LD_LIBRARY_PATH
pip install warpedlmm

#running WarpedLMM
python -m warpedlmm AA208.08Jul.binary phenotype_file

#didn't work

# E. Separating zero and non-zero values

#methods from (Goodman et al., 2019) doi:10.1002/gepi.22162:
#“A more sophisticated approach is to adapt previously proposed variance component tests, using a hybrid binary-continuous approach. Here as a first stage, the outcome is dichotomized as zeros vs non-zeros and tested using a logistic model. Then in a second stage, the data is reduced to the non-zeros and tested using a linear model, i.e. the non-zero outcome is treated as continuous with identity link. We have implemented these methods using the previously proposed SKAT logistic and linear model methods (Lee et al., 2012; Wu et al., 2010), while adapting the methods to combine the p-values from the two stages, and term this the SKAT hurdle model.”

#Didn’t really work because the non-zero part is still non-normal


# F. Subsample individuals for easier data transformation for GEMMA

# F.a. randomly subsample 50 individuals
#*shuf -n [x] would display x randomly drawn lines from the subject, in this case the “file”

for i in {1..10}; do shuf -n 50 AA208_data_table > samples_$i; done
for i in {1..10}; do awk '{print $1}' samples_$i > name_$i; done
for i in {1..10}; do vcftools  --vcf AA208.MAF0.01.missing0.5.08Jul.recode.vcf  --out AA50.sample.$i --keep name_$i --recode
vcftools  --vcf AA208.MAF0.01.missing0.5.08Jul.recode.vcf  --out AA50.sample1 –keep name_1 --recode
#prepare phenotype file
for i in {1..10}; do awk '{print $10}' samples_$i > phenotype_$i; done
#prepare covariates
for i in {1..10}; do awk '{print $4, $6, $8}' samples_$i > covariate_$i; done
#transpose columns of covariates into rows to be added to the bimbam file
for i in {1..10}; do awk  '{printf( "%s ", $1); } END { printf( "\n" ); }' covariate_$i > cov_$i
for i in {1..10}; do awk  '{printf( "%s ", $2); } END { printf( "\n" ); }' covariate_$i >> cov_$i
for i in {1..10}; do awk  '{printf( "%s ", $3); } END { printf( "\n" ); }' covariate_$i >> cov_$i
cat cov_1 >> AA50.sample1.cov.bimbam.txt

$gemma -g AA50.sample1.cov.bimbam.txt -notsnp -p phenotype_1 -bslmm 1 -o AA50.sample1.24Jul
$gemma -g GM.02Aug.bimbam.txt -notsnp -p GM.phenoype.02Aug.txt -bslmm 1 -o GM.02Aug

# F.b. downsample the inflated GE-extreme to make the whole distribution more normal (done in R)


# G. Single plant preference and right-skewed data transformation
#**Transformation method for the right-skewed GM- and RR-only preference data
#-Here I separated GM and RR preference, using (GM/total eggs) ratio as a measurement for GM preference, and (RR/total eggs = the original phenotypic ratio) for RR preference, and excluding 0 values in both.
#-The resultant phenotypic datasets hence have extremely right-skewed distribution.
#-Log and square-root transformations didn’t work
#-End up using this transformation and it worked, resulting in smooth near-normal transformations:
https://stats.stackexchange.com/questions/513598/how-do-i-normalise-severely-right-skewed-data
#-Note that this transformation does not work on the total AA208 data when data at both extremes are included.



# Check results from GWAS
#**refer to distribution of p values and their interpretation: http://varianceexplained.org/statistics/interpreting-pvalue-histogram/

# See a tutorial of what they did with GEMMA BSLMM and their outputs:
https://github.com/visoca/popgenomworkshop-gwas_gemma

