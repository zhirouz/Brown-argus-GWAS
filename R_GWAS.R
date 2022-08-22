#read output tables
##param.txt: the posterior mean estimates for the effect size parameters (¦Á, ¦Â|¦Ã == 1 and ¦Ã)
test <- read.table(file="AAsub.sample3.29Jul.param.txt", header = TRUE)
head(test)
tail(SNP_param)
dim(SNP_param)

##"the effect size estimate for the additional effect is ¦Âi¦Ãi"
test$additional_effect <- with(test,beta*gamma)

##hyp.txt: posterior samples for the hyper-parameters (h, PVE, ¦Ñ, PGE, ¦Ð and |¦Ã|), for every 10th iteration
hyper <- read.table(file="RR.02Aug.hyp.txt", header = TRUE)
head(hyper)
dim(hyper) # n_row = 100,000 = number of MCMC runs

##gamma.txt: contains the posterior samples for the gamma, for every 10th iteration. 
##Each row lists the SNPs included in the model in that iteration, 
##and those SNPs are represented by their row numbers (+1) in the prefix.param.txt file
gamma <- read.table(file="AAsub.sample3.29Jul.gamma.txt", header = TRUE)
head(gamma)
dim(gamma) # n_row = 100,000 = number of MCMC runs


#look at SNP beta values
hist(SNP_param$beta)
summary(SNP_param$beta)
hist(SNP_param$additional_effect)
summary(SNP_param$additional_effect)

library(ggplot2)
pdf("AAsub.sample3.29Jul.betaplots.pdf",paper = "USr")
ggplot(test, aes(x=pos, y=beta,colour=CHR))+geom_point(size=0.5)
ggplot(test, aes(x=pos, y=additional_effect,colour=CHR))+geom_point(size=0.5)+facet_wrap(~CHR, nrow =2)
dev.off()

#make a vector of frequencies of SNP-row-numbers appearing in the gamma table
##the SNP row numbers are row numbers in the SNP parameter table

dim(test) # total number of rows = total number of SNPs = 38666
SNP_row_number <- c()
for (i in 1:33821)
  {
    SNP_row_number[i]<-sum(gamma== i)
  }
##check output vector
length(SNP_row_number) # this number should equal to 31451
sum(SNP_row_number) # this number should equal to sum(hyper$n_gamma)
sum(hyper$n_gamma)

#add the SNP frequency vector to the SNP table
test$MCMC_frequency <- SNP_row_number
test$MCMC_prop <- with(test,MCMC_frequency/100000)
max(SNP_param$MCMC_frequency)
max(SNP_param$MCMC_prop)
head(test)

##this MCMC_prop should be the simplest form of PIP?

#save the new SNP_param table with SNP MCMC frequencies
write.table(test,"AAsub.sample3.29Jul.txt",row.names=FALSE)
test<-read.table(file="AA164.binarytrait.txt", header = TRUE)
head(test)
tail(test)

install.packages("tidyr")
library(dplyr)
library(tidyr)
?separate()
test <- test %>% separate(col=rs, into=c('na','chr', 'pos'),sep=":")

CHR <- c()
CHR[grep("SUPER_9",test$chr)] <- "Chr 9"
CHR[grep("SUPER_1",test$chr)] <- "Chr 1"
CHR[grep("SUPER_2",test$chr)] <- "Chr 2"
CHR[grep("SUPER_3",test$chr)] <- "Chr 3"
CHR[grep("SUPER_4",test$chr)] <- "Chr 4"
CHR[grep("SUPER_5",test$chr)] <- "Chr 5"
CHR[grep("SUPER_6",test$chr)] <- "Chr 6"
CHR[grep("SUPER_7",test$chr)] <- "Chr 7"
CHR[grep("SUPER_8",test$chr)] <- "Chr 8"
CHR[grep("SUPER_10",test$chr)] <- "Chr 10"
CHR[grep("SUPER_11",test$chr)] <- "Chr 11"
CHR[grep("SUPER_12",test$chr)] <- "Chr 12"
CHR[grep("SUPER_13",test$chr)] <- "Chr 13"
CHR[grep("SUPER_14",test$chr)] <- "Chr 14"
CHR[grep("SUPER_15",test$chr)] <- "Chr 15"
CHR[grep("SUPER_16",test$chr)] <- "Chr 16"
CHR[grep("SUPER_17",test$chr)] <- "Chr 17"
CHR[grep("SUPER_18",test$chr)] <- "Chr 18"
CHR[grep("SUPER_19",test$chr)] <- "Chr 19"
CHR[grep("SUPER_20",test$chr)] <- "Chr 20"
CHR[grep("SUPER_21",test$chr)] <- "Chr 21"
CHR[grep("SUPER_22",test$chr)] <- "Chr 22"
CHR[grep("SUPER_Z",test$chr)] <- "Chr Z"
CHR[grep("SCA",test$chr)] <- "Scaffold"
CHR[grep(".year",test$na)] <- "year"
CHR[grep(".site.type",test$na)] <- "site.type"
CHR[grep(".site.plant",test$na)] <- "site.plant"
length(CHR)
test$CHR <- CHR



  
#plot PIP across SNPs
library(ggplot2)
pdf("AAsub.sample3.29Jul.SNPplots.pdf",paper = "USr")
#pooled plot coloured by chr
ggplot(test, aes(x=pos, y=MCMC_prop,colour=CHR))+geom_point(size=0.7)
#plots separated by chr
ggplot(test, aes(x=pos, y=MCMC_prop))+geom_point(size=0.7)+facet_wrap(~CHR, nrow =1)
#pooled plot with top points (PIP>0.025?) labelled
ggplot(test, aes(x=pos, y=MCMC_prop,colour=CHR,label=as.character(pos)))+geom_point(size=0.8)+geom_text(aes(label=ifelse(MCMC_prop>0.5,as.character(pos),'')),hjust=0, vjust=0)
#plots separated by chr with top points (PIP>0.025?) labelled
ggplot(test, aes(x=pos, y=MCMC_prop,label=as.character(pos)))+geom_point(size=0.6)+facet_wrap(~CHR, nrow =2)+geom_text(aes(label=ifelse(MCMC_prop>0.5,as.character(pos),'')),hjust=0, vjust=0)
dev.off()


#deal with hyperparameters
head(hyper)
pdf("AA239.test.hyperplots.pdf",paper = "USr")
ggplot(hyper, aes(x=pve))+geom_density()+geom_vline(aes(xintercept=mean(pve)),color="blue", linetype="dashed", size=1)
ggplot(hyper, aes(x=pge))+geom_density()+geom_vline(aes(xintercept=mean(pge)),color="blue", linetype="dashed", size=1)
dev.off()

summary(hyper$pve)
summary(hyper$pge)


#check relatedness matrix
##Read matrix files
RM <- read.table(file="plink2.rel", header = FALSE)
RM <- as.matrix(RM)
head(RM)
KING <- read.table(file="plink2.king", header = FALSE)
KING <- as.matrix(KING)
##check values
summary(KING)
##plots
help(heatmap)
heatmap(RM,main="heat.map.relatedness.matrix")
heatmap(KING,main="heat.map.KING.kinship",col=heat.colors(3))
legend(x="right",legend=c("min", "med", "max"),fill=heat.colors(3))

##check king kinship estimator values
value_cat <- c("negative","0","<=0.125","<=0.177","<=0.354","<0.5","0.5")
count <- c()
count[1]<-sum(KING<0)
count[2]<-sum(KING==0)
count[3]<-sum(KING>0 & KING<=0.125)
count[4]<-sum(KING>0.125 & KING<=0.177)
count[5]<-sum(KING>0.177 & KING<=0.354)
count[6]<-sum(KING>0.354 & KING<0.5)
count[7]<-sum(KING==0.5)
relatedness_values <- data.frame(value_cat, count)
relatedness_values

library(dplyr)
KING <- as.data.frame(KING)
KING_positive <- KING %>% filter(KING>0)
KING_positive

#check phenotype data
library("readxl")
##just phenotype table
pheno <- read_excel("GWAS BA data table.xlsx", sheet = 2)
head(pheno)
ggplot(pheno, aes(x=phenotypic.ratio))+geom_density()+geom_vline(aes(xintercept=mean(phenotypic.ratio)),color="blue", linetype="dashed", size=1)
##phenotyope + covariates
PnC <- read_excel("GWAS BA data table.xlsx", sheet = 4)
head(PnC)

pheno <- read.table("phenotype_1")
pheno <- as.data.frame(pheno)
head(pheno)
ggplot(pheno, aes(x=V1))+geom_density()+geom_vline(aes(xintercept=mean(V1)),color="blue", linetype="dashed", size=1)


##phenotypic data transformation
###Code for Inverse normal transformation from (Yang et al., 2012):
####INT per se:
y <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))

##1) Directly transforming phenotyoic data:
pheno$trans.ratio <- qnorm((rank(pheno$phenotypic.ratio)-0.5)/sum(!is.na(pheno$phenotypic.ratio)))
head(pheno)
ggplot(pheno, aes(x=trans.ratio))+geom_density()+geom_vline(aes(xintercept=mean(trans.ratio)),color="blue", linetype="dashed", size=1)

##2)regress phenotype on covariates, get residuals, and transform residuals:
attach(PnC)
LM <- lm(phenotypic.ratio ~ site_C + year_C + site.type_C + site.plant_C)
summary(LM)
PnC$res <- residuals(LM)
PnC$trans.res <- qnorm((rank(PnC$res,na.last="keep")-0.5)/sum(!is.na(PnC$res)))
head(PnC$trans.res)
ggplot(PnC, aes(x=trans.res))+geom_density()+geom_vline(aes(xintercept=mean(trans.res)),color="blue", linetype="dashed", size=1)

LM <- lm(phenotypic.ratio ~ site.plant_C)
PnC$res <- residuals(LM)
PnC$trans.res <- qnorm((rank(PnC$res,na.last="keep")-0.5)/sum(!is.na(PnC$res)))
ggplot(PnC, aes(x=trans.res))+geom_density()+geom_vline(aes(xintercept=mean(trans.res)),color="blue", linetype="dashed", size=1)

LM <- lm(phenotypic.ratio ~ site.type_C)
PnC$res <- residuals(LM)
PnC$trans.res <- qnorm((rank(PnC$res,na.last="keep")-0.5)/sum(!is.na(PnC$res)))
ggplot(PnC, aes(x=trans.res))+geom_density()+geom_vline(aes(xintercept=mean(trans.res)),color="blue", linetype="dashed", size=1)


##the distribution of the transformed residuals is standard normal with mean = 0!!!!
library("writexl")
write_xlsx(pheno,"table of transformed data.xlsx")
write_xlsx(PnC,"table of transformed residual data.xlsx")

###The R package for INT
install.packages("RNOmni")
library("RNOmni")
z <- RankNorm(pheno$phenotypic.ratio)
z <- RankNorm(PnC$res)
head(z)
plot(density(z))


# distribution of number of eggs
head(AA208)
ggplot(AA208, aes(x=total.RR))+geom_density()+geom_vline(aes(xintercept=mean(total.RR)),color="blue", linetype="dashed", size=1)


#make data files for subsampling GEMMA & separation of GM vs RR preference
AA208 <- read_excel("AA208 with full categorical traits.xlsx", sheet = 1)
head(AA208)
GM <- AA208 %>% filter(GM_pref!=0)
ggplot(GM, aes(x=GM_pref))+geom_density()+geom_vline(aes(xintercept=mean(GM_pref)),color="blue", linetype="dashed", size=1) + ggtitle("distribution of GM preference")
GM$trans_GM_pref <- qnorm(ecdf(GM$GM_pref)(GM$GM_pref))
ggplot(GM, aes(x=trans_GM_pref))+geom_density()+geom_vline(aes(xintercept=mean(trans_GM_pref)),color="blue", linetype="dashed", size=1)+ ggtitle("distribution of transformed GM preference")

RR <- AA208 %>% filter(phenotypic.ratio!=0)
ggplot(RR, aes(x=phenotypic.ratio))+geom_density()+geom_vline(aes(xintercept=mean(phenotypic.ratio)),color="blue", linetype="dashed", size=1)+ ggtitle("distribution of RR preference")
RR$trans_RR_pref <- qnorm(ecdf(RR$phenotypic.ratio)(RR$phenotypic.ratio))
ggplot(RR, aes(x=trans_RR_pref))+geom_density()+geom_vline(aes(xintercept=mean(trans_RR_pref)),color="blue", linetype="dashed", size=1)+ ggtitle("distribution of transformed RR preference")

pdf("distributions of phenotypic & transformed phenotypic data.pdf")
dev.off()


library("writexl")
write_xlsx(GM,"GM preference.xlsx")
write_xlsx(RR,"RR preference.xlsx")
GM <- read_excel("GM preference.xlsx")
RR <- read_excel("RR preference.xlsx")
head(RR)
library(ggplot2)

length(which(AA208$pheno.group.0.2=="1"))
length(which(AA208$pheno.group.0.2=="2"))
length(which(AA208$pheno.group.0.2=="3"))
length(which(AA208$pheno.group.0.2=="4"))
length(which(AA208$pheno.group.0.2=="5"))

AA208_GMextreme <- subset(AA208,AA208$pheno.group.0.2=="1")
write_xlsx(AA208_GMextreme,"AA208_GMextreme.xlsx")
AA208_GMextreme <- subset(AA208,AA208$pheno.group.0.2!="1")
write_xlsx(AA208_GMextreme,"AA208_without_GMextreme.xlsx")



#he binary-continuous hybrid method - separating zeros and non-zero values
##check distribution of non-zero values and transform
AA159 <- read_excel("GWAS_hybrid_binary-continuous.xlsx", sheet = 3)
head(AA159)
dim(AA159)
library(ggplot2)
ggplot(AA159, aes(x=phenotypic.ratio))+geom_density()+geom_vline(aes(xintercept=mean(phenotypic.ratio)),color="blue", linetype="dashed", size=1)
AA159$trans.ratio <- log(AA159$phenotypic.ratio)
ggplot(AA159, aes(x=trans.ratio))+geom_density()+geom_vline(aes(xintercept=mean(trans.ratio)),color="blue", linetype="dashed", size=1)

##the package by (Goodman et al., 2019)
library("devtools")
install.packages("CompQuadForm")
library("BiocManager")
BiocManager::install("EmpiricalBrownsMethod")
devtools::install_github("MatthewOGoodman/ZeroInflatedVCtest")
library(ZeroInflatedVCtest)

#GWAS on categorical data - POLMM
##installing POLMM in R
install.packages("devtools")
library("devtools")  # author version: 2.3.0
install_github("WenjianBi/POLMM")
library(POLMM)

#fit regression model
?POLMM_Null_Model() 
help(POLMM)

#test for associations
?POLMM
?POLMM.plink()

#read data table
library("readxl")
AA208 <- read_excel("GWAS BA data table.xlsx", sheet = 2)
AA208 <- read_excel("AA208 with full categorical traits.xlsx")
head(AA208)
dim(AA208)
head(AA208$phenotypic.ratio)

#bin phenotypic data into categories _ bin size 0.1

library(dplyr)

AA208 <- AA208 %>% mutate(pheno.group = case_when(phenotypic.ratio >= 0  & phenotypic.ratio < 0.1 ~ '1',
                                                  phenotypic.ratio >= 0.1  & phenotypic.ratio < 0.2 ~ '2',
                                                  phenotypic.ratio >= 0.2  & phenotypic.ratio < 0.3 ~ '3',
                                                  phenotypic.ratio >= 0.3  & phenotypic.ratio < 0.4 ~ '4',
                                                  phenotypic.ratio >= 0.4  & phenotypic.ratio < 0.5 ~ '5',
                                                  phenotypic.ratio >= 0.5  & phenotypic.ratio < 0.6 ~ '6',
                                                  phenotypic.ratio >= 0.6  & phenotypic.ratio < 0.7 ~ '7',
                                                  phenotypic.ratio >= 0.7  & phenotypic.ratio < 0.8 ~ '8',
                                                  phenotypic.ratio >= 0.8  & phenotypic.ratio < 0.9 ~ '9',
                                                  phenotypic.ratio >= 0.9  & phenotypic.ratio <= 1 ~ '10')) # end function

length(which(AA208$pheno.group=="1"))
length(which(AA208$pheno.group=="2"))
length(which(AA208$pheno.group=="3"))
length(which(AA208$pheno.group=="4"))
length(which(AA208$pheno.group=="5"))
length(which(AA208$pheno.group=="6"))
length(which(AA208$pheno.group=="7"))
length(which(AA208$pheno.group=="8"))
length(which(AA208$pheno.group=="9"))
length(which(AA208$pheno.group=="10"))

length(which(AA208$site.type=="new" & AA208$pheno.group=="1"))
length(which(AA208$site.type=="new" & AA208$pheno.group=="10"))
length(which(AA208$site.type=="old" & AA208$pheno.group=="1"))
length(which(AA208$site.type=="old" & AA208$pheno.group=="10"))

length(which(AA208$year=="2013" & AA208$pheno.group=="1"))
length(which(AA208$year=="2013" & AA208$pheno.group=="10"))
length(which(AA208$year=="2014" & AA208$pheno.group=="1"))
length(which(AA208$year=="2014" & AA208$pheno.group=="10"))
length(which(AA208$year=="2014"))

# bin into only 2 (what Kristian did) or 3 (plus intermediate) or 0.2-bin-size categories 

AA208 <- AA208 %>% mutate(pheno.cat = case_when(phenotypic.ratio >= 0  & phenotypic.ratio <= 0.3 ~ 'GM',
                                                  phenotypic.ratio >= 0.7  & phenotypic.ratio <= 1 ~ 'RR',
                                                  phenotypic.ratio >= 0.35  & phenotypic.ratio <=0.65 ~ 'INT')) # end function
head(AA208$pheno.cat)
tail(AA208$pheno.cat)
length(which(AA208$pheno.cat=="GM"))
length(which(AA208$pheno.cat=="RR"))
length(which(AA208$pheno.cat=="INT"))

AA208 <- AA208 %>% mutate(pheno.binary = case_when(phenotypic.ratio >= 0  & phenotypic.ratio < 0.3 ~ 'GM',
                                                  phenotypic.ratio > 0.7  & phenotypic.ratio <= 1 ~ 'RR')) # end function
head(AA208$pheno.binary)
length(which(AA208$pheno.binary=="GM"))
length(which(AA208$pheno.binary=="RR"))

AA208 <- AA208 %>% mutate(pheno.group.0.2 = case_when(phenotypic.ratio >= 0  & phenotypic.ratio < 0.2 ~ '1',
                                                  phenotypic.ratio >= 0.2  & phenotypic.ratio < 0.4 ~ '2',
                                                  phenotypic.ratio >= 0.4  & phenotypic.ratio < 0.6 ~ '3',
                                                  phenotypic.ratio >= 0.6  & phenotypic.ratio < 0.8 ~ '4',
                                                  phenotypic.ratio >= 0.8  & phenotypic.ratio <= 1 ~ '5')) # end function

head(AA208$pheno.group.0.2)
length(which(AA208$pheno.group.0.2=="1"))
length(which(AA208$pheno.group.0.2=="2"))
length(which(AA208$pheno.group.0.2=="3"))
length(which(AA208$pheno.group.0.2=="4"))
length(which(AA208$pheno.group.0.2=="5"))
#min # = 14

#write these tables
head(AA208)
library("writexl")
write_xlsx(AA208,"AA208 with full categorical traits.xlsx")

AA208_bi <- AA208 %>% filter(AA208$pheno.binary!="NA")
head(AA208_bi$pheno.binary)
dim(AA208_bi)
write_xlsx(AA208_bi,"AA164 with only binary traits.xlsx")

AA208_cat <- AA208 %>% filter(AA208$pheno.cat!="NA")
head(AA208_cat$pheno.cat)
dim(AA208_cat)
write_xlsx(AA208_cat,"AA199 with only 3-category traits.xlsx")

# separate 2013 and 2014 years

AA2013 <- AA208 %>% filter(AA208$year_C=="2013")
head(AA2013)
AA2014 <- AA208 %>% filter(AA208$year_C=="2014")
head(AA2014)
write_xlsx(AA2013,"AA208 with only year 2013.xlsx")
write_xlsx(AA2014,"AA208 with only year 2014.xlsx")

AA208_bi_2013 <- AA208_bi %>% filter(AA208_bi$year_C=="2013")
AA208_bi_2014 <- AA208_bi %>% filter(AA208_bi$year_C=="2014")
write_xlsx(AA208_bi_2013,"AA164 binary with only year 2013.xlsx")
write_xlsx(AA208_bi_2014,"AA164 binary with only year 2014.xlsx")

pdf("distribution of phenotypic ratio in 2013 vs 2014.pdf")
ggplot(AA2013, aes(x=phenotypic.ratio))+geom_density()+geom_vline(aes(xintercept=mean(phenotypic.ratio)),color="blue", linetype="dashed", size=1)
ggplot(AA2014, aes(x=phenotypic.ratio))+geom_density()+geom_vline(aes(xintercept=mean(phenotypic.ratio)),color="blue", linetype="dashed", size=1)
dev.off()

#subsampling GM-end categories
head(AA208)
length(which(AA208$pheno.group.0.2=="1"))
length(which(AA208$pheno.group.0.2=="2"))
length(which(AA208$pheno.group.0.2=="3"))
length(which(AA208$pheno.group.0.2=="4"))
length(which(AA208$pheno.group.0.2=="5"))
#median number of individuals/0.2 category ==24
length(which(AA208$pheno.cat=="GM"))
length(which(AA208$pheno.cat=="RR"))
length(which(AA208$pheno.cat=="INT"))
#GM=134; RR=30; INT=35

#subsampling
?subset()
?sample()
AA208_sub0.2cat <- subset(AA208,AA208$pheno.group.0.2=="1")
AA208_sub0.2cat <- AA208_sub0.2cat[sample(nrow(AA208_sub0.2cat), 24), ]
nrow(AA208_sub0.2cat)
AA208_sub0.2cat <- bind_rows(AA208_sub0.2cat,AA208[which(AA208$pheno.group.0.2!="1"),])
length(which(AA208_sub0.2cat$pheno.group.0.2=="1"))
length(which(AA208_sub0.2cat$pheno.group.0.2=="2"))
length(which(AA208_sub0.2cat$pheno.group.0.2=="3"))
length(which(AA208_sub0.2cat$pheno.group.0.2=="4"))
length(which(AA208_sub0.2cat$pheno.group.0.2=="5"))
ggplot(AA208_sub0.2cat, aes(x=phenotypic.ratio))+geom_density()+geom_vline(aes(xintercept=mean(phenotypic.ratio)),color="blue", linetype="dashed", size=1)
write_xlsx(AA208_sub0.2cat,"AA208 with 0.2 cat subsampled.xlsx")

head(AA2014)
length(which(AA2014$pheno.group.0.2=="1"))
length(which(AA2014$pheno.group.0.2=="2"))
length(which(AA2014$pheno.group.0.2=="3"))
length(which(AA2014$pheno.group.0.2=="4"))
length(which(AA2014$pheno.group.0.2=="5")) #median = 6
AA2014_sub0.2cat <- subset(AA2014,AA2014$pheno.group.0.2=="1")
AA2014_sub0.2cat <- AA2014_sub0.2cat[sample(nrow(AA2014_sub0.2cat), 6), ]
nrow(AA2014_sub0.2cat)
AA2014_sub0.2cat <- bind_rows(AA2014_sub0.2cat,AA2014[which(AA2014$pheno.group.0.2!="1"),])
length(which(AA2014_sub0.2cat$pheno.group.0.2=="1"))
length(which(AA2014_sub0.2cat$pheno.group.0.2=="2"))
length(which(AA2014_sub0.2cat$pheno.group.0.2=="3"))
length(which(AA2014_sub0.2cat$pheno.group.0.2=="4"))
length(which(AA2014_sub0.2cat$pheno.group.0.2=="5"))
ggplot(AA2014_sub0.2cat, aes(x=phenotypic.ratio))+geom_density()+geom_vline(aes(xintercept=mean(phenotypic.ratio)),color="blue", linetype="dashed", size=1)
write_xlsx(AA2014_sub0.2cat,"AA2014 with 0.2 cat subsampled.xlsx")

#subsampling all categories to 10 samples each

AA208_cat1 <- subset(AA208,AA208$pheno.group.0.2=="1")
AA208_cat2 <- subset(AA208,AA208$pheno.group.0.2=="2")
AA208_cat3 <- subset(AA208,AA208$pheno.group.0.2=="3")
AA208_cat4 <- subset(AA208,AA208$pheno.group.0.2=="4")
AA208_cat5 <- subset(AA208,AA208$pheno.group.0.2=="5")

AA208_cat1sub <- AA208_cat1[sample(nrow(AA208_cat1), 10), ]
AA208_cat2sub <- AA208_cat2[sample(nrow(AA208_cat2), 10), ]
AA208_cat3sub <- AA208_cat3[sample(nrow(AA208_cat3), 10), ]
AA208_cat4sub <- AA208_cat4[sample(nrow(AA208_cat4), 10), ]
AA208_cat5sub <- AA208_cat5[sample(nrow(AA208_cat5), 10), ]

AA208_sub <- bind_rows(AA208_cat1sub,AA208_cat2sub,AA208_cat3sub,AA208_cat4sub,AA208_cat5sub)
ggplot(AA208_sub, aes(x=phenotypic.ratio))+geom_density()+geom_vline(aes(xintercept=mean(phenotypic.ratio)),color="blue", linetype="dashed", size=1)
head(AA208_sub)
write_xlsx(AA208_sub,"subsampled AA208 10 samples per cat.xlsx")

#POLMM
library(POLMM)

#build null model
NM = POLMM_Null_Model(as.factor(pheno.group)~as.factor(year)+site.type+site.plant,
                           data=AA208, PlinkFile = "AA208.MAF0.01.15Jul.binary",subjData = AA208$Sample.ID)

NM_cat = POLMM_Null_Model(as.factor(pheno.cat)~as.factor(year)+site.type+site.plant,
                      data=AA208_cat, PlinkFile = "AA208.MAF0.01.15Jul.binary",subjData = AA208_cat$Sample.ID)

NM_bi = POLMM_Null_Model(as.factor(pheno.binary)~as.factor(year)+site.type+site.plant,
                          data=AA208_bi, PlinkFile = "AA208.MAF0.01.15Jul.binary",subjData = AA208_bi$Sample.ID)

NM_2013 = POLMM_Null_Model(as.factor(pheno.group.0.2)~site.type_N+site.plant_N,
                         data=AA2013, PlinkFile = "AA208.year2013.23Jul.binary",subjData = AA2013$Sample.ID)

NM_2013cat = POLMM_Null_Model(as.factor(pheno.cat)~site.type_N+site.plant_N,
                           data=AA2013, PlinkFile = "AA208.year2013.23Jul.binary",subjData = AA2013$Sample.ID)

NM_2014 = POLMM_Null_Model(as.factor(pheno.group)~site.type+site.plant,
                           data=AA2014, PlinkFile = "AA208.MAF0.01.15Jul.binary",subjData = AA2014$Sample.ID)

NM_2014sub = POLMM_Null_Model(as.factor(pheno.group)~site.type+site.plant,
                           data=AA2014_sub0.2cat, PlinkFile = "AA208.MAF0.01.15Jul.binary",subjData = AA2014_sub0.2cat$Sample.ID)

NM_208sub = POLMM_Null_Model(as.factor(pheno.group)~as.factor(year)+site.type+site.plant,
                              data=AA208_sub0.2cat, PlinkFile = "AA208.MAF0.01.15Jul.binary",subjData = AA208_sub0.2cat$Sample.ID)

NM_sub = POLMM_Null_Model(as.factor(pheno.group.0.2)~year_N+site.type_N+site.plant_N,
                          data=AA208_sub, PlinkFile = "AA208.MAF0.01.15Jul.binary",subjData = AA208_sub$Sample.ID)


#association testing
POLMM.plink(objNull=NM,PlinkFile = "AA208.MAF0.01.15Jul.binary",output.file = "AA208.categorical.MAF0.01",minMAF = 0.01)
POLMM.plink(objNull=NM,PlinkFile = "AA208.MAF0.01.15Jul.binary",output.file = "AA208.categorical.MAF0.05",minMAF = 0.05)
POLMM.plink(objNull=NM,PlinkFile = "AA208.MAF0.01.15Jul.binary",output.file = "AA208.categorical.MAF0.1",minMAF = 0.1)

POLMM.plink(objNull=NM_bi,PlinkFile = "AA208.MAF0.01.15Jul.binary",output.file = "A.binary.MAF0.01",minMAF = 0.01)
POLMM.plink(objNull=NM_cat,PlinkFile = "AA208.MAF0.01.15Jul.binary",output.file = "B.3-category.MAF0.01",minMAF = 0.01)
POLMM.plink(objNull=NM_2013,PlinkFile = "AA208.year2013.23Jul.binary",output.file = "C.AA2013.0.2cat.MAF0.01",minMAF = 0.01)
POLMM.plink(objNull=NM_2013cat,PlinkFile = "AA208.year2013.23Jul.binary",output.file = "C.AA2013.3cat.MAF0.01",minMAF = 0.01)
POLMM.plink(objNull=NM_2014,PlinkFile = "AA208.MAF0.01.15Jul.binary",output.file = "D.AA2014.categorical.MAF0.01",minMAF = 0.01)

POLMM.plink(objNull=NM_2014sub,PlinkFile = "AA208.MAF0.01.15Jul.binary",output.file = "F.AA2014.subsampled.MAF0.01",minMAF = 0.01)
POLMM.plink(objNull=NM_208sub,PlinkFile = "AA208.MAF0.01.15Jul.binary",output.file = "E.AA208.subsampled.MAF0.01",minMAF = 0.01)

POLMM.plink(objNull=NM_sub,PlinkFile = "AA208.MAF0.01.15Jul.binary",output.file = "G.AA208sub.n50.MAF0.01",minMAF = 0.01)

#read output file

outPOLMM = read.table("AA208.categorical.MAF0.01", header = T)
head(outPOLMM)
tail(outPOLMM) #61209 rows
out2 = read.table("AA208.categorical.MAF0.05", header = T)
head(out2)
tail(out2)
out3 = read.table("AA208.categorical.MAF0.1", header = T)
head(out3)

out_cat = read.table("B.3-category.MAF0.01", header = T)
head(out_cat)
out_2014 = read.table("D.AA2014.categorical.MAF0.01", header = T)
head(out_2014)

out_208sub = read.table("E.AA208.subsampled.MAF0.01", header = T)
head(out_208sub)

out_208sub = read.table("G.AA208sub.n50.MAF0.01", header = T)
head(out_208sub)

out_2013 = read.table("C.AA2013.0.2cat.MAF0.01",header=T)
out_2013cat = read.table("C.AA2013.3cat.MAF0.01",header=T)

#plot output

##get SNP position info
SNP<- read.table("SNP_pos.txt",header=T)
head(SNP)
tail(SNP) #61209 rows
##add SNP pos to output tables
outPOLMM$SNP_POS <- SNP$POS
out2$SNP_POS <- SNP$POS
out3$SNP_POS <- SNP$POS
out_cat$SNP_POS <- SNP$POS
out_2014$SNP_POS <- SNP$POS
out_2013cat$SNP_POS <- SNP$POS

##filter out2 by NA (MAF 0.05)
library(dplyr)
out2_new <- out2 %>% filter(out2$MAF>0.05)
head(out2_new)

out3_new <- out3 %>% filter(out3$MAF>0.1)
head(out3_new)

##plot p values and beta values by position and chromosome
library(ggplot2)
#MAF0.01
pdf("AA208.POLMM.MAF0.01.pdf",paper = "USr")
#pooled plot coloured by chr
ggplot(outPOLMM, aes(x=SNP_POS, y=beta,colour=chr))+geom_point(size=0.6)
ggplot(outPOLMM, aes(x=SNP_POS, y=pval.norm,colour=chr))+geom_point(size=0.6)
ggplot(outPOLMM, aes(x=SNP_POS, y=pval.spa,colour=chr))+geom_point(size=0.6)
#plots separated by chr
ggplot(outPOLMM, aes(x=SNP_POS, y=beta))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(outPOLMM, aes(x=SNP_POS, y=pval.norm))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(outPOLMM, aes(x=pval.norm))+geom_density()+geom_vline(aes(xintercept=mean(pval.norm)),color="blue", linetype="dashed", size=1)
dev.off()

#MAF0.05
pdf("AA208.POLMM.MAF0.05.pdf",paper = "USr")
#pooled plot coloured by chr
ggplot(out2_new, aes(x=SNP_POS, y=beta,colour=chr))+geom_point(size=0.6)
ggplot(out2_new, aes(x=SNP_POS, y=pval.norm,colour=chr))+geom_point(size=0.6)
#plots separated by chr
ggplot(out2_new, aes(x=SNP_POS, y=beta))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(out2_new, aes(x=SNP_POS, y=pval.norm))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(out2_new, aes(x=pval.norm))+geom_density()+geom_vline(aes(xintercept=mean(pval.norm)),color="blue", linetype="dashed", size=1)
dev.off()

#MAF0.1
pdf("AA208.POLMM.MAF0.1.pdf",paper = "USr")
#pooled plot coloured by chr
ggplot(out3_new, aes(x=SNP_POS, y=beta,colour=chr))+geom_point(size=0.6)
ggplot(out3_new, aes(x=SNP_POS, y=pval.norm,colour=chr))+geom_point(size=0.6)
#plots separated by chr
ggplot(out3_new, aes(x=SNP_POS, y=beta))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(out3_new, aes(x=SNP_POS, y=pval.norm))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(out3_new, aes(x=pval.norm))+geom_density()+geom_vline(aes(xintercept=mean(pval.norm)),color="blue", linetype="dashed", size=1)
dev.off()

#3-category trait MAF0.01
pdf("3-category.MAF0.01.pdf",paper = "USr")
#pooled plot coloured by chr
ggplot(out_cat, aes(x=SNP_POS, y=beta,colour=chr))+geom_point(size=0.6)
ggplot(out_cat, aes(x=SNP_POS, y=pval.norm,colour=chr))+geom_point(size=0.6)
#plots separated by chr
ggplot(out_cat, aes(x=SNP_POS, y=beta))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(out_cat, aes(x=SNP_POS, y=pval.norm))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(out_cat, aes(x=pval.norm))+geom_density()+geom_vline(aes(xintercept=mean(pval.norm)),color="blue", linetype="dashed", size=1)
dev.off()

#only year 2014 MAF0.01
pdf("year 2014.MAF0.01.pdf",paper = "USr")
#pooled plot coloured by chr
ggplot(out_2014, aes(x=SNP_POS, y=beta,colour=chr))+geom_point(size=0.6)
ggplot(out_2014, aes(x=SNP_POS, y=pval.norm,colour=chr))+geom_point(size=0.6)
#plots separated by chr
ggplot(out_2014, aes(x=SNP_POS, y=beta))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(out_2014, aes(x=SNP_POS, y=pval.norm))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(out_2014, aes(x=pval.norm))+geom_density()+geom_vline(aes(xintercept=mean(pval.norm)),color="blue", linetype="dashed", size=1)
dev.off()

#only year 2013 MAF0.01
pdf("year 2013.3cat.pdf",paper = "USr")
#pooled plot coloured by chr
ggplot(out_2013cat, aes(x=SNP_POS, y=beta,colour=chr))+geom_point(size=0.6)
ggplot(out_2013cat, aes(x=SNP_POS, y=pval.norm,colour=chr))+geom_point(size=0.6)
#plots separated by chr
ggplot(out_2013cat, aes(x=SNP_POS, y=beta))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(out_2013cat, aes(x=SNP_POS, y=pval.norm))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(out_2013cat, aes(x=pval.norm))+geom_density()+geom_vline(aes(xintercept=mean(pval.norm)),color="blue", linetype="dashed", size=1)
dev.off()


#subsampled AA208 - 10 per cat
pdf("AA208sub.n50.MAF0.01.pdf",paper = "USr")
#pooled plot coloured by chr
ggplot(out_208sub, aes(x=SNP_POS, y=beta,colour=chr))+geom_point(size=0.6)
ggplot(out_208sub, aes(x=SNP_POS, y=pval.norm,colour=chr))+geom_point(size=0.6)
#plots separated by chr
ggplot(out_208sub, aes(x=SNP_POS, y=beta))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(out_208sub, aes(x=SNP_POS, y=pval.norm))+geom_point(size=0.6)+facet_wrap(~chr, nrow =1)
ggplot(out_208sub, aes(x=pval.norm))+geom_density()+geom_vline(aes(xintercept=mean(pval.norm)),color="blue", linetype="dashed", size=1)
dev.off()


#check p value distribution using QQ plots
length(outPOLMM$pval.norm) #61209
expected <- seq(1/61209, 1, length.out=61209)
head(expected)
qqplot(-log(expected), -log(outPOLMM$pval.norm), pch = 20, xlab = "-log (expected)", ylab = "-log (observed)")

length(out2_new$pval.norm) #35462
expected2 <- seq(1/35462, 1, length.out=35462)
head(expected2)
qqplot(-log(expected2), -log(out2_new$pval.norm), pch = 20, xlab = "-log (expected)", ylab = "-log (observed)")

length(out3_new$pval.norm) #24091
expected3 <- seq(1/24091, 1, length.out=24091)
head(expected3)
qqplot(-log(expected3), -log(out3_new$pval.norm), pch = 20, xlab = "-log (expected)", ylab = "-log (observed)")



#permGWAS
pval <- read.csv("p_values_phenotype_value.csv")
head(pval)
dim(pval)
max(pval$p_value)

library(ggplot2)
pdf("permGWAS.pval.plots.pdf",paper = "USr")
#pooled plot coloured by chr
ggplot(pval, aes(x=POS, y=p_value,colour=CHR))+geom_point(size=0.7)
#plots separated by chr
ggplot(pval, aes(x=POS, y=p_value))+geom_point(size=0.7)+facet_wrap(~CHR, nrow =1)
#density distribution of p-values
ggplot(pval, aes(x=p_value))+geom_density()+geom_vline(aes(xintercept=mean(p_value)),color="blue", linetype="dashed", size=1)
dev.off()

pdf("permGWAS.effectsize.plots.pdf",paper = "USr")
ggplot(pval, aes(x=POS, y=effect_size,colour=CHR))+geom_point(size=0.5)
ggplot(pval, aes(x=POS, y=effect_size,colour=CHR))+geom_point(size=0.5)+facet_wrap(~CHR, nrow =2)
dev.off()



# [Supplementary codes - not part of my GWAS work]

## -------- ABANDONED CODE:quantile normalisation to a randomly generated reference normal dsitribution
ref.dis <- rnorm(nrow(pheno), mean(pheno$phenotypic.ratio), sd(pheno$phenotypic.ratio))
head(ref.dis)
head(pheno$phenotypic.ratio)
RANK <- rank(pheno$phenotypic.ratio)
head(RANK)
RANK2 <- rank(ref.dis)
head(RANK2)
df <- data.frame(ref.dis,RANK2)
head(df)
df <- df[order(df$RANK2), ]

head(pheno)
pheno$RANK <- RANK
pheno$QT <- NA
is.integer(pheno$RANK[1])
for (i in 1:239)
{
  if (is.integer(pheno$RANK[i])=="TRUE") {
    pheno$QT[i] <- df$ref.dis[which(df$RANK2==pheno$RANK[i])]
  } else {
    pheno$QT[i] <- mean(df$ref.dis[which(df$RANK2==35)],df$ref.dis[which(df$RANK2==36)])
  }
}
## --------till here: ABANDONED CODE ABOVE -------------


## -------- set up after updating Rtools
Sys.which("make")
writeLines('PATH="${RTOOLS42_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
