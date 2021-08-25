Reference-based versus reference-free cell type estimation in DNA
methylation studies using human placental tissue - Script 1,
Pre-Filtering
================
Linda Dieckmann
02-07/2021

  - [preparation](#preparation)
      - [loading packages](#loading-packages)
      - [define function(s)](#define-functions)
      - [R setup](#r-setup)
      - [save packages info](#save-packages-info)
  - [loading data](#loading-data)
  - [Look for sample outliers in PC1
    methylation](#look-for-sample-outliers-in-pc1-methylation)
  - [Filtering for non-variable CpGs in placenta from EPIC
    Array](#filtering-for-non-variable-cpgs-in-placenta-from-epic-array)
      - [Identify the probes present in all studies of a
        tissue](#identify-the-probes-present-in-all-studies-of-a-tissue)
      - [Take a look at the data and start
        QC](#take-a-look-at-the-data-and-start-qc)
      - [Call Non-variable CpGs](#call-non-variable-cpgs)

Analysis referring to Edgar et al. (2017) to identify low-varying CpGs
in placenta methylation data
<https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-017-0320-z>

# preparation

## loading packages

``` r
library(RCurl)
library(ggplot2)
library(gplots)
```

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
library(RColorBrewer)
library(data.table)
library(rstatix)
```

    ## 
    ## Attaching package: 'rstatix'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

## define function(s)

``` r
intersect_all <- function(a, b, ...){
  Reduce(intersect, list(a, b, ...))
}
```

## R setup

## save packages info

``` r
writeLines(capture.output(sessionInfo()), "sessionInfo1.txt")
```

# loading data

``` r
load("Input_Data_prepared/Beta_Placenta_FullInfo_ExprSet_ITU.rda")
load("Input_Data_prepared/Beta_CVS_FullInfo_ExprSet.rda")
load("Input_Data_prepared/Beta_Placenta_FullInfo_ExprSet_PREDO.rda")
load("Input_Data_prepared/Beta_Placenta_FullInfo_ExprSet_BET.rda")
```

get Betas and Phenos

``` r
Beta_Placenta_ITU <- exprs(Beta_Placenta_FullInfo_ExprSet_ITU)
dim(Beta_Placenta_ITU)
Pheno_Placenta_ITU <- pData(Beta_Placenta_FullInfo_ExprSet_ITU)

Beta_CVS_ITU <- exprs(Beta_CVS_FullInfo_ExprSet)
dim(Beta_CVS_ITU)
Pheno_CVS_ITU <- pData(Beta_CVS_FullInfo_ExprSet)

Beta_Placenta_PREDO <- exprs(Beta_Placenta_FullInfo_ExprSet_PREDO)
dim(Beta_Placenta_PREDO)
Pheno_Placenta_PREDO <- pData(Beta_Placenta_FullInfo_ExprSet_PREDO)

Beta_Placenta_BET <- exprs(Beta_Placenta_FullInfo_ExprSet_BET)
dim(Beta_Placenta_BET)
Pheno_Placenta_BET <- pData(Beta_Placenta_FullInfo_ExprSet_BET)
```

# Look for sample outliers in PC1 methylation

``` r
# CVS
names_outliers_cvs <- Pheno_CVS_ITU[is_outlier(Pheno_CVS_ITU$PC_methB1, coef = 3), "Sample_Name"]
# Placenta ITU
names_outliers_placenta_itu <- Pheno_Placenta_ITU[is_outlier(Pheno_Placenta_ITU$PC_methB1, coef = 3), "Sample_Name"]
# Placenta PREDO
names_outliers_placenta_predo <- Pheno_Placenta_PREDO[is_outlier(Pheno_Placenta_PREDO$PC_methB1, coef = 3), "ID"]
# Placenta BET
names_outliers_placenta_bet <- Pheno_Placenta_BET[is_outlier(Pheno_Placenta_BET$PC_methB1, coef = 3), "Sample_Name"]

# only for ITU Placenta 16 outliers

# values outliers IQR 3
samples_different_placenta_itu_meth <- Pheno_Placenta_ITU[is_outlier(Pheno_Placenta_ITU$PC_methB1, coef = 3), "Sample_Name"]
```

``` r
 save(samples_different_placenta_itu_meth, file="Input_Data_prepared/samples_different_placenta_itu_meth.RData")
# sample names of those that are extreme outliers in PC1 methylation
```

# Filtering for non-variable CpGs in placenta from EPIC Array

## Identify the probes present in all studies of a tissue

Term & CVS

``` r
cpgs_in_all_all <- intersect_all(rownames(Beta_CVS_ITU), rownames(Beta_Placenta_ITU), rownames(Beta_Placenta_PREDO), rownames(Beta_Placenta_BET))

length(cpgs_in_all_all)
# 652341
```

beta values for only the cpgs in all data sets

``` r
Betas_CVS_sub_all <- Beta_CVS_ITU[rownames(Beta_CVS_ITU) %in% cpgs_in_all_all, ]
Betas_Placenta_ITU_sub_all <- Beta_Placenta_ITU[rownames(Beta_Placenta_ITU) %in% cpgs_in_all_all, ]
Betas_Placenta_PREDO_sub_all <- Beta_Placenta_PREDO[rownames(Beta_Placenta_PREDO) %in% cpgs_in_all_all, ]
Betas_Placenta_BET_sub_all <- Beta_Placenta_BET[rownames(Beta_Placenta_BET) %in% cpgs_in_all_all, ]
```

``` r
dim(Betas_CVS_sub_all)
dim(Betas_Placenta_ITU_sub_all)
dim(Betas_Placenta_PREDO_sub_all)
dim(Betas_Placenta_BET_sub_all)

# be sure rows are all the same
all(sapply(list(rownames(Betas_CVS_sub_all), rownames(Betas_Placenta_BET_sub_all), rownames(Betas_Placenta_PREDO_sub_all)), FUN = identical, rownames(Betas_Placenta_ITU_sub_all)))
```

``` r
Betas_Placenta_all <- cbind(Betas_CVS_sub_all, Betas_Placenta_ITU_sub_all, Betas_Placenta_PREDO_sub_all, Betas_Placenta_BET_sub_all)
Betas_Placenta_data_all <- as.data.frame(Betas_Placenta_all)
```

all in one data set

``` r
# trannsposed
tBetas_CVS_sub_all <- data.frame(t(Betas_CVS_sub_all))
tBetas_CVS_sub_all$cohort <- rep("CVS-ITU", 264)

tBetas_Placenta_ITU_sub_all <- data.frame(t(Betas_Placenta_ITU_sub_all))
tBetas_Placenta_ITU_sub_all$cohort <- rep("ITU", 486)

tBetas_Placenta_PREDO_sub_all <- data.frame(t(Betas_Placenta_PREDO_sub_all))
tBetas_Placenta_PREDO_sub_all$cohort <- rep("PREDO", 139)

tBetas_Placenta_BET_sub_all <- data.frame(t(Betas_Placenta_BET_sub_all))
tBetas_Placenta_BET_sub_all$cohort <- rep("BET", 137)
```

``` r
tBetas_Placenta_all <- rbind(tBetas_CVS_sub_all, tBetas_Placenta_ITU_sub_all, tBetas_Placenta_PREDO_sub_all, tBetas_Placenta_BET_sub_all)# 1026 652342
tBetas_Placenta_all_data <- as.data.frame(tBetas_Placenta_all)
tBetas_Placenta_all_data$name <- row.names(tBetas_Placenta_all) # 1026 x 652343
```

``` r
long_tBetas_Placenta_all_data <- melt(setDT(tBetas_Placenta_all_data), id.vars = c("cohort", "name"), variable.name = "CpG")
```

``` r
 save(long_tBetas_Placenta_all_data, file="Input_Data_prepared/long_tBetas_Placenta_all_data.rda")
 save(Betas_Placenta_data_all, file="Input_Data_prepared/Betas_Placenta_data_all.rda")
 save(tBetas_Placenta_all_data, file="Input_Data_prepared/tBetas_Placenta_all_data.rda")
```

## Take a look at the data and start QC

``` r
load("Input_Data_prepared/long_tBetas_Placenta_all_data.rda")
load("Input_Data_prepared/Betas_Placenta_data_all.rda")
load("Input_Data_prepared/tBetas_Placenta_all_data.rda")
```

``` r
ggplot(long_tBetas_Placenta_all_data, aes(value, group=cohort, color=cohort))+
geom_density()+theme_bw()
```

*QC I)*

``` r
# sample filtering: NA probes per sample must be < 2.5% of CpGs
na_count_sample_all <-sapply(Betas_Placenta_data_all, function(y) sum(length(which(is.na(y)))))
na_count_sample_good_all <-which(na_count_sample_all<(nrow(Betas_Placenta_data_all)*0.025))
# all samples ok
```

``` r
# probe filtering: NA of probe over samples must be < 5% of samples
na_count_probe_all <-sapply(1:nrow(Betas_Placenta_data_all), function(y) length(which(is.na(Betas_Placenta_data_all[y,]))))
na_count_probe_good_all<-which(na_count_probe_all<(ncol(Betas_Placenta_data_all)*0.05))
# all probes ok
```

*Sample-Sample Correlation*

``` r
corr_placenta_all <-cor(Betas_Placenta_data_all, use="complete.obs", method="spearman")
```

``` r
save(corr_placenta_all, file="Input_Data_prepared/corr_placenta_all.Rdata")
```

``` r
cohort_samples_all <- data.frame(tBetas_Placenta_all_data[ ,c("cohort","name")])
cohort_samples_all$cohort <- as.factor(cohort_samples_all$cohort)

cohort_samples_all <- cohort_samples_all[match(rownames(corr_placenta_all), cohort_samples_all$name), ]
identical(colnames(corr_placenta_all), cohort_samples_all$name)
```

``` r
hmcols_all <-  colorRampPalette(brewer.pal(9,"Blues"))
selcol_all <- colorRampPalette(brewer.pal(8,"Set2"))
clustcol.height_placenta_all = selcol_all(length(unique(cohort_samples_all$cohort)))
cols_placenta_all <-clustcol.height_placenta_all[as.numeric(cohort_samples_all$cohort)]
```

*heatmap*

``` r
heatmap.2(as.matrix(corr_placenta_all),Colv=FALSE,Rowv=FALSE,
          density.info="none",scale="none",ColSideColors=cols_placenta_all,RowSideColors=cols_placenta_all,
          margin=c(10,10),col=hmcols_all,symkey=F,trace="none",dendrogram="none",
          labRow=NA,
          labCol=NA, breaks=seq(0.8, 1, 0.01))

legend("topright",      
       legend = levels(cohort_samples_all$cohort),
       col = clustcol.height_placenta_all, 
       lty= 1,             
       lwd = 4,           
       cex=1
)
```

``` r
## Average sample correlation
average_sample_correlation_placenta_all<-colMeans(corr_placenta_all)
```

*Plot with PC1 methylation outliers marked*

``` r
load("../Input_Data_prepared/samples_different_placenta_itu_meth.RData")
# sample names of those that are extreme outliers in PC1 methylation
# n = 16
```

``` r
cohort_samples_all$ave_samp_cor<-average_sample_correlation_placenta_all
cohort_samples_all$outlier <- ifelse(cohort_samples_all$name %in% samples_different_placenta_itu_meth, "yes", "no")
cohort_samples_all$cohort <- gsub("CVS-ITU", "ITU-CVS", cohort_samples_all$cohort)
```

``` r
samplesample_cor_plot <- ggplot(cohort_samples_all, aes(name, ave_samp_cor, color=cohort))+
  geom_point(shape=19, size=1)+
  geom_point(data=cohort_samples_all[cohort_samples_all$outlier =="yes",], aes(x=name, y=ave_samp_cor), color="red")+
  theme_bw()+scale_color_manual(values=clustcol.height_placenta_all, name="Cohort")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y=element_text(size=10),axis.title.x = element_blank(), axis.title.y = element_text(size=10))+
  xlab("Sample")+ylab("mean sample-sample correlation (DNAm)")+
  ylim(0.7,1)+ 
  guides(colour = guide_legend(override.aes = list(size = 7, shape=15)))

samplesample_cor_plot
```

``` r
ggsave("Results/meth_sample_sample_cor.pdf",
samplesample_cor_plot, width=84, height=50, units="mm", dpi=600, scale=2, device = cairo_pdf)
```

``` r
save(samplesample_cor_plot, file="Results/RData/samplesample_cor_plot.Rdata")
```

*QC II)*

``` r
sample_names_low_cor_all <- cohort_samples_all[rstatix::is_outlier(cohort_samples_all$ave_samp_cor, coef = 3), ]
# which samples are (extreme) outliers (3*IQR) in average sample correlation?
# n = 10 

#sample_names_low_cor_all
# these are no good samples -> delete for probe variation calculation
# n = 1026 - 10 
```

``` r
# exclude the 10 samples with low average sample-sample correlation
Betas_Placenta_data_good_all <- Betas_Placenta_data_all[ ,-c(294,330,361,546,574,592,617,620,652,177)]
#  652341   1016
```

``` r
save(Betas_Placenta_data_good_all, file="../Input_Data_prepared/Betas_Placenta_data_good_all.rda")
save(corr_placenta_all, file="../Input_Data_prepared/corr_placenta_all.rda")
save(cohort_samples_all, file="../Input_Data_prepared/cohort_samples_all.rda")
```

## Call Non-variable CpGs

Calculate reference range variability (without low correlation studies,
and high NA probes)

``` r
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
```

``` r
ref_range_placenta_all<-lapply(1:nrow(Betas_Placenta_data_good_all), function(x) Variation(Betas_Placenta_data_good_all[x,]))
```

Count the non-variable CpGs

``` r
ref_range_placenta_all<-unlist(ref_range_placenta_all)
length(which(ref_range_placenta_all<0.05))
```

in EDgar et al.: placenta 101,367 (21%)

``` r
120548/652341
```

``` r
nonvariable_placenta_all <-Betas_Placenta_data_good_all[which(ref_range_placenta_all<0.05),]
```

``` r
Placenta_nonvariable_cpgs_all <- data.frame(CpG=rownames(nonvariable_placenta_all), RefRange = ref_range_placenta_all[which(ref_range_placenta_all < 0.05)])
```

``` r
write.csv(Placenta_nonvariable_cpgs_all, file="Results/Nonvariable_Placenta_CpGs_all.csv")
write.xlsx2(Placenta_nonvariable_cpgs_all, "../Results_Main/Nonvariable_Placenta_CpGs_all.xlsx", sheetName = "List", col.names = TRUE, row.names = FALSE, append = FALSE)

save(Placenta_nonvariable_cpgs_all, file="Results/RData/Placenta_nonvariable_cpgs_all.Rdata")
```

defining CpG as non-variable: threshold of 5% range in beta values (DNAm
level ranging from 0 to 1) between the 10th and 90th percentile was used
= CpGs with less than 5% reference range of beta values in were
considered non-variable
