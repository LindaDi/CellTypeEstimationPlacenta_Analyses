Reliability of a novel approach for reference-based cell type estimation
in human placental DNA methylation studies - Script 2, Main_Analyses
================
Linda Dieckmann
02-07/2021

-   [preparation](#preparation)
    -   [loading packages](#loading-packages)
    -   [define function(s)](#define-functions)
    -   [R setup](#r-setup)
    -   [save packages info](#save-packages-info)
-   [loading data](#loading-data)
-   [Descriptives](#descriptives)
    -   [methylation](#methylation)
        -   [1st PC](#1st-pc)
    -   [Cell Types](#cell-types)
        -   [plots reference-based cell
            types](#plots-reference-based-cell-types)
    -   [Phenos](#phenos)
-   [Correlation reference-free and reference-based estimated cell
    types](#correlation-reference-free-and-reference-based-estimated-cell-types)
-   [Correlation reference-free estimated cell types and
    phenotypes](#correlation-reference-free-estimated-cell-types-and-phenotypes)
    -   [CVS](#cvs)
    -   [Placenta ITU](#placenta-itu)
    -   [Placenta PREDO](#placenta-predo)
    -   [Placenta BET](#placenta-bet)
-   [Check if methylation can be predicted using cell types
    (using CV)](#check-if-methylation-can-be-predicted-using-cell-types-using-cv)
    -   [CVS](#cvs-1)
        -   [predict PC1 methylation](#predict-pc1-methylation)
    -   [Placenta ITU](#placenta-itu-1)
        -   [predict PC1 methylation](#predict-pc1-methylation-1)
    -   [Placenta PREDO](#placenta-predo-1)
        -   [predict PC1 methylation](#predict-pc1-methylation-2)
    -   [Placenta BET](#placenta-bet-1)
        -   [predict PC1 methylation](#predict-pc1-methylation-3)
        -   [Make outlier plot for BET CV
            Model](#make-outlier-plot-for-bet-cv-model)
    -   [arrange plots together](#arrange-plots-together)
-   [Check how methylation of single CpGs can be predicted using cell
    types](#check-how-methylation-of-single-cpgs-can-be-predicted-using-cell-types)
    -   [CVS](#cvs-2)
    -   [Placenta ITU](#placenta-itu-2)
    -   [Placenta PREDO](#placenta-predo-2)
    -   [Placenta BET](#placenta-bet-2)
    -   [arrange Model figures](#arrange-model-figures)
-   [CpGs most influenced by cell types
    (reference-based)](#cpgs-most-influenced-by-cell-types-reference-based)
    -   [extract CpGs](#extract-cpgs)
    -   [genes of all CpGs that overlap between data
        sets](#genes-of-all-cpgs-that-overlap-between-data-sets)
    -   [extract CpGs with R2 > 30% predicted by reference-based cell
        types that are in all data
        sets](#extract-cpgs-with-r2--30-predicted-by-reference-based-cell-types-that-are-in-all-data-sets)
    -   [Tissue-specific Gene
        Enrichment](#tissue-specific-gene-enrichment)
        -   [genes with R2 > .30 that overlap between the data sets
            against all
            genes](#genes-with-r2--30-that-overlap-between-the-data-sets-against-all-genes)
    -   [Placenta cell enrichment](#placenta-cell-enrichment)
-   [CpGs most influenced by cell types
    (reference-free)](#cpgs-most-influenced-by-cell-types-reference-free)
    -   [extract CpGs](#extract-cpgs-1)
    -   [extract CpGs with R2 > 30% that are in all data sets
        (ref-free)](#extract-cpgs-with-r2--30-that-are-in-all-data-sets-ref-free)
    -   [Tissue-specific Gene
        Enrichment](#tissue-specific-gene-enrichment-1)
        -   [genes with R2 > .30 that overlap between the data sets
            against all
            genes](#genes-with-r2--30-that-overlap-between-the-data-sets-against-all-genes-1)
    -   [Placenta cell enrichment](#placenta-cell-enrichment-1)
-   [save Tissue Enrichment Plot](#save-tissue-enrichment-plot)
-   [save Cell Enrichment Plot](#save-cell-enrichment-plot)
-   [Plot Cell types (RPC,
    reference-based)](#plot-cell-types-rpc-reference-based)
-   [Compare data sets in RPC cell type
    proportions](#compare-data-sets-in-rpc-cell-type-proportions)
    -   [Placenta: CVS vs. term Placenta in
        ITU](#placenta-cvs-vs-term-placenta-in-itu)
    -   [Term Placentas](#term-placentas)
-   [Phenotype relationships](#phenotype-relationships)
    -   [CVS (ITU)](#cvs-itu)
    -   [Placenta (ITU)](#placenta-itu-3)
    -   [Placenta (PREDO)](#placenta-predo-3)
    -   [Placenta (BET)](#placenta-bet-3)
    -   [arrange plots together](#arrange-plots-together-1)

{#top}

# preparation

## loading packages

``` r
library(plyr)
library(dplyr)
library(devtools)
library(ggplot2)
library(corrplot)
library(gridExtra)
library(Hmisc)
library(tidyr)
library(viridis)
library(psych)
library(ggpubr)
library(ggdendro)
library(reshape2)
library(planet)
library(GGally)
library(xvalglms)
library(cowplot)
library(rcompanion)
library(ggcorrplot)
library(TissueEnrich)
library(readxl)
library(tidyr)
library(factoextra)
library(npmv)
library(plyr)
library(psych)
library(matrixStats)
library(rcompanion)
library(planet)
library(Biobase)
library(rstatix)
library(gridGraphics)
library(rsq)
library(xlsx)
library(scales)
```

## define function(s)

``` r
scaleFUN <- function(x) sprintf("%.2f", x)
```

``` r
boxplot_fun = function(x, id, dat) {
  dat[ ,c(x, id)] %>%   
  pivot_longer(cols=-id) %>% 
  rename_at(1, ~"id") %>%
  ggplot(aes(x = .data[["name"]], y = .data[["value"]])) +
    geom_boxplot() +
    theme_bw()+
    labs(x="", y="")+
    theme(axis.text.x = element_text(size=6, angle=30, hjust=1), axis.title.x=element_text(size=6), axis.title.y = element_text(size=6), axis.text.y = element_text(size=6))
}
```

## R setup

## save packages info

``` r
writeLines(capture.output(sessionInfo()), "sessionInfo2.txt")
```

# loading data

``` r
# data with filtered CpGs & without outliers for ITU placenta
load("Input_Data_prepared/Beta_Placenta_FullInfo_ExprSet_ITU_reduced_filtered.rda")
load("Input_Data_prepared/Beta_Placenta_FullInfo_ExprSet_ITU_filtered.rda")
load("Input_Data_prepared/Beta_CVS_FullInfo_ExprSet_filtered.rda")
load("Input_Data_prepared/Beta_Placenta_FullInfo_ExprSet_PREDO_filtered.rda")
load("Input_Data_prepared/Beta_Placenta_FullInfo_ExprSet_BET_filtered.rda")
```

``` r
Beta_Placenta_ITU_reduced_filtered <-exprs(Beta_Placenta_FullInfo_ExprSet_ITU_reduced_filtered)
dim(Beta_Placenta_ITU_reduced_filtered)
Pheno_Placenta_ITU_reduced_filtered <-pData(Beta_Placenta_FullInfo_ExprSet_ITU_reduced_filtered)

Beta_Placenta_ITU_filtered <-exprs(Beta_Placenta_FullInfo_ExprSet_ITU_filtered)
dim(Beta_Placenta_ITU_filtered)
Pheno_Placenta_ITU_filtered <-pData(Beta_Placenta_FullInfo_ExprSet_ITU_filtered)

Beta_CVS_ITU_filtered <-exprs(Beta_CVS_FullInfo_ExprSet_filtered)
dim(Beta_CVS_ITU_filtered)
Pheno_CVS_ITU_filtered <-pData(Beta_CVS_FullInfo_ExprSet_filtered)

Beta_Placenta_PREDO_filtered <-exprs(Beta_Placenta_FullInfo_ExprSet_PREDO_filtered)
dim(Beta_Placenta_PREDO_filtered)
Pheno_Placenta_PREDO_filtered <-pData(Beta_Placenta_FullInfo_ExprSet_PREDO_filtered)

Beta_Placenta_BET_filtered <-exprs(Beta_Placenta_FullInfo_ExprSet_BET_filtered)
dim(Beta_Placenta_BET_filtered)
Pheno_Placenta_BET_filtered <-pData(Beta_Placenta_FullInfo_ExprSet_BET_filtered)
```

``` r
names_refbased_cells <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC", "Syncytiotrophoblast")
```

``` r
data("plColors") # get color from publication from Yuan et al., to make it easier comparable
```

``` r
# load 600 reference CpGs from Yuan et al.
data("plCellCpGsFirst")
cpgs_1yuan <- rownames(plCellCpGsFirst)
data("plCellCpGsThird")
cpgs_3yuan <- rownames(plCellCpGsThird)
```

load 16 samples different in ITU placenta (names)

``` r
load("Input_Data_prepared/samples_different_placenta_itu_meth.RData")
# sample names of those that are extreme outliers in PC1 methylation
```

load list of non-variable CpGs

``` r
load("Results/RData/Placenta_nonvariable_cpgs_all.Rdata")
```

# Descriptives

## methylation

### 1st PC

``` r
load("Input_Data_prepared/pc_m_cvs_itu_filtered.RData")
load("Input_Data_prepared/pc_m_placenta_itu_reduced_filtered.RData")
load("Input_Data_prepared/pc_m_placenta_itu_filtered.RData")
load("Input_Data_prepared/pc_m_placenta_predo_filtered.RData")
load("Input_Data_prepared/pc_m_placenta_BET_filtered.RData")
```

``` r
PC_plot_cvs <- fviz_eig(pc_m_cvs_itu_filtered, barfill = "grey", barcolor = "grey", main="")+ 
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title=element_text(size=10),
        axis.text = element_text(size = 10), 
        plot.margin= margin(t = 5.5, r = 10.5, b = 5.5, l = 5.5))+ 
  labs(title="CVS (ITU)", size=10)+
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits=c(0,20), labels=scaleFUN) 

PC_plot_placenta_itu <- fviz_eig(pc_m_placenta_itu_reduced_filtered, barfill = "grey",barcolor = "grey", main="")+ 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 10), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), 
        plot.margin= margin(t = 5.5, r = 10.5, b = 5.5, l = 20.5))+
  labs(title="Placenta (ITU)", size=10)+
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits=c(0,20), labels=scaleFUN) 

PC_plot_placenta_predo <- fviz_eig(pc_m_placenta_predo_filtered, barfill = "grey",barcolor = "grey", main="")+ 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 10), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), 
        plot.margin= margin(t = 5.5, r = 15.5, b = 5.5, l = 20.5))+ 
  labs(title="Placenta (PREDO)", size=10)+
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits=c(0,20), labels=scaleFUN)

PC_plot_placenta_BET <-fviz_eig(pc_m_placenta_BET_filtered, barfill = "grey",barcolor = "grey", main="")+ 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 10), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), 
        plot.margin= margin(t = 5.5, r = 5.5, b = 5.5, l = 20.5))+ 
  labs(title="Placenta (BET)", size=10)+
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits=c(0,20), labels=scaleFUN) 

# See how it looks including our outliers in ITU placenta
pc1_placenta_itu_allpersons <- Pheno_Placenta_ITU_filtered$PC_methB1
group_placenta_itu_allpersons <- replicate(486, "Placenta (ITU)")
outlier_placenta_itu_allpersons <- ifelse(Pheno_Placenta_ITU_filtered$Sample_Name %in% Pheno_Placenta_ITU_reduced_filtered$Sample_Name, "no", "yes")
PC1_Placenta_ITU_allpersons <- data.frame(pc1_placenta_itu_allpersons, group_placenta_itu_allpersons, outlier_placenta_itu_allpersons, stringsAsFactors = FALSE)
colnames(PC1_Placenta_ITU_allpersons) <- c("PC1", "group", "outlier")

#PC1_Placenta_ITU_allpersons
PC1_Placenta_ITU_allpersons_plot <- ggplot(PC1_Placenta_ITU_allpersons, aes(x=group, y=PC1, fill=outlier))+
  geom_boxplot(fill='#A4A4A4', color="black")+
  geom_point(data=PC1_Placenta_ITU_allpersons[PC1_Placenta_ITU_allpersons$outlier =="yes",], aes(x=group, y=PC1, color="red"))+
  theme_minimal()+
  labs(y="PC1 (DNAm)", title="Placenta (ITU)", x="")+
  theme(axis.text.x = element_blank(), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_blank(), legend.position = "none", plot.title = element_text(size=10, hjust = 0.5)) 
```

``` r
M_Scree <- ggarrange(
          PC_plot_cvs +
           theme(legend.position="none", plot.title=element_text(size=10)), 
          PC_plot_placenta_itu +
               theme(legend.position="none",axis.title.y = element_blank(), axis.text.y = element_blank(), plot.title=element_text(size=10)), #, plot.margin = margin(r = 0.2, l = 0.2)
          #, axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank()
          PC_plot_placenta_predo +
               theme(legend.position="none",axis.title.y = element_blank(), axis.text.y = element_blank(), plot.title=element_text(size=10)), #, plot.margin = margin(l = 0.2)
          PC_plot_placenta_BET +
               theme(legend.position="none",axis.title.y = element_blank(), axis.text.y = element_blank(), plot.title=element_text(size=10)), #, plot.margin = margin(l = 0.2)
          nrow = 1,
          labels = c("a", "b", "c","d"),
          align = "h", widths = c(1.3, 1, 1, 1))

annotate_figure(M_Scree,bottom = text_grob("Principal Component",size = 12))
```

PC1 CVS: 10.9% (2: 4.4%, 3: 2.8%) PC1 Placenta ITU: 12.7% (2: 3.9%, 3:
2.6%) PC1 Placenta PREDO: 11.6% (2: 3.6%, 3 2.7%) PC1 Placenta BET: 8.5%
(2: 3.7%, 3: 3.2%)

``` r
ggsave("Results/Methylation_PC_Screeplots.pdf",
annotate_figure(M_Scree,bottom = text_grob("Principal Component",size = 12)), width=84, height=50, units="mm", dpi=600, scale=2, device = cairo_pdf)
```

<!-- see which samples are outliers in PC1 -->
<!-- ```{r} -->
<!-- # use rstatix package -->
<!-- # values outliers IQR 1.5 -->
<!-- #is_outlier(Pheno_Placenta_ITU$PC_methB1, coef = 1.5) -->
<!-- samples_IQR_outliers_placenta_itu_meth <- Pheno_Placenta_ITU[is_outlier(Pheno_Placenta_ITU$PC_methB1, coef = 1.5), "Sample_Name"] -->
<!-- # values outliers IQR 3 -->
<!-- #is_outlier(Pheno_Placenta_ITU$PC_methB1, coef = 3) -->
<!-- samples_different_placenta_itu_meth <- Pheno_Placenta_ITU[is_outlier(Pheno_Placenta_ITU$PC_methB1, coef = 3), "Sample_Name"] -->
<!-- ``` -->
<!-- ```{r} -->
<!-- save(samples_different_placenta_itu_meth,file="Input_Data_prepared/samples_different_placenta_itu_meth.RData") -->
<!-- # sample names of those that are extreme outliers in PC1 methylation -->
<!-- ``` -->
<!-- these outliers are excluded in the _reduced data -->

*take a look at the outlier samples from ITU placenta*

``` r
load("Results/RData/samplesample_cor_plot.Rdata")
```

``` r
PC1_Placenta_ITU_allpersons_plot
samplesample_cor_plot
```

``` r
ggsave("Results/Methylation_PC1_PlacentaITU_allpersons_boxplot.pdf",
PC1_Placenta_ITU_allpersons_plot, width=84, height=50, units="mm", dpi=600, scale=2, device = cairo_pdf)
```

[to the top](#top)

## Cell Types

``` r
desc_placenta_itu_rpc <- psych::describe(Pheno_Placenta_ITU_reduced_filtered[,names_refbased_cells])
desc_cvs_itu_rpc <- psych::describe(Pheno_CVS_ITU_filtered[,names_refbased_cells])
desc_placenta_predo_rpc <- psych::describe(Pheno_Placenta_PREDO_filtered[,names_refbased_cells])
desc_placenta_BET_rpc <- psych::describe(Pheno_Placenta_BET_filtered[,names_refbased_cells])

desc_cvs_itu_rpc
desc_placenta_itu_rpc
desc_placenta_predo_rpc
desc_placenta_BET_rpc
```

``` r
desc_placenta_itu_rf <- psych::describe(Pheno_Placenta_ITU_reduced_filtered[,c("C1", "C2", "C3", "C4", "C5","C6","C7", "C8")])
desc_cvs_itu_rf <- psych::describe(Pheno_CVS_ITU_filtered[,c("C1", "C2", "C3", "C4", "C5")])
desc_placenta_predo_rf <- psych::describe(Pheno_Placenta_PREDO_filtered[,c("C1", "C2")])
desc_placenta_BET_rf <- psych::describe(Pheno_Placenta_BET_filtered[,c("C1", "C2", "C3")])

desc_cvs_itu_rf
desc_placenta_itu_rf
desc_placenta_predo_rf
desc_placenta_BET_rf
```

### plots reference-based cell types

barplots

``` r
barplot_cells_placenta_itu <- 
  ggplot(gather(Pheno_Placenta_ITU_reduced_filtered[,names_refbased_cells]), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x', ncol=6)+
  scale_x_continuous(labels=scaleFUN)+
  theme_minimal()

barplot_cells_cvs_itu <- 
  ggplot(gather(Pheno_CVS_ITU_filtered[,names_refbased_cells]), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x', ncol=6)+
  scale_x_continuous(labels=scaleFUN)+
  theme_minimal()

barplot_cells_placenta_predo <- 
  ggplot(gather(Pheno_Placenta_PREDO_filtered[,names_refbased_cells]), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x', ncol=6)+
  scale_x_continuous(labels=scaleFUN)+
  theme_minimal()

barplot_cells_placenta_BET <- 
  ggplot(gather(Pheno_Placenta_BET_filtered[,names_refbased_cells]), aes(value)) +
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x', ncol=6)+
  scale_x_continuous(labels=scaleFUN)+
  theme_minimal()
```

``` r
barplot_cells_cvs_itu

barplot_cells_placenta_itu

barplot_cells_placenta_predo

barplot_cells_placenta_BET
```

save barplots

``` r
png(file="Results/barplot_rpc_placenta.png", width= 6000, height=2000, res=600)
barplot_cells_placenta_itu
dev.off()

png(file="Results/barplot_rpc_cvs.png", width= 6000, height=2000, res=600)
barplot_cells_cvs_itu
dev.off()

png(file="Results/barplot_rpc_placenta.png", width= 6000, height=2000, res=600)
barplot_cells_placenta_predo
dev.off()

png(file="Results/barplot_rpc_placenta.png", width= 6000, height=2000, res=600)
barplot_cells_placenta_BET
dev.off()
```

boxplots

``` r
ggplot(gather(Pheno_CVS_ITU_filtered[,names_refbased_cells]), aes(x=key, y= value)) + 
  geom_boxplot(fill='#A4A4A4', color="black") +
   ylab("proportion (%)")+
  theme_minimal() +
  theme(axis.text.x = element_text(size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_blank()) 

ggplot(gather(Pheno_Placenta_ITU_reduced_filtered[,names_refbased_cells]), aes(x=key, y= value)) + 
  geom_boxplot(fill='#A4A4A4', color="black") +
   ylab("proportion (%)")+
  theme_minimal() +
  theme(axis.text.x = element_text(size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_blank()) 

ggplot(gather(Pheno_Placenta_PREDO_filtered[,names_refbased_cells]), aes(x=key, y= value)) + 
  geom_boxplot(fill='#A4A4A4', color="black") +
   ylab("proportion (%)")+
  theme_minimal() +
  theme(axis.text.x = element_text(size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_blank()) 

ggplot(gather(Pheno_Placenta_BET_filtered[,names_refbased_cells]), aes(x=key, y= value)) + 
  geom_boxplot(fill='#A4A4A4', color="black") +
   ylab("proportion (%)")+
  theme_minimal() +
  theme(axis.text.x = element_text(size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_blank()) 
```

outlier plots

``` r
Pheno_Placenta_ITU_filtered$outlier <- ifelse(Pheno_Placenta_ITU_filtered$Sample_Name %in% samples_different_placenta_itu_meth, "yes", "no")
```

``` r
Placenta_refbasedCells_allpersons_plot <- ggplot(gather(Pheno_Placenta_ITU_filtered[,names_refbased_cells]), aes(x=key, y= value)) + 
  geom_boxplot(fill='#A4A4A4', color="black") +
  geom_point(data=gather(Pheno_Placenta_ITU_filtered[Pheno_Placenta_ITU_filtered$outlier =="yes",names_refbased_cells]), aes(x=key, y=value), color="red")+
  ylab("cell types proportion (%)")+
  labs(title="Placenta (ITU)")+
  theme_minimal() +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_blank(), legend.position="none", plot.title = element_text(size=10, hjust = 0.5)) 
   
Placenta_refbasedCells_allpersons_plot
```

group and save the plots that show how placenta methylation outliers are
different

``` r
ITU_Placenta_persons_different_plot <- ggarrange(
          PC1_Placenta_ITU_allpersons_plot,
          samplesample_cor_plot,
          Placenta_refbasedCells_allpersons_plot,
          nrow = 1,
          labels = c("a", "b", "c"),
          align = "hv")

ITU_Placenta_persons_different_plot
```

``` r
ggsave("Results/Methylation_ITU_Placenta_Samples_different.pdf",
ITU_Placenta_persons_different_plot, width=174, height=70, units="mm", dpi=600, scale=2, device = cairo_pdf)
```

## Phenos

``` r
desc_placenta_itu_pheno <- psych::describe(Pheno_Placenta_ITU_reduced_filtered[,c("Gestational_Age_Weeks", "PC1_ethnicity", "PC2_ethnicity")])
desc_cvs_itu_pheno <- psych::describe(Pheno_CVS_ITU_filtered[,c("gestage_at_CVS_weeks", "PC1_ethnicity", "PC2_ethnicity")])
desc_placenta_predo_pheno <- psych::describe(Pheno_Placenta_PREDO_filtered[,c("Gestational_Age", "PC1_ethnicity", "PC2_ethnicity")])
desc_placenta_BET_pheno <- psych::describe(Pheno_Placenta_BET_filtered[,c("gestage_weeks", "PC1_ethnicity", "PC2_ethnicity", "PC3_ethnicity", "PC4_ethnicity")])

desc_cvs_itu_pheno
prop.table(table(Pheno_CVS_ITU_filtered$Child_Sex))
desc_placenta_itu_pheno
prop.table(table(Pheno_Placenta_ITU_reduced_filtered$Child_Sex))
desc_placenta_predo_pheno
prop.table(table(Pheno_Placenta_PREDO_filtered$Child_Sex))
desc_placenta_BET_pheno
prop.table(table(Pheno_Placenta_BET_filtered$Sex))
prop.table(table(Pheno_Placenta_BET_filtered$group))
```

# Correlation reference-free and reference-based estimated cell types

*Fig. 1*

To make publication-ready plot for corrplot need to adjust the function
to have significance stars:

change the place_points function within the corrplot function:

``` r
#trace(corrplot, edit=TRUE)
```

Then replace on line 443

``` r
# place_points = function(sig.locs, point) {
#   text(pos.pNew[, 1][sig.locs], pos.pNew[, 2][sig.locs], 
#        labels = point, col = pch.col, cex = pch.cex, 
#        lwd = 2)
```

with:

``` r
# place_points = function(sig.locs, point) {
#       text(pos.pNew[, 1][sig.locs]+0.40, (pos.pNew[, 2][sig.locs]), 
#            labels = point, col = pch.col, cex = pch.cex, 
#            lwd = 2)
```

here 0.3 0.4 0.2 0.25 is added to the x-position, depending on the plot

data prepared for correlation

``` r
reffree_rpc_placenta_itu <- Pheno_Placenta_ITU_reduced_filtered[,c("C1","C2","C3","C4","C5","C6","C7","C8", names_refbased_cells)]
reffree_rpc_cvs_itu <- Pheno_CVS_ITU_filtered[,c("C1","C2","C3","C4","C5", names_refbased_cells)]
reffree_rpc_placenta_predo <- Pheno_Placenta_PREDO_filtered[,c("C1","C2",names_refbased_cells)]
reffree_rpc_placenta_BET <- Pheno_Placenta_BET_filtered[,c("C1","C2","C3",names_refbased_cells)]

phenor_reffree_placenta_itu <- reffree_rpc_placenta_itu[,c("C1","C2","C3","C4","C5","C6","C7","C8")]
phenor_reffree_cvs_itu <- reffree_rpc_cvs_itu[,c("C1","C2","C3","C4","C5")]
phenor_reffree_placenta_predo <- reffree_rpc_placenta_predo[,c("C1","C2")]
phenor_reffree_placenta_BET <- reffree_rpc_placenta_BET[,c("C1","C2","C3")]

phenor_rpc_placenta_itu <- reffree_rpc_placenta_itu[,c(names_refbased_cells)]
phenor_rpc_cvs_itu <- reffree_rpc_cvs_itu[,c(names_refbased_cells)]
phenor_rpc_placenta_predo <- reffree_rpc_placenta_predo[,c(names_refbased_cells[ names_refbased_cells != "Hofbauer"])]
phenor_rpc_placenta_BET <- reffree_rpc_placenta_BET[,c(names_refbased_cells)]
```

``` r
colnames(phenor_rpc_cvs_itu) <- c("Trophoblasts        ", "Stromal       ", "Hofbauer       ", "Endothelial       ", "nRBC       ", "Synctio-\n trophoblasts        ")
colnames(phenor_rpc_placenta_itu) <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC", "Synctio-\n trophoblasts")
colnames(phenor_rpc_placenta_predo) <- c("Trophoblasts      ", "Stromal      ", "Endothelial      ", "nRBC      ", "Synctio-\n trophoblasts        \n")
colnames(phenor_rpc_placenta_BET) <- c("Trophoblasts       ", "Stromal       ", "Hofbauer       ", "Endothelial       ", "nRBC       ", "Synctio-\n trophoblasts               ")
```

CVS

``` r
r_cor_reffree_rpc_cvs_itu <- psych::corr.test(phenor_reffree_cvs_itu, phenor_rpc_cvs_itu, method="spearman", adjust = "bonferroni")
```

``` r
par(xpd=TRUE)
plot_r_cor_reffree_rpc_cvs_itu <- corrplot(r_cor_reffree_rpc_cvs_itu$r,
     method="color", number.digits = 1,
     mar=c(1,1,1,1),
     addCoef.col = "black",
     p.mat = r_cor_reffree_rpc_cvs_itu$p,
     insig = "label_sig",
     sig.level = c(.001, .01),
     pch.cex = 0.7,
     pch.col = "black",
     tl.col="black",
     tl.cex=0.75,
     tl.srt = 45,
     number.cex=0.7,
     cl.pos='n')

grid.echo()
P1_corr_refbased_reffree <- grid.grab()
```

``` r
matrix.colors <- getGrob(P1_corr_refbased_reffree, gPath("square"), grep = TRUE)[["gp"]][["fill"]]
P1_corr_refbased_reffree <- editGrob(P1_corr_refbased_reffree,
               gPath("square"), grep = TRUE,
               gp = gpar(col = NA,
                         fill = NA))

# apply the saved colours to the underlying matrix grob
P1_corr_refbased_reffree<- editGrob(P1_corr_refbased_reffree,
               gPath("symbols-rect-1"), grep = TRUE,
               gp = gpar(fill = matrix.colors))

# convert the background fill from white to transparent, while we are at it
P1_corr_refbased_reffree <- editGrob(P1_corr_refbased_reffree,
               gPath("background"), grep = TRUE,
               gp = gpar(fill = NA))
```

``` r
r_cor_reffree_rpc_placenta_itu <- psych::corr.test(phenor_reffree_placenta_itu, phenor_rpc_placenta_itu, method="spearman", adjust = "bonferroni")
```

``` r
plot_r_cor_reffree_rpc_placenta_itu <- corrplot(r_cor_reffree_rpc_placenta_itu$r,
     method="color", number.digits=1,
     mar=c(1,0,1,0),
     addCoef.col = "black",
     p.mat = r_cor_reffree_rpc_placenta_itu$p,
     insig = "label_sig",
     sig.level = c(.001, .01),
     pch.cex = 0.7,
     pch.col = "black",
     tl.col="black",
     tl.cex=0.75,
     tl.srt = 45,
     number.cex=0.7,
     cl.pos='n')

grid.echo()
P2_corr_refbased_reffree <- grid.grab()
```

``` r
matrix.colors <- getGrob(P2_corr_refbased_reffree, gPath("square"), grep = TRUE)[["gp"]][["fill"]]
P2_corr_refbased_reffree <- editGrob(P2_corr_refbased_reffree,
               gPath("square"), grep = TRUE,
               gp = gpar(col = NA,
                         fill = NA))

P2_corr_refbased_reffree<- editGrob(P2_corr_refbased_reffree,
               gPath("symbols-rect-1"), grep = TRUE,
               gp = gpar(fill = matrix.colors))

P2_corr_refbased_reffree <- editGrob(P2_corr_refbased_reffree,
               gPath("background"), grep = TRUE,
               gp = gpar(fill = NA))
```

``` r
r_cor_reffree_rpc_placenta_predo <- psych::corr.test(phenor_reffree_placenta_predo, phenor_rpc_placenta_predo, method="spearman", adjust = "bonferroni")
```

``` r
plot_r_cor_reffree_rpc_placenta_predo <- corrplot(r_cor_reffree_rpc_placenta_predo$r,
     method="color", number.digits=1,
     mar=c(1,6,1,6),
     addCoef.col = "black",
     p.mat = r_cor_reffree_rpc_placenta_predo$p,
     insig = "label_sig",
     sig.level = c(.001, .01),
     pch.cex = 0.75,
     pch.col = "black",
     tl.col="black",
     tl.cex=0.7,
     tl.srt = 45,
     number.cex=0.7,
     cl.pos='n')

grid.echo()
P3_corr_refbased_reffree <- grid.grab()
```

``` r
matrix.colors <- getGrob(P3_corr_refbased_reffree, gPath("square"), grep = TRUE)[["gp"]][["fill"]]
P3_corr_refbased_reffree <- editGrob(P3_corr_refbased_reffree,
               gPath("square"), grep = TRUE,
               gp = gpar(col = NA,
                         fill = NA))

P3_corr_refbased_reffree<- editGrob(P3_corr_refbased_reffree,
               gPath("symbols-rect-1"), grep = TRUE,
               gp = gpar(fill = matrix.colors))

P3_corr_refbased_reffree <- editGrob(P3_corr_refbased_reffree,
               gPath("background"), grep = TRUE,
               gp = gpar(fill = NA))
```

``` r
r_cor_reffree_rpc_placenta_bet <- psych::corr.test(phenor_reffree_placenta_BET, phenor_rpc_placenta_BET, method="spearman", adjust = "bonferroni")
```

``` r
par(xpd=TRUE)
plot_r_cor_reffree_rpc_placenta_bet <- corrplot(r_cor_reffree_rpc_placenta_bet$r,
     method="color", number.digits=1,
     mar=c(1,6,0,6),
     addCoef.col = "black",
     p.mat = r_cor_reffree_rpc_placenta_bet$p,
     insig = "label_sig",
     sig.level = c(.001, .01),
     pch.cex = 0.75,
     pch.col = "black",
     tl.col="black",
     tl.cex=0.7,
     tl.srt = 45,
     number.cex=0.7,
     cl.pos='n')

grid.echo()
P4_corr_refbased_reffree <- grid.grab()
```

``` r
matrix.colors <- getGrob(P4_corr_refbased_reffree, gPath("square"), grep = TRUE)[["gp"]][["fill"]]
P4_corr_refbased_reffree <- editGrob(P4_corr_refbased_reffree,
               gPath("square"), grep = TRUE,
               gp = gpar(col = NA,
                         fill = NA))

P4_corr_refbased_reffree<- editGrob(P4_corr_refbased_reffree,
               gPath("symbols-rect-1"), grep = TRUE,
               gp = gpar(fill = matrix.colors))

P4_corr_refbased_reffree <- editGrob(P4_corr_refbased_reffree,
               gPath("background"), grep = TRUE,
               gp = gpar(fill = NA))
```

``` r
ggsave("Results/Cor_Refbased_Reffree_combined.pdf",
       ggarrange(P1_corr_refbased_reffree, P2_corr_refbased_reffree, P3_corr_refbased_reffree, P4_corr_refbased_reffree, ncol=1, labels = c("a", "b", "c", "d"), heights = c(2,4,0.8,1), widths=c(2,4,0.8,1)), width=84, height=160, units="mm", dpi=600, scale=2)
```

[to the top](#top)

# Correlation reference-free estimated cell types and phenotypes

*Fig. 2*

again need to adapt corrplot function by adding +0.25 / +0.35 to the
x-position

## CVS

``` r
CVS_ITU_filtered_corPhReffree <- na.omit(Pheno_CVS_ITU_filtered[ ,c("gestage_at_CVS_weeks", "PC1_ethnicity", "PC2_ethnicity", "Child_Sex", "C1", "C2", "C3", "C4", "C5")])
# n = 200
names(CVS_ITU_filtered_corPhReffree) <- c("gestational age\n(weeks)", "PC1 ethnicity", "PC2 ethnicity", "child sex", "C1  ", "C2  ", "C3  ", "C4  ", "C5  ")
cells_CVS_ITU_filtered_corPhReffree <- CVS_ITU_filtered_corPhReffree[ ,c("C1  ", "C2  ", "C3  ", "C4  ", "C5  ")]
phenos_CVS_ITU_filtered_corPhReffree <- CVS_ITU_filtered_corPhReffree[ ,c("gestational age\n(weeks)", "child sex", "PC1 ethnicity", "PC2 ethnicity")]
phenos_CVS_ITU_filtered_corPhReffree$`child sex` = as.numeric(phenos_CVS_ITU_filtered_corPhReffree$`child sex`)
```

``` r
r_CVS_ITU_filtered_corPhReffree <- psych::corr.test(phenos_CVS_ITU_filtered_corPhReffree, cells_CVS_ITU_filtered_corPhReffree, method="spearman", adjust = "bonferroni")
```

``` r
corrplot(r_CVS_ITU_filtered_corPhReffree$r,
     method="color", number.digits=1,
     mar=c(1,3,2,2),
     addCoef.col = "black",
     p.mat = r_CVS_ITU_filtered_corPhReffree$p,
     insig = "label_sig",
     sig.level = c(.001, .01),
     pch.cex = 0.7,
     pch.col = "black",
     tl.col="black",
     tl.cex=0.75,
     number.cex= 0.7,
     outline=TRUE,
     cl.pos='n')

grid.echo()
P1_corr_pheno_reffree <- grid.grab()
```

``` r
matrix.colors <- getGrob(P1_corr_pheno_reffree, gPath("square"), grep = TRUE)[["gp"]][["fill"]]
P1_corr_pheno_reffree <- editGrob(P1_corr_pheno_reffree,
               gPath("square"), grep = TRUE,
               gp = gpar(col = NA,
                         fill = NA))

# apply the saved colours to the underlying matrix grob
P1_corr_pheno_reffree <- editGrob(P1_corr_pheno_reffree,
               gPath("symbols-rect-1"), grep = TRUE,
               gp = gpar(fill = matrix.colors))

# convert the background fill from white to transparent, while we are at it
P1_corr_pheno_reffree <- editGrob(P1_corr_pheno_reffree,
               gPath("background"), grep = TRUE,
               gp = gpar(fill = NA))
```

## Placenta ITU

``` r
Placenta_ITU_reduced_filtered_corPhReffree <- na.omit(Pheno_Placenta_ITU_reduced_filtered[ ,c("Gestational_Age_Weeks", "PC1_ethnicity", "PC2_ethnicity", "Child_Sex", "caseVScontrol", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")])
# n = 425
names(Placenta_ITU_reduced_filtered_corPhReffree) <- c("    gestational age\n(weeks)", "PC1 ethnicity", "PC2 ethnicity", "child sex", "    fetal chromosomal\n testing (yes/no)", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
cells_Placenta_ITU_reduced_filtered_corPhReffree <- Placenta_ITU_reduced_filtered_corPhReffree[ ,c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")]
phenos_Placenta_ITU_reduced_filtered_corPhReffree <- Placenta_ITU_reduced_filtered_corPhReffree[ ,c("    gestational age\n(weeks)", "child sex", "PC1 ethnicity", "PC2 ethnicity", "    fetal chromosomal\n testing (yes/no)")]
phenos_Placenta_ITU_reduced_filtered_corPhReffree$`child sex` = as.numeric(Placenta_ITU_reduced_filtered_corPhReffree$`child sex`)
phenos_Placenta_ITU_reduced_filtered_corPhReffree$`    fetal chromosomal\n testing (yes/no)` = as.numeric(Placenta_ITU_reduced_filtered_corPhReffree$`    fetal chromosomal\n testing (yes/no)`)
```

``` r
r_Placenta_ITU_reduced_filtered_corPhReffree <- psych::corr.test(phenos_Placenta_ITU_reduced_filtered_corPhReffree, cells_Placenta_ITU_reduced_filtered_corPhReffree, method="spearman", adjust = "bonferroni")
```

``` r
corrplot(r_Placenta_ITU_reduced_filtered_corPhReffree$r,
     method="color", number.digits=1,
     mar=c(1,3,1,3),
     addCoef.col = "black",
     p.mat = r_Placenta_ITU_reduced_filtered_corPhReffree$p,
     insig = "label_sig",
     sig.level = c(.001, .01),
     pch.cex = 0.7,
     pch.col = "black",
     tl.col="black",
     tl.cex=0.75,
     number.cex= 0.7,
     outline=TRUE,
     cl.pos='n')

grid.echo()
P2_corr_pheno_reffree <- grid.grab()
```

``` r
matrix.colors <- getGrob(P2_corr_pheno_reffree, gPath("square"), grep = TRUE)[["gp"]][["fill"]]
P2_corr_pheno_reffree <- editGrob(P2_corr_pheno_reffree,
               gPath("square"), grep = TRUE,
               gp = gpar(col = NA,
                         fill = NA))

# apply the saved colours to the underlying matrix grob
P2_corr_pheno_reffree <- editGrob(P2_corr_pheno_reffree,
               gPath("symbols-rect-1"), grep = TRUE,
               gp = gpar(fill = matrix.colors))

# convert the background fill from white to transparent, while we are at it
P2_corr_pheno_reffree <- editGrob(P2_corr_pheno_reffree,
               gPath("background"), grep = TRUE,
               gp = gpar(fill = NA))
```

## Placenta PREDO

``` r
Placenta_PREDO_filtered_corPhReffree <- na.omit(Pheno_Placenta_PREDO_filtered[ ,c("Gestational_Age", "PC1_ethnicity", "PC2_ethnicity", "Child_Sex", "C1", "C2")])
# n = 118
names(Placenta_PREDO_filtered_corPhReffree) <- c("gestational age\n(weeks)", "PC1 ethnicity", "PC2 ethnicity", "child sex", "C1 ", "C2 ")
cells_Placenta_PREDO_filtered_corPhReffree <- Placenta_PREDO_filtered_corPhReffree[ ,c("C1 ", "C2 ")]
phenos_Placenta_PREDO_filtered_corPhReffree <- Placenta_PREDO_filtered_corPhReffree[ ,c("gestational age\n(weeks)", "child sex", "PC1 ethnicity", "PC2 ethnicity")]
phenos_Placenta_PREDO_filtered_corPhReffree$`child sex` = as.numeric(Placenta_PREDO_filtered_corPhReffree$`child sex`)
```

``` r
r_Placenta_PREDO_filtered_corPhReffree <- psych::corr.test(phenos_Placenta_PREDO_filtered_corPhReffree, cells_Placenta_PREDO_filtered_corPhReffree, method="spearman", adjust = "bonferroni")
```

``` r
corrplot(r_Placenta_PREDO_filtered_corPhReffree$r,
     method="color", number.digits=1,
     mar=c(1,0,1,1),
     addCoef.col = "black",
     p.mat = r_Placenta_PREDO_filtered_corPhReffree$p,
     insig = "label_sig",
     sig.level = c(.001, .01),
     pch.cex = 0.7,
     pch.col = "black",
     tl.col="black",
     tl.cex=0.75,
     number.cex=0.7,
     outline=TRUE,
     cl.pos='n')

grid.echo()
P3_corr_pheno_reffree <- grid.grab()
```

``` r
matrix.colors <- getGrob(P3_corr_pheno_reffree, gPath("square"), grep = TRUE)[["gp"]][["fill"]]
P3_corr_pheno_reffree <- editGrob(P3_corr_pheno_reffree,
               gPath("square"), grep = TRUE,
               gp = gpar(col = NA,
                         fill = NA))

# apply the saved colours to the underlying matrix grob
P3_corr_pheno_reffree <- editGrob(P3_corr_pheno_reffree,
               gPath("symbols-rect-1"), grep = TRUE,
               gp = gpar(fill = matrix.colors))

# convert the background fill from white to transparent, while we are at it
P3_corr_pheno_reffree <- editGrob(P3_corr_pheno_reffree,
               gPath("background"), grep = TRUE,
               gp = gpar(fill = NA))
```

## Placenta BET

``` r
Placenta_BET_filtered_corPhReffree <- na.omit(Pheno_Placenta_BET_filtered[ ,c("gestage_weeks", "PC1_ethnicity", "PC2_ethnicity","PC3_ethnicity", "PC4_ethnicity", "Sex", "group", "C1", "C2", "C3")])
# n = 136
names(Placenta_BET_filtered_corPhReffree) <- c("gestational age\n(weeks)", "PC1 ethnicity", "PC2 ethnicity", "PC3 ethnicity", "PC4 ethnicity", "child sex", "BET status", "C1", "C2", "C3")
cells_Placenta_BET_filtered_corPhReffree <- Placenta_BET_filtered_corPhReffree[ ,c("C1", "C2", "C3")]
phenos_Placenta_BET_filtered_corPhReffree <- Placenta_BET_filtered_corPhReffree[ ,c("gestational age\n(weeks)", "child sex", "BET status", "PC1 ethnicity", "PC2 ethnicity", "PC3 ethnicity", "PC4 ethnicity")]
phenos_Placenta_BET_filtered_corPhReffree$`child sex` = as.numeric(Placenta_BET_filtered_corPhReffree$`child sex`)
phenos_Placenta_BET_filtered_corPhReffree$`BET status` = as.numeric(Placenta_BET_filtered_corPhReffree$`BET status`)
```

``` r
r_Placenta_BET_filtered_corPhReffree <- psych::corr.test(phenos_Placenta_BET_filtered_corPhReffree, cells_Placenta_BET_filtered_corPhReffree, method="spearman", adjust = "bonferroni")
```

``` r
corrplot(r_Placenta_BET_filtered_corPhReffree$r,
     method="color", number.digits=1,
     mar=c(1,0,1,1),
     addCoef.col = "black",
     p.mat = r_Placenta_BET_filtered_corPhReffree$p,
     insig = "label_sig",
     sig.level = c(.001, .01),
     pch.cex = 0.7,
     pch.col = "black",
     tl.col="black",
     tl.cex=0.75,
     number.cex=0.7,
     outline=TRUE,
     cl.pos='n')

grid.echo()
P4_corr_pheno_reffree <- grid.grab()
```

``` r
matrix.colors <- getGrob(P4_corr_pheno_reffree, gPath("square"), grep = TRUE)[["gp"]][["fill"]]
P4_corr_pheno_reffree <- editGrob(P4_corr_pheno_reffree,
               gPath("square"), grep = TRUE,
               gp = gpar(col = NA,
                         fill = NA))

# apply the saved colours to the underlying matrix grob
P4_corr_pheno_reffree <- editGrob(P4_corr_pheno_reffree,
               gPath("symbols-rect-1"), grep = TRUE,
               gp = gpar(fill = matrix.colors))

# convert the background fill from white to transparent, while we are at it
P4_corr_pheno_reffree <- editGrob(P4_corr_pheno_reffree,
               gPath("background"), grep = TRUE,
               gp = gpar(fill = NA))
```

``` r
ggsave("Results/Cor_Pheno_Reffree_combined.pdf",
       ggarrange(P1_corr_pheno_reffree, P2_corr_pheno_reffree, P3_corr_pheno_reffree, P4_corr_pheno_reffree, ncol=1, labels = c("a", "b", "c", "d"), heights = c(1.5,2,1.5,2), widths = c(2,2,1,1), align='v'), width=84, height=180, units="mm", dpi=600, scale=2)
```

[to the top](#top)

# Check if methylation can be predicted using cell types (using CV)

## CVS

``` r
Reg_Pheno_CVS_ITU_filtered <- na.omit(Pheno_CVS_ITU_filtered[ ,c("gestage_at_CVS_weeks", "PC1_ethnicity", "PC2_ethnicity", "Child_Sex", "Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "Syncytiotrophoblast", "nRBC","PC_cell1", "PC_cell2", "PC_cell3","PC_cell4", "C1", "C2", "C3", "C4", "C5", "PC_methB1", "Sample_Name")])
# n = 200
dim(Reg_Pheno_CVS_ITU_filtered)

save(Reg_Pheno_CVS_ITU_filtered, file="Input_Data_prepared/Reg_Pheno_CVS_ITU_filtered.Rdata")
```

### predict PC1 methylation

plot cell & PC 1

``` r
cell1_reg <- Reg_Pheno_CVS_ITU_filtered[,c(5:10,20)]
cell1_reg %>%
  gather(-PC_methB1, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y =PC_methB1)) +
  facet_wrap(~ var, scales = "free") +
  geom_point() +
  theme_bw()+
  stat_smooth()

cell2_reg <- Reg_Pheno_CVS_ITU_filtered[,c(15:19,20)]
cell2_reg %>%
  gather(-PC_methB1, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y =PC_methB1)) +
  facet_wrap(~ var, scales = "free") +
  geom_point() +
  theme_bw()+
  stat_smooth()
```

-RMSE-

``` r
models = vector(mode = "list", length = 6) 
models[[1]] = PC_methB1 ~ 1
models[[2]] = PC_methB1 ~ gestage_at_CVS_weeks + Child_Sex + PC1_ethnicity + PC2_ethnicity

models[[3]] = PC_methB1 ~  Trophoblasts +Stromal + Hofbauer + Endothelial + Syncytiotrophoblast + nRBC
models[[4]] = PC_methB1 ~  Trophoblasts +Stromal + Hofbauer + Endothelial + Syncytiotrophoblast + nRBC+ gestage_at_CVS_weeks + Child_Sex + PC1_ethnicity + PC2_ethnicity

models[[5]] = PC_methB1 ~ C1 + C2 + C3 + C4 + C5
models[[6]] = PC_methB1 ~ C1 + C2 + C3 + C4 + C5+ gestage_at_CVS_weeks + Child_Sex+ PC1_ethnicity + PC2_ethnicity

CV_output_CVS_Meth = xval.glm(data = Reg_Pheno_CVS_ITU_filtered, models, folds = 10, repeats = 500, seed=200)
```

``` r
save(CV_output_CVS_Meth, file="Results/RData/CV_output_CVS_Meth.Rdata")
```

get adjusted R2 from winning model

``` r
rsq(glm(models[[3]], data= Reg_Pheno_CVS_ITU_filtered),adj=TRUE,type='v')
```

``` r
summary(glm(models[[3]], data= Reg_Pheno_CVS_ITU_filtered))
```

``` r
# difference null and residual variance
12220566-1177114
199-194
#p-value = 1 - pchisq(deviance, degrees of freedom)
1-pchisq(11043452,5)
# -> significant
```

*Plots*

``` r
load("Results/RData/CV_output_CVS_Meth.Rdata")
```

``` r
# this is how plots were defined in original package function
p2c <- CV_output_CVS_Meth[["den.plot"]]
p1c <- CV_output_CVS_Meth[["stab.plot"]]
pc <- CV_output_CVS_Meth[["box.plot"]]
```

``` r
p1c <- p1c +
  theme_minimal()+
  theme(axis.text.x = element_text(size=8), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5), legend.text=element_text(size=8), legend.title=element_text(size=8), legend.key.size = unit(0.2, 'cm'))+
  labs(tag = "a ", title="CVS (ITU)")+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25))+
  scale_color_discrete(name="Models",labels=c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6"))
```

``` r
pc <- pc +
  theme_minimal()+
  theme(axis.text.x = element_text(size=10), axis.title.x=element_text(size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), legend.position="none", axis.text.x.top = element_text(size=10))+
  labs(tag = "    ")
```

``` r
p2c <- p2c +
  theme_minimal()+
  theme(legend.position="none", axis.text.x = element_text(size=10))
```

``` r
# define parameters we need for plot
my.ylab <- "RMSE" 
K <- 10
repeats <- 500
wins <- CV_output_CVS_Meth$wins
```

``` r
# titleplot
titleplot <- ggplot() + geom_point(aes(1,1), colour="white") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank()) + annotate("text", x=1, y=1,
                                                                               label=paste0(my.ylab,'\n (',K,'-fold, ',repeats,' repeats) \nModel ',which.max(wins),' wins.'))
```

``` r
ggCV_CVS_M <- grid.arrange(p1c, titleplot, pc, p2c, ncol=2, nrow=2, widths=c(5, 2), heights=c(2.5, 5))
```

``` r
ggsave("Results/CV_CVS_M.pdf",
       grid.arrange(p1c, titleplot, pc, p2c, ncol=2, nrow=2, widths=c(5, 2), heights=c(2, 5)),
       width=84, height=50, units="mm", dpi=600, scale=2, device = cairo_pdf)
```

[to the top](#top)

## Placenta ITU

``` r
Reg_Pheno_Placenta_ITU_reduced_filtered <- na.omit(Pheno_Placenta_ITU_reduced_filtered[ ,c("Gestational_Age_Weeks", "PC1_ethnicity", "PC2_ethnicity", "Child_Sex", "Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "Syncytiotrophoblast", "nRBC", "PC_cell1", "PC_cell2", "PC_cell3", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "PC_methB1", "Sample_Name", "caseVScontrol")])

dim(Reg_Pheno_Placenta_ITU_reduced_filtered)
save(Reg_Pheno_Placenta_ITU_reduced_filtered, file="Input_Data_prepared/Reg_Pheno_Placenta_ITU_reduced_filtered.Rdata")
```

### predict PC1 methylation

plot PC and cell types

``` r
cell1_reg <- Reg_Pheno_Placenta_ITU_reduced_filtered[,c(5:10,22)]
cell1_reg %>%
  gather(-PC_methB1, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y =PC_methB1)) +
  facet_wrap(~ var, scales = "free") +
  geom_point() +
  theme_bw()+
  stat_smooth()

cell2_reg <- Reg_Pheno_Placenta_ITU_reduced_filtered[,c(14:21, 22)]
cell2_reg %>%
  gather(-PC_methB1, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y =PC_methB1)) +
  facet_wrap(~ var, scales = "free") +
  geom_point() +
  theme_bw()+
  stat_smooth()
```

-RMSE-

``` r
models = vector(mode = "list", length = 6) 
models[[1]] = PC_methB1 ~ 1
models[[2]] = PC_methB1 ~ Gestational_Age_Weeks + Child_Sex + PC1_ethnicity + PC2_ethnicity

models[[3]] = PC_methB1 ~ Trophoblasts + Stromal + Hofbauer + Endothelial + Syncytiotrophoblast + nRBC
models[[4]] = PC_methB1 ~  Trophoblasts + Stromal + Hofbauer + Endothelial + Syncytiotrophoblast + nRBC + Gestational_Age_Weeks + Child_Sex + PC1_ethnicity + PC2_ethnicity

models[[5]] = PC_methB1 ~ C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8
models[[6]] = PC_methB1 ~ C1 + C2 + C3 + C4 + C5+ C6 + C7 + C8+ Gestational_Age_Weeks + Child_Sex + PC1_ethnicity + PC2_ethnicity

CV_output_Placenta_Meth = xval.glm(data = Reg_Pheno_Placenta_ITU_reduced_filtered, models, folds = 10, repeats = 500, seed=200)
```

``` r
save(CV_output_Placenta_Meth, file="Results/RData/CV_output_Placenta_Meth.Rdata")
```

get adjusted R2 from winning model

``` r
rsq(glm(models[[6]], data= Reg_Pheno_Placenta_ITU_reduced_filtered),adj=TRUE,type='v')
```

``` r
summary(glm(models[[6]], data= Reg_Pheno_Placenta_ITU_reduced_filtered))
```

``` r
# difference null and residual variance
28342062-2080515
424-411
#p-value = 1 - pchisq(deviance, degrees of freedom)
1-pchisq(26261547, 13)
# -> significant
```

-RMSE- see if inclusion of case/control changes something

``` r
models = vector(mode = "list", length = 6) 
models[[1]] = PC_methB1 ~ 1
models[[2]] = PC_methB1 ~ Gestational_Age_Weeks + Child_Sex + PC1_ethnicity + PC2_ethnicity + caseVScontrol

models[[3]] = PC_methB1 ~ Trophoblasts + Stromal + Hofbauer + Endothelial + Syncytiotrophoblast + nRBC
models[[4]] = PC_methB1 ~  Trophoblasts + Stromal + Hofbauer + Endothelial + Syncytiotrophoblast + nRBC + Gestational_Age_Weeks + Child_Sex + PC1_ethnicity + PC2_ethnicity + caseVScontrol

models[[5]] = PC_methB1 ~ C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8
models[[6]] = PC_methB1 ~ C1 + C2 + C3 + C4 + C5+ C6 + C7 + C8+ Gestational_Age_Weeks + Child_Sex + PC1_ethnicity + PC2_ethnicity + caseVScontrol

CV_output_Placenta_Meth_casecontrol = xval.glm(data = Reg_Pheno_Placenta_ITU_reduced_filtered, models, folds = 10, repeats = 500, seed=200)
```

``` r
save(CV_output_Placenta_Meth_casecontrol, file="Results/RData/CV_output_Placenta_Meth_casecontrol.Rdata")
```

*Plots*

``` r
load("Results/RData/CV_output_Placenta_Meth.Rdata")
```

``` r
# this is how plots were defined in original package function
pp2 <- CV_output_Placenta_Meth[["den.plot"]]
pp1 <- CV_output_Placenta_Meth[["stab.plot"]]
pp <- CV_output_Placenta_Meth[["box.plot"]]
```

``` r
pp1 <- pp1 +
  theme_minimal()+
  theme(axis.text.x = element_text(size=8), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5), legend.text=element_text(size=8), legend.title=element_text(size=8), legend.key.size = unit(0.2, 'cm'))+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25))+
  labs(tag = "b", title="Placenta (ITU)")+
  scale_color_discrete(name="Models",labels=c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6"))
```

``` r
pp <- pp +
  theme_minimal()+
  theme(axis.text.x = element_text(size=10), axis.title.x=element_text(size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), legend.position="none", axis.text.x.top = element_text(size=10))+
  labs(tag= "    ")
```

``` r
pp2 <- pp2 +
  theme_minimal()+
  theme(legend.position="none", axis.text.x = element_text(size=10))
```

``` r
# define parameters we need for plot
my.ylab <- "RMSE" 
K <- 10
repeats <- 500
wins <- CV_output_Placenta_Meth$wins
```

``` r
# titleplot
titleplot <- ggplot() + geom_point(aes(1,1), colour="white") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank()) + annotate("text", x=1, y=1,
                                                                               label=paste0(my.ylab,'\n (',K,'-fold, ',repeats,' repeats) \nModel ',which.max(wins),' wins.'))
```

``` r
ggCV_PlacentaITU_M <- grid.arrange(pp1, titleplot, pp, pp2, ncol=2, nrow=2, widths=c(5, 2), heights=c(2.5, 5))
```

``` r
ggsave("Results/CV_PlacentaITU_M.pdf",
       ggCV_PlacentaITU_M,
       width=100, height=80, units="mm", dpi=600, scale=2)
```

## Placenta PREDO

### predict PC1 methylation

``` r
Reg_Pheno_Placenta_PREDO <- na.omit(Pheno_Placenta_PREDO_filtered[ ,c("Gestational_Age", "PC1_ethnicity", "PC2_ethnicity", "Child_Sex", "Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "Syncytiotrophoblast", "nRBC", "PC_cell1", "PC_cell2", "PC_cell3", "C1", "C2", "PC_methB1", "ID")])
dim(Reg_Pheno_Placenta_PREDO)

save(Reg_Pheno_Placenta_PREDO, file="Input_Data_prepared/Reg_Pheno_Placenta_PREDO.Rdata")
```

1st PC methylation and cell types

``` r
cell1_reg <- Reg_Pheno_Placenta_PREDO[,c(5:10,16)]
cell1_reg %>%
  gather(-PC_methB1, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y =PC_methB1)) +
  facet_wrap(~ var, scales = "free") +
  geom_point() +
  theme_bw()+
  stat_smooth()

cell2_reg <- Reg_Pheno_Placenta_PREDO[,c(14:15,16)]
cell2_reg %>%
  gather(-PC_methB1, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y =PC_methB1)) +
  facet_wrap(~ var, scales = "free") +
  geom_point() +
  theme_bw()+
  stat_smooth()
```

-RMSE-

``` r
# Cross-validation
models = vector(mode = "list", length = 6) 
models[[1]] = PC_methB1 ~ 1
models[[2]] = PC_methB1 ~ Gestational_Age + Child_Sex +PC1_ethnicity + PC2_ethnicity

models[[3]] = PC_methB1 ~ Trophoblasts + Stromal + Hofbauer + Endothelial + Syncytiotrophoblast+ nRBC
models[[4]] = PC_methB1 ~ Trophoblasts + Stromal + Hofbauer + Endothelial + Syncytiotrophoblast+ nRBC+ Gestational_Age + Child_Sex + PC1_ethnicity + PC2_ethnicity 

models[[5]] = PC_methB1 ~ C1 + C2 
models[[6]] = PC_methB1 ~ C1 + C2 + Gestational_Age + Child_Sex + PC1_ethnicity + PC2_ethnicity 

CV_output_Placenta_predo_meth = xval.glm(data = Reg_Pheno_Placenta_PREDO, models, folds = 10, repeats = 500, seed = 200)
```

get adjusted R2 from winning model

``` r
rsq(glm(models[[4]], data= Reg_Pheno_Placenta_PREDO),adj=TRUE,type='v')
```

``` r
summary(glm(models[[4]], data= Reg_Pheno_Placenta_PREDO))
```

``` r
# difference null and residual variance
9079814-1190407
117-109
#p-value = 1 - pchisq(deviance, degrees of freedom)
1-pchisq(7889407, 8)
# -> significant
```

``` r
save(CV_output_Placenta_predo_meth, file="Results/RData/CV_output_Placenta_predo_meth.Rdata")
```

*Plots*

``` r
load("Results/RData/CV_output_Placenta_predo_meth.Rdata")
```

``` r
# this is how plots were defined in original package function
p2 <- CV_output_Placenta_predo_meth[["den.plot"]]
p1 <- CV_output_Placenta_predo_meth[["stab.plot"]]
p <- CV_output_Placenta_predo_meth[["box.plot"]]
```

``` r
p1 <- p1 +
  theme_minimal()+
  theme(axis.text.x = element_text(size=8), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5), legend.text=element_text(size=8), legend.title=element_text(size=8), legend.key.size = unit(0.2, 'cm'))+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25))+
  labs(tag = "c ", title="Placenta (PREDO)")+
    scale_color_discrete(name="Models",labels=c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6"))
```

``` r
p <- p +
  theme_minimal()+
  theme(axis.text.x = element_text(size=10), axis.title.x=element_text(size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), legend.position="none", axis.text.x.top = element_text(size=10))+
  labs(tag= "    ")
```

``` r
p2 <- p2 +
  theme_minimal()+
  theme(legend.position="none", axis.text.x = element_text(size=10))
```

``` r
# define parameters we need for plot
my.ylab <- "RMSE" 
K <- 10
repeats <- 500
wins <- CV_output_Placenta_predo_meth$wins
```

``` r
# titleplot
titleplot <- ggplot() + geom_point(aes(1,1), colour="white") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank()) + annotate("text", x=1, y=1,
                                                                               label=paste0(my.ylab,'\n (',K,'-fold, ',repeats,' repeats) \nModel ',which.max(wins),' wins.'))
```

``` r
ggCV_PlacentaPREDO_M <- grid.arrange(p1, titleplot, p, p2, ncol=2, nrow=2, widths=c(5, 2), heights=c(2.5, 5))
```

``` r
ggsave("Results/CV_PlacentaPREDO_M.pdf",
       ggCV_PlacentaPREDO_M,
       width=84, height=50, units="mm", dpi=600, scale=2)
```

## Placenta BET

### predict PC1 methylation

``` r
Reg_Pheno_Placenta_BET <- na.omit(Pheno_Placenta_BET_filtered[ ,c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "Syncytiotrophoblast", "nRBC", "C1", "C2","C3", "PC_methB1", "PC_methB2", "PC_methB3", "Sample_Name", "Sex", "PC1_ethnicity", "PC2_ethnicity", "PC3_ethnicity", "PC4_ethnicity", "gestage_weeks", "group")])

dim(Reg_Pheno_Placenta_BET)

save(Reg_Pheno_Placenta_BET, file="Input_Data_prepared/Reg_Pheno_Placenta_BET.Rdata")
```

1st PC methylation and cell types

``` r
cell1_reg <- Reg_Pheno_Placenta_BET[,c(1:6,10)]
cell1_reg %>%
  gather(-PC_methB1, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y =PC_methB1)) +
  facet_wrap(~ var, scales = "free") +
  geom_point() +
  stat_smooth()+
  theme_bw()

cell2_reg <- Reg_Pheno_Placenta_BET[,c(7:9,10)]
cell2_reg %>%
  gather(-PC_methB1, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y =PC_methB1)) +
  facet_wrap(~ var, scales = "free") +
  geom_point() +
  stat_smooth()+
  theme_bw()
```

-RMSE, model with all samples-

``` r
# Cross-validation
models = vector(mode = "list", length = 6) 
models[[1]] = PC_methB1 ~ 1
models[[2]] = PC_methB1 ~ gestage_weeks + Sex +PC1_ethnicity + PC2_ethnicity + PC3_ethnicity + PC4_ethnicity 

models[[3]] = PC_methB1 ~ Trophoblasts + Stromal + Hofbauer + Endothelial + Syncytiotrophoblast+ nRBC
models[[4]] = PC_methB1 ~ Trophoblasts + Stromal + Hofbauer + Endothelial + Syncytiotrophoblast+ nRBC+ gestage_weeks + Sex +PC1_ethnicity + PC2_ethnicity + PC3_ethnicity + PC4_ethnicity 

models[[5]] = PC_methB1 ~ C1 + C2 + C3
models[[6]] = PC_methB1 ~ C1 + C2+ C3 + gestage_weeks + Sex +PC1_ethnicity + PC2_ethnicity + PC3_ethnicity + PC4_ethnicity

CV_output_Placenta_BET_meth = xval.glm(data = Reg_Pheno_Placenta_BET, models, folds = 10, repeats = 500, seed = 200)
```

``` r
save(CV_output_Placenta_BET_meth, file="Results/RData/CV_output_Placenta_BET_meth.Rdata")
```

*identify conspicious samples* extract ref-based cell types

``` r
refbased_cells_Pheno_BET_filtered <- Pheno_Placenta_BET_filtered[ ,c(names_refbased_cells, "Sample_Name")]
```

IQR \> 3

``` r
IQR_outliers3_refbased_cells_Pheno_BET_filtered <- cbind(data.frame(sapply(refbased_cells_Pheno_BET_filtered[,!names(refbased_cells_Pheno_BET_filtered) %in% "Sample_Name"], function(x) is_outlier(x, coef=3))), refbased_cells_Pheno_BET_filtered$Sample_Name)

names(IQR_outliers3_refbased_cells_Pheno_BET_filtered)[names(IQR_outliers3_refbased_cells_Pheno_BET_filtered) == "refbased_cells_Pheno_BET_filtered$Sample_Name"] <- "Sample_Name"
```

``` r
IQR_outliers3_names <- unique(c(IQR_outliers3_refbased_cells_Pheno_BET_filtered[IQR_outliers3_refbased_cells_Pheno_BET_filtered$nRBC == TRUE, "Sample_Name"], IQR_outliers3_refbased_cells_Pheno_BET_filtered[IQR_outliers3_refbased_cells_Pheno_BET_filtered$Hofbauer == TRUE, "Sample_Name"]))
```

exclude the 5 outliers in Hofbauer & nRBC

``` r
Reg_Pheno_Placenta_BET_cleaned <- Reg_Pheno_Placenta_BET[! Reg_Pheno_Placenta_BET$Sample_Name %in% IQR_outliers3_names, ]
```

-RMSE- with ‘cleaned’ model

``` r
# Cross-validation
models = vector(mode = "list", length = 6) 
models[[1]] = PC_methB1 ~ 1
models[[2]] = PC_methB1 ~ gestage_weeks + Sex +PC1_ethnicity + PC2_ethnicity + PC3_ethnicity + PC4_ethnicity 

models[[3]] = PC_methB1 ~ Trophoblasts + Stromal + Hofbauer + Endothelial + Syncytiotrophoblast+ nRBC
models[[4]] = PC_methB1 ~ Trophoblasts + Stromal + Hofbauer + Endothelial + Syncytiotrophoblast+ nRBC+ gestage_weeks + Sex +PC1_ethnicity + PC2_ethnicity + PC3_ethnicity + PC4_ethnicity 

models[[5]] = PC_methB1 ~ C1 + C2 + C3
models[[6]] = PC_methB1 ~ C1 + C2+ C3 + gestage_weeks + Sex +PC1_ethnicity + PC2_ethnicity + PC3_ethnicity + PC4_ethnicity

CV_output_Placenta_BET_meth_cleaned = xval.glm(data = Reg_Pheno_Placenta_BET_cleaned, models, folds = 10, repeats = 500, seed = 200)
```

``` r
save(CV_output_Placenta_BET_meth_cleaned, file="Results/RData/CV_output_Placenta_BET_meth_cleaned.Rdata")
```

-RMSE- with ‘cleaned model’ and BET status included

``` r
# Cross-validation
models = vector(mode = "list", length = 6) 
models[[1]] = PC_methB1 ~ 1
models[[2]] = PC_methB1 ~ gestage_weeks + Sex +PC1_ethnicity + PC2_ethnicity + PC3_ethnicity + PC4_ethnicity + group

models[[3]] = PC_methB1 ~ Trophoblasts + Stromal + Hofbauer + Endothelial + Syncytiotrophoblast+ nRBC
models[[4]] = PC_methB1 ~ Trophoblasts + Stromal + Hofbauer + Endothelial + Syncytiotrophoblast+ nRBC+ gestage_weeks + Sex +PC1_ethnicity + PC2_ethnicity + PC3_ethnicity + PC4_ethnicity + group

models[[5]] = PC_methB1 ~ C1 + C2 + C3
models[[6]] = PC_methB1 ~ C1 + C2+ C3 + gestage_weeks + Sex +PC1_ethnicity + PC2_ethnicity + PC3_ethnicity + PC4_ethnicity+ group

CV_output_Placenta_BET_meth_group_cleaned = xval.glm(data = Reg_Pheno_Placenta_BET_cleaned, models, folds = 10, repeats = 500, seed = 200)
```

get adjusted R2 from winning model

``` r
rsq(glm(models[[3]], data= Reg_Pheno_Placenta_BET_cleaned),adj=TRUE,type='v')
```

``` r
summary(glm(models[[3]], data= Reg_Pheno_Placenta_BET_cleaned))
```

``` r
# difference null and residual variance
6621201-932287
130-126
#p-value = 1 - pchisq(deviance, degrees of freedom)
1-pchisq(5688914, 4)
# -> significant
```

*Plots*

``` r
load("Results/RData/CV_output_Placenta_BET_meth_cleaned.Rdata")
```

``` r
# this is how plots were defined in original package function
pe2 <- CV_output_Placenta_BET_meth_cleaned[["den.plot"]]
pe1 <- CV_output_Placenta_BET_meth_cleaned[["stab.plot"]]
pe <- CV_output_Placenta_BET_meth_cleaned[["box.plot"]]
```

``` r
pe1 <- pe1 +
  theme_minimal()+
  theme(axis.text.x = element_text(size=8), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5),legend.text=element_text(size=8), legend.title=element_text(size=8), legend.key.size = unit(0.2, 'cm'))+
  labs(tag = "d ", title="Placenta (BET)")+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25))+
  scale_color_discrete(name="Models",labels=c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6"))
```

``` r
pe <- pe +
  theme_minimal()+
  theme(axis.text.x = element_text(size=10), axis.title.x=element_text(size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), legend.position="none", axis.text.x.top = element_text(size=10))+
  labs(tag= "    ")
```

``` r
pe2 <- pe2 +
  theme_minimal()+
  theme(legend.position="none", axis.text.x = element_text(size=10))
```

``` r
# define parameters we need for plot
my.ylab <- "RMSE" 
K <- 10
repeats <- 500
wins <- CV_output_Placenta_BET_meth_cleaned$wins
```

``` r
# titleplot
titleplot <- ggplot() + geom_point(aes(1,1), colour="white") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank()) + annotate("text", x=1, y=1,
                                                                               label=paste0(my.ylab,'\n (',K,'-fold, ',repeats,' repeats) \nModel ',which.max(wins),' wins.'))
```

``` r
ggCV_PlacentaBET_M <- grid.arrange(pe1, titleplot, pe, pe2, ncol=2, nrow=2, widths=c(5, 2), heights=c(2.5, 5))
```

``` r
ggsave("Results/CV_PlacentaBET_M.pdf",
       ggCV_PlacentaBET_M,
       width=84, height=50, units="mm", dpi=600, scale=2)
```

### Make outlier plot for BET CV Model

this is an adaption to the plotting with facets I need only for this
plot:

``` r
Facet <- ggproto(
  init_scales = function(layout, x_scale = NULL, y_scale = NULL, params) {
    scales <- list()
    if (!is.null(x_scale)) {
      scales$x <- plyr::rlply(max(layout$SCALE_X), x_scale$clone())
    }
    if (!is.null(y_scale)) {
      scales$y <- plyr::rlply(max(layout$SCALE_Y), y_scale$clone())
    }
    scales
  },
)

scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
    shrink = facet_super$shrink,
    params = facet_super$params
  )
}
```

``` r
plot_scat_refbased_bet <- Reg_Pheno_Placenta_BET[,c("PC_methB1", names_refbased_cells)] %>%
  gather(-PC_methB1, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y =PC_methB1)) +
  facet_wrap(~ var, scales = "free_x") +
  geom_point() +
  theme_bw()+
  labs(x="", y="PC1 methylation")+
    theme(axis.text.x = element_text(size=6, angle=30, hjust=1), axis.title.x=element_text(size=6), axis.title.y = element_text(size=6), axis.text.y = element_text(size=6), plot.margin=unit(c(8,8,8,8),"points"), strip.text = element_text(size = 6))
```

``` r
pe_o <- CV_output_Placenta_BET_meth[["box.plot"]]
```

``` r
pe_o <- pe_o +
  theme_minimal()+
  theme(axis.text.x = element_text(size=8), axis.title.x=element_text(size=8), axis.title.y = element_text(size=8), axis.text.y = element_text(size=8), legend.position="none", axis.text.x.top = element_text(size=8))
```

``` r
plot_box_refbased_bet <- boxplot_fun(x = names_refbased_cells, id = "Sample_Name", dat= Reg_Pheno_Placenta_BET)
```

``` r
plot_scat_refbased_bet_custom <- plot_scat_refbased_bet +
  facet_wrap_custom(~var, scales = "free", ncol = 3, scale_overrides = list(
    scale_override(1, scale_x_continuous(breaks = c(0,0.05,0.10,0.15,0.20), limits=c(0,0.20), labels = label_number(accuracy = .01))),
    scale_override(2, scale_x_continuous(breaks = c(0,0.01,0.02,0.03,0.04), limits=c(0,0.045), labels = label_number(accuracy = .01))),
    scale_override(3, scale_x_continuous(breaks = c(0,0.01,0.02,0.03,0.04), limits=c(0,0.045), labels = label_number(accuracy = .01))),
    scale_override(4, scale_x_continuous(breaks = c(0,0.05,0.10,0.15,0.20), limits=c(0,0.21), labels = label_number(accuracy = .01))),
    scale_override(5, scale_x_continuous(breaks = c(0,0.20,0.40,0.60,0.80), limits=c(0,0.85), labels = label_number(accuracy = .01))),
    scale_override(6, scale_x_continuous(breaks = c(0,0.10,0.20,0.30,0.40), limits=c(0,0.4), labels = label_number(accuracy = .01)))
  ))
```

``` r
ggsave("Results/BET_CV_Model_Outlier_Plots.pdf",
       ggarrange(pe_o, plot_scat_refbased_bet_custom, plot_box_refbased_bet, ncol=2, nrow=2, labels = c("a", "b", "c"), heights = c(1,1,1), widths=c(1,1.2,1)), width=84, height=100, units="mm", dpi=600, scale=2)
```

## arrange plots together

*Fig. 3* Cross-Validation Models

``` r
ggA_M <- grid.arrange(arrangeGrob(ggCV_CVS_M, ggCV_PlacentaITU_M, ggCV_PlacentaPREDO_M, ggCV_PlacentaBET_M, nrow=4))
```

``` r
ggsave("Results/Models_Meth_combined.pdf",
       grid.arrange(arrangeGrob(ggCV_CVS_M, ggCV_PlacentaITU_M, ggCV_PlacentaPREDO_M, ggCV_PlacentaBET_M, nrow=4)), width=84, height=180, units="mm", dpi=600, scale=2)
```

[to the top](#top)

# Check how methylation of single CpGs can be predicted using cell types

## CVS

``` r
results_cvs_itu_lm_cpg_rb_filtered <-matrix(ncol=1,nrow=dim(Beta_CVS_ITU_filtered)[1])
results_cvs_itu_lm_cpg_rb_filtered[,1]<-rownames(Beta_CVS_ITU_filtered)

results_cvs_itu_lm_cpg_rf_filtered <-matrix(ncol=1,nrow=dim(Beta_CVS_ITU_filtered)[1])
results_cvs_itu_lm_cpg_rf_filtered[,1]<-rownames(Beta_CVS_ITU_filtered)
```

``` r
# ref-based
for (i in 1:dim(Beta_CVS_ITU_filtered)[1]) {
results_cvs_itu_lm_cpg_rb_filtered[i, 1] <- summary(lm(Beta_CVS_ITU_filtered[i,]~ Pheno_CVS_ITU_filtered$Trophoblasts + Pheno_CVS_ITU_filtered$Syncytiotrophoblast + Pheno_CVS_ITU_filtered$nRBC + Pheno_CVS_ITU_filtered$Endothelial + Pheno_CVS_ITU_filtered$Stromal + Pheno_CVS_ITU_filtered$Hofbauer))$adj.r.squared
}
```

``` r
# ref-free
for (i in 1:dim(Beta_CVS_ITU_filtered)[1]) {
  results_cvs_itu_lm_cpg_rf_filtered[i, 1] <- summary(lm(Beta_CVS_ITU_filtered[i,]~ Pheno_CVS_ITU_filtered$C1 + Pheno_CVS_ITU_filtered$C2 + Pheno_CVS_ITU_filtered$C3 + Pheno_CVS_ITU_filtered$C4 + Pheno_CVS_ITU_filtered$C5))$adj.r.squared
}
```

``` r
colnames(results_cvs_itu_lm_cpg_rb_filtered) <- c("adj.Rsquared")
results_cvs_itu_lm_cpg_rb_filtered <- data.frame(results_cvs_itu_lm_cpg_rb_filtered)
rownames(results_cvs_itu_lm_cpg_rb_filtered) <- rownames(Beta_CVS_ITU_filtered)

colnames(results_cvs_itu_lm_cpg_rf_filtered) <- c("adj.Rsquared")
results_cvs_itu_lm_cpg_rf_filtered <- data.frame(results_cvs_itu_lm_cpg_rf_filtered)
rownames(results_cvs_itu_lm_cpg_rf_filtered) <- rownames(Beta_CVS_ITU_filtered)
```

``` r
results_cvs_itu_lm_cpg_rb_filtered$adj.Rsquared[results_cvs_itu_lm_cpg_rb_filtered$adj.Rsquared <0] <- 0
results_cvs_itu_lm_cpg_rb_filtered$adj.Rsquared <- as.numeric(results_cvs_itu_lm_cpg_rb_filtered$adj.Rsquared)

results_cvs_itu_lm_cpg_rf_filtered$adj.Rsquared[results_cvs_itu_lm_cpg_rf_filtered$adj.Rsquared <0] <- 0
results_cvs_itu_lm_cpg_rf_filtered$adj.Rsquared <- as.numeric(results_cvs_itu_lm_cpg_rf_filtered$adj.Rsquared)
```

``` r
save(results_cvs_itu_lm_cpg_rb_filtered, file="Results/RData/results_cvs_itu_lm_cpg_rb_filtered.Rdata")
save(results_cvs_itu_lm_cpg_rf_filtered, file="Results/RData/results_cvs_itu_lm_cpg_rf_filtered.Rdata")
```

*check mean and sd in r2*

``` r
load("Results/RData/results_cvs_itu_lm_cpg_rb_filtered.Rdata")
load("Results/RData/results_cvs_itu_lm_cpg_rf_filtered.Rdata")
```

``` r
# mean and sd reference-based
mean(results_cvs_itu_lm_cpg_rb_filtered$adj.Rsquared)
sd(results_cvs_itu_lm_cpg_rb_filtered$adj.Rsquared)

cpgs_hist_cvs_rb_filtered <- ggplot(data=results_cvs_itu_lm_cpg_rb_filtered, aes(adj.Rsquared)) + 
  geom_histogram()+
  geom_vline(xintercept = 0.3, color="red", linetype="dotted")+
  xlab(as.expression(bquote("adjusted" ~ R^2)))+ ylab("count (CpGs)")+
  theme_bw()+ scale_x_continuous(breaks = seq(0, 0.8, 0.2))+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title = element_text(size=10), axis.title.y=element_text(size=10))+
  labs(title="reference-based")

cpgs_hist_cvs_rb_filtered
```

``` r
# mean and sd reference-free
mean(results_cvs_itu_lm_cpg_rf_filtered$adj.Rsquared)
sd(results_cvs_itu_lm_cpg_rf_filtered$adj.Rsquared)

cpgs_hist_cvs_rf_filtered <- ggplot(data=results_cvs_itu_lm_cpg_rf_filtered, aes(adj.Rsquared)) + 
  geom_histogram()+
  geom_vline(xintercept = 0.3, color="red", linetype="dotted")+
  xlab(as.expression(bquote("adjusted" ~ R^2)))+ ylab("count (CpGs)")+
  theme_bw()+ scale_x_continuous(breaks = seq(0, 0.8, 0.2))+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title = element_text(size=10), axis.title.y=element_text(size=10))+
  labs(title="reference-free")

cpgs_hist_cvs_rf_filtered
```

``` r
M_CpGs_CVS <- ggarrange(cpgs_hist_cvs_rb_filtered, cpgs_hist_cvs_rf_filtered, 
                        nrow = 1,
                        labels = c("a"),
                        align = "hv")

M_CpGs_CVS
```

[to the top](#top)

## Placenta ITU

``` r
results_placenta_itu_lm_cpg_rb_filtered <-matrix(ncol=1,nrow=dim(Beta_Placenta_ITU_reduced_filtered)[1])
results_placenta_itu_lm_cpg_rb_filtered[,1]<-rownames(Beta_Placenta_ITU_reduced_filtered)

results_placenta_itu_lm_cpg_rf_filtered <-matrix(ncol=1,nrow=dim(Beta_Placenta_ITU_reduced_filtered)[1])
results_placenta_itu_lm_cpg_rf_filtered[,1]<-rownames(Beta_Placenta_ITU_reduced_filtered)
```

``` r
# ref-based
for (i in 1:dim(Beta_Placenta_ITU_reduced_filtered)[1]) {
results_placenta_itu_lm_cpg_rb_filtered[i, 1] <- summary(lm(Beta_Placenta_ITU_reduced_filtered[i,]~ Pheno_Placenta_ITU_reduced_filtered$Trophoblasts + Pheno_Placenta_ITU_reduced_filtered$Syncytiotrophoblast + Pheno_Placenta_ITU_reduced_filtered$nRBC + Pheno_Placenta_ITU_reduced_filtered$Endothelial + Pheno_Placenta_ITU_reduced_filtered$Stromal + Pheno_Placenta_ITU_reduced_filtered$Hofbauer))$adj.r.squared
}
# note: without the outlier samples
```

``` r
# ref-free
for (i in 1:dim(Beta_Placenta_ITU_reduced_filtered)[1]) {
results_placenta_itu_lm_cpg_rf_filtered[i, 1] <- summary(lm(Beta_Placenta_ITU_reduced_filtered[i,]~ Pheno_Placenta_ITU_reduced_filtered$C1 + Pheno_Placenta_ITU_reduced_filtered$C2 + Pheno_Placenta_ITU_reduced_filtered$C3 + Pheno_Placenta_ITU_reduced_filtered$C4 + Pheno_Placenta_ITU_reduced_filtered$C5 + Pheno_Placenta_ITU_reduced_filtered$C6 + Pheno_Placenta_ITU_reduced_filtered$C7 + Pheno_Placenta_ITU_reduced_filtered$C8))$adj.r.squared
}
```

``` r
colnames(results_placenta_itu_lm_cpg_rb_filtered) <- c("adj.Rsquared")
results_placenta_itu_lm_cpg_rb_filtered <- data.frame(results_placenta_itu_lm_cpg_rb_filtered)
rownames(results_placenta_itu_lm_cpg_rb_filtered) <- rownames(Beta_Placenta_ITU_reduced_filtered)

colnames(results_placenta_itu_lm_cpg_rf_filtered) <- c("adj.Rsquared")
results_placenta_itu_lm_cpg_rf_filtered <- data.frame(results_placenta_itu_lm_cpg_rf_filtered)
rownames(results_placenta_itu_lm_cpg_rf_filtered) <- rownames(Beta_Placenta_ITU_reduced_filtered)
```

``` r
results_placenta_itu_lm_cpg_rb_filtered$adj.Rsquared[results_placenta_itu_lm_cpg_rb_filtered$adj.Rsquared <0] <- 0
results_placenta_itu_lm_cpg_rb_filtered$adj.Rsquared <- as.numeric(results_placenta_itu_lm_cpg_rb_filtered$adj.Rsquared)

results_placenta_itu_lm_cpg_rf_filtered$adj.Rsquared[results_placenta_itu_lm_cpg_rf_filtered$adj.Rsquared <0] <- 0
results_placenta_itu_lm_cpg_rf_filtered$adj.Rsquared <- as.numeric(results_placenta_itu_lm_cpg_rf_filtered$adj.Rsquared)
```

``` r
save(results_placenta_itu_lm_cpg_rb_filtered, file="Results/RData/results_placenta_itu_lm_cpg_rb_filtered.Rdata")
save(results_placenta_itu_lm_cpg_rf_filtered, file="Results/RData/results_placenta_itu_lm_cpg_rf_filtered.Rdata")
```

*check mean and sd in r2*

``` r
load("Results/RData/results_placenta_itu_lm_cpg_rb_filtered.Rdata")
load("Results/RData/results_placenta_itu_lm_cpg_rf_filtered.Rdata")
```

``` r
# mean and sd referenbece based
mean(results_placenta_itu_lm_cpg_rb_filtered$adj.Rsquared)
sd(results_placenta_itu_lm_cpg_rb_filtered$adj.Rsquared)

cpgs_hist_placenta_rb_filtered <- ggplot(data=results_placenta_itu_lm_cpg_rb_filtered, aes(adj.Rsquared)) + 
  geom_histogram()+
  geom_vline(xintercept = 0.3, color="red", linetype="dotted")+
  xlab(as.expression(bquote("adjusted" ~ R^2)))+ ylab("count (CpGs)")+
  theme_bw()+ scale_x_continuous(breaks = seq(0, 0.8, 0.2))+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title = element_text(size=10), axis.title.y=element_text(size=10))+
  labs(title="reference-based")

cpgs_hist_placenta_rb_filtered
```

``` r
# mean and sd reference-free
mean(results_placenta_itu_lm_cpg_rf_filtered$adj.Rsquared)
sd(results_placenta_itu_lm_cpg_rf_filtered$adj.Rsquared)

cpgs_hist_placenta_rf_filtered <- ggplot(data=results_placenta_itu_lm_cpg_rf_filtered, aes(adj.Rsquared)) + 
  geom_histogram()+
  geom_vline(xintercept = 0.3, color="red", linetype="dotted")+
  xlab(as.expression(bquote("adjusted" ~ R^2)))+ ylab("count (CpGs)")+
  theme_bw()+ scale_x_continuous(breaks = seq(0, 0.8, 0.2))+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title = element_text(size=10), axis.title.y=element_text(size=10))+
  labs(title="reference-free")

cpgs_hist_placenta_rf_filtered
```

``` r
M_CpGs_Placenta_ITU <- ggarrange(cpgs_hist_placenta_rb_filtered, cpgs_hist_placenta_rf_filtered, 
                                 nrow = 1,
                                 labels = c("b"),
                                 align = "hv")

M_CpGs_Placenta_ITU
```

[to the top](#top)

## Placenta PREDO

``` r
results_placenta_predo_lm_cpg_rb_filtered <-matrix(ncol=1,nrow=dim(Beta_Placenta_PREDO_filtered)[1])
results_placenta_predo_lm_cpg_rb_filtered[,1]<-rownames(Beta_Placenta_PREDO_filtered)

results_placenta_predo_lm_cpg_rf_filtered <-matrix(ncol=1,nrow=dim(Beta_Placenta_PREDO_filtered)[1])
results_placenta_predo_lm_cpg_rf_filtered[,1]<-rownames(Beta_Placenta_PREDO_filtered)
```

``` r
# ref-based
for (i in 1:dim(Beta_Placenta_PREDO_filtered)[1]) {
  results_placenta_predo_lm_cpg_rb_filtered[i, 1] <- summary(lm(Beta_Placenta_PREDO_filtered[i,]~ Pheno_Placenta_PREDO_filtered$Trophoblasts + Pheno_Placenta_PREDO_filtered$Syncytiotrophoblast + Pheno_Placenta_PREDO_filtered$nRBC + Pheno_Placenta_PREDO_filtered$Endothelial + Pheno_Placenta_PREDO_filtered$Stromal + Pheno_Placenta_PREDO_filtered$Hofbauer))$adj.r.squared
}
```

``` r
# ref-free
for (i in 1:dim(Beta_Placenta_PREDO_filtered)[1]) {
  results_placenta_predo_lm_cpg_rf_filtered[i, 1] <- summary(lm(Beta_Placenta_PREDO_filtered[i,]~ Pheno_Placenta_PREDO_filtered$C1 + Pheno_Placenta_PREDO_filtered$C2))$adj.r.squared
}
```

``` r
colnames(results_placenta_predo_lm_cpg_rb_filtered) <- c("adj.Rsquared")
results_placenta_predo_lm_cpg_rb_filtered <- data.frame(results_placenta_predo_lm_cpg_rb_filtered)
rownames(results_placenta_predo_lm_cpg_rb_filtered) <- rownames(Beta_Placenta_PREDO_filtered)

colnames(results_placenta_predo_lm_cpg_rf_filtered) <- c("adj.Rsquared")
results_placenta_predo_lm_cpg_rf_filtered <- data.frame(results_placenta_predo_lm_cpg_rf_filtered)
rownames(results_placenta_predo_lm_cpg_rf_filtered) <- rownames(Beta_Placenta_PREDO_filtered)
```

``` r
results_placenta_predo_lm_cpg_rb_filtered$adj.Rsquared[results_placenta_predo_lm_cpg_rb_filtered$adj.Rsquared <0] <- 0
results_placenta_predo_lm_cpg_rb_filtered$adj.Rsquared <- as.numeric(results_placenta_predo_lm_cpg_rb_filtered$adj.Rsquared)

results_placenta_predo_lm_cpg_rf_filtered$adj.Rsquared[results_placenta_predo_lm_cpg_rf_filtered$adj.Rsquared <0] <- 0
results_placenta_predo_lm_cpg_rf_filtered$adj.Rsquared <- as.numeric(results_placenta_predo_lm_cpg_rf_filtered$adj.Rsquared)
```

``` r
save(results_placenta_predo_lm_cpg_rb_filtered, file="Results/RData/results_placenta_predo_lm_cpg_rb_filtered.Rdata")
save(results_placenta_predo_lm_cpg_rf_filtered, file="Results/RData/results_placenta_predo_lm_cpg_rf_filtered.Rdata")
```

*check mean and sd in r2*

``` r
load("Results/RData/results_placenta_predo_lm_cpg_rb_filtered.Rdata")
load("Results/RData/results_placenta_predo_lm_cpg_rf_filtered.Rdata")
```

``` r
# mean and sd reference-based
mean(results_placenta_predo_lm_cpg_rb_filtered$adj.Rsquared)
sd(results_placenta_predo_lm_cpg_rb_filtered$adj.Rsquared)

cpgs_hist_placenta_predo_rb_filtered <- ggplot(data=results_placenta_predo_lm_cpg_rb_filtered, aes(adj.Rsquared)) + 
  geom_histogram()+
  geom_vline(xintercept = 0.3, color="red", linetype="dotted")+
  xlab(as.expression(bquote("adjusted" ~ R^2)))+ ylab("count (CpGs)")+
  theme_bw()+ scale_x_continuous(breaks = seq(0, 0.8, 0.2))+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title = element_text(size=10), axis.title.y=element_text(size=10))+
  labs(title="reference-based")

cpgs_hist_placenta_predo_rb_filtered
```

``` r
# mean and sd reference-free
mean(results_placenta_predo_lm_cpg_rf_filtered$adj.Rsquared)
sd(results_placenta_predo_lm_cpg_rf_filtered$adj.Rsquared)

cpgs_hist_placenta_predo_rf_filtered <- ggplot(data=results_placenta_predo_lm_cpg_rf_filtered, aes(adj.Rsquared)) + 
  geom_histogram()+
  geom_vline(xintercept = 0.3, color="red", linetype="dotted")+
  xlab(as.expression(bquote("adjusted" ~ R^2)))+ ylab("count (CpGs)")+
  theme_bw()+ scale_x_continuous(breaks = seq(0, 0.8, 0.2))+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title = element_text(size=10), axis.title.y=element_text(size=10))+
  labs(title="reference-free")

cpgs_hist_placenta_predo_rf_filtered
```

``` r
M_CpGs_Placenta_PREDO <- ggarrange(cpgs_hist_placenta_predo_rb_filtered, cpgs_hist_placenta_predo_rf_filtered, 
                                   nrow = 1,
                                   labels = c("c"),
                                   align = "hv")

M_CpGs_Placenta_PREDO
```

[to the top](#top)

## Placenta BET

``` r
results_placenta_BET_lm_cpg_rb_filtered <-matrix(ncol=1,nrow=dim(Beta_Placenta_BET_filtered)[1])
results_placenta_BET_lm_cpg_rb_filtered[,1]<-rownames(Beta_Placenta_BET_filtered)

results_placenta_BET_lm_cpg_rf_filtered <-matrix(ncol=1,nrow=dim(Beta_Placenta_BET_filtered)[1])
results_placenta_BET_lm_cpg_rf_filtered[,1]<-rownames(Beta_Placenta_BET_filtered)
```

``` r
for (i in 1:dim(Beta_Placenta_BET_filtered)[1]) {
  results_placenta_BET_lm_cpg_rb_filtered[i, 1] <- summary(lm(Beta_Placenta_BET_filtered[i,]~ Pheno_Placenta_BET_filtered$Trophoblasts + Pheno_Placenta_BET_filtered$Syncytiotrophoblast + Pheno_Placenta_BET_filtered$nRBC + Pheno_Placenta_BET_filtered$Endothelial + Pheno_Placenta_BET_filtered$Stromal + Pheno_Placenta_BET_filtered$Hofbauer))$adj.r.squared
}
```

``` r
for (i in 1:dim(Beta_Placenta_BET_filtered)[1]) {
  results_placenta_BET_lm_cpg_rf_filtered[i, 1] <- summary(lm(Beta_Placenta_BET_filtered[i,]~ Pheno_Placenta_BET_filtered$C1 + Pheno_Placenta_BET_filtered$C2 + Pheno_Placenta_BET_filtered$C3))$adj.r.squared
}
```

``` r
colnames(results_placenta_BET_lm_cpg_rb_filtered) <- c("adj.Rsquared")
results_placenta_BET_lm_cpg_rb_filtered <- data.frame(results_placenta_BET_lm_cpg_rb_filtered)
rownames(results_placenta_BET_lm_cpg_rb_filtered) <- rownames(Beta_Placenta_BET_filtered)

colnames(results_placenta_BET_lm_cpg_rf_filtered) <- c("adj.Rsquared")
results_placenta_BET_lm_cpg_rf_filtered <- data.frame(results_placenta_BET_lm_cpg_rf_filtered)
rownames(results_placenta_BET_lm_cpg_rf_filtered) <- rownames(Beta_Placenta_BET_filtered)
```

``` r
results_placenta_BET_lm_cpg_rb_filtered$adj.Rsquared[results_placenta_BET_lm_cpg_rb_filtered$adj.Rsquared <0] <- 0
results_placenta_BET_lm_cpg_rb_filtered$adj.Rsquared <- as.numeric(results_placenta_BET_lm_cpg_rb_filtered$adj.Rsquared)

results_placenta_BET_lm_cpg_rf_filtered$adj.Rsquared[results_placenta_BET_lm_cpg_rf_filtered$adj.Rsquared <0] <- 0
results_placenta_BET_lm_cpg_rf_filtered$adj.Rsquared <- as.numeric(results_placenta_BET_lm_cpg_rf_filtered$adj.Rsquared)
```

``` r
save(results_placenta_BET_lm_cpg_rb_filtered, file="Results/RData/results_placenta_BET_lm_cpg_rb_filtered.Rdata")
save(results_placenta_BET_lm_cpg_rf_filtered, file="Results/RData/results_placenta_BET_lm_cpg_rf_filtered.Rdata")
```

*check mean and sd in r2*

``` r
load("Results/RData/results_placenta_BET_lm_cpg_rb_filtered.Rdata")
load("Results/RData/results_placenta_BET_lm_cpg_rf_filtered.Rdata")
```

``` r
# mean and sd reference-based
mean(results_placenta_BET_lm_cpg_rb_filtered$adj.Rsquared)
sd(results_placenta_BET_lm_cpg_rb_filtered$adj.Rsquared)

cpgs_hist_placenta_BET_rb_filtered <- ggplot(data=results_placenta_BET_lm_cpg_rb_filtered, aes(adj.Rsquared)) + 
  geom_histogram()+
  geom_vline(xintercept = 0.3, color="red", linetype="dotted")+
  xlab(as.expression(bquote("adjusted" ~ R^2)))+ ylab("count (CpGs)")+
  theme_bw()+ scale_x_continuous(breaks = seq(0, 0.8, 0.2))+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title = element_text(size=10), axis.title.y=element_text(size=10))+
  labs(title="reference-based")

cpgs_hist_placenta_BET_rb_filtered
```

``` r
# mean and sd reference-free
mean(results_placenta_BET_lm_cpg_rf_filtered$adj.Rsquared)
sd(results_placenta_BET_lm_cpg_rf_filtered$adj.Rsquared)

cpgs_hist_placenta_BET_rf_filtered <- ggplot(data=results_placenta_BET_lm_cpg_rf_filtered, aes(adj.Rsquared)) + 
  geom_histogram()+
  geom_vline(xintercept = 0.3, color="red", linetype="dotted")+
  xlab(as.expression(bquote("adjusted" ~ R^2)))+ ylab("count (CpGs)")+
  theme_bw()+ scale_x_continuous(breaks = seq(0, 0.8, 0.2))+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title = element_text(size=10), axis.title.y=element_text(size=10))+
  labs(title="reference-free")

cpgs_hist_placenta_BET_rf_filtered
```

reference-based explains more

``` r
M_CpGs_Placenta_BET <- ggarrange(cpgs_hist_placenta_BET_rb_filtered, cpgs_hist_placenta_BET_rf_filtered, 
                                       nrow = 1,
                                       labels = c("d"),
                                       align = "hv")

M_CpGs_Placenta_BET
```

[to the top](#top)

## arrange Model figures

single CpG prediction models

``` r
M_CpGs_R2 <- ggarrange(M_CpGs_CVS, M_CpGs_Placenta_ITU, M_CpGs_Placenta_PREDO, M_CpGs_Placenta_BET, nrow = 4, align="hv")
M_CpGs_R2
```

``` r
ggsave("Results/CpGs_R2_combined.pdf",
       ggarrange(M_CpGs_CVS, M_CpGs_Placenta_ITU, M_CpGs_Placenta_PREDO, M_CpGs_Placenta_BET, nrow = 4, align="hv"), width=84, height=180, units="mm", dpi=600, scale=2)
```

[to the top](#top)

# CpGs most influenced by cell types (reference-based)

## extract CpGs

load R2 (predicted by RPC cell types)

``` r
# reference-based
load("Results/RData/results_cvs_itu_lm_cpg_rb_filtered.Rdata")
load("Results/RData/results_placenta_itu_lm_cpg_rb_filtered.Rdata")
load("Results/RData/results_placenta_predo_lm_cpg_rb_filtered.Rdata")
load("Results/RData/results_placenta_BET_lm_cpg_rb_filtered.Rdata")
```

``` r
q9_cvs_rb <- quantile(results_cvs_itu_lm_cpg_rb_filtered$adj.Rsquared, probs=0.9)
q9_placenta_itu_rb <- quantile(results_placenta_itu_lm_cpg_rb_filtered$adj.Rsquared, probs=0.9)
q9_placenta_predo_rb <- quantile(results_placenta_predo_lm_cpg_rb_filtered$adj.Rsquared, probs=0.9)
q9_placenta_bet_rb <- quantile(results_placenta_BET_lm_cpg_rb_filtered$adj.Rsquared, probs=0.9)
```

## genes of all CpGs that overlap between data sets

*mapping CpGs to genes*

``` r
# load matching EPIC CpGs Genes
genes_cpgs_map <- read.table("Input_Data_prepared/genes_CpGs_EPIC.txt", header=T)
# load gene information
load("Input_Data_prepared/genes_CpGs_EPIC.Rdata") #epic_genes
```

``` r
identical(genes_cpgs_map$gene, as.character(epic_genes$name))
genes_cpgs_fullinfo <- cbind(genes_cpgs_map, epic_genes)
```

``` r
write.table(genes_cpgs_fullinfo, file="Input_Data_prepared/genes_cpgs_fullinfo.txt")
save(genes_cpgs_fullinfo, file="Input_Data_prepared/genes_cpgs_fullinfo.Rdata")
```

``` r
cpgs_all_filtered <- list(rownames(results_cvs_itu_lm_cpg_rb_filtered), rownames(results_placenta_itu_lm_cpg_rb_filtered), rownames(results_placenta_predo_lm_cpg_rb_filtered), rownames(results_placenta_predo_lm_cpg_rb_filtered))
# note: although I use the reference-based results here, I only extract CpG names which are the same for ref-based and ref-free - so this is the same for ref-based and ref-free and only for simplicity and object loading 
```

``` r
# which are the cpgs that overlap between the data sets
overlap_cpgs_all_filtered <- Reduce(intersect,cpgs_all_filtered)
length(overlap_cpgs_all_filtered)
```

``` r
overlap_cpgs_genes_all_filtered <- genes_cpgs_fullinfo[genes_cpgs_fullinfo$CpG %in% overlap_cpgs_all_filtered, ]
dim(overlap_cpgs_genes_all_filtered)
```

``` r
genes_all_Symbol_filtered <- as.character(unique(overlap_cpgs_genes_all_filtered$gene)) #

length(genes_all_Symbol_filtered)
```

``` r
write.table(genes_all_Symbol_filtered, file="Results/RData/genes_all_Symbol_filtered.txt", row.names = FALSE,col.names = FALSE)
save(genes_all_Symbol_filtered, file="Results/RData/genes_all_Symbol_filtered.Rdata")
```

## extract CpGs with R2 \> 30% predicted by reference-based cell types that are in all data sets

``` r
cpgs_cvs_itu_lmrb_filtered <- subset(results_cvs_itu_lm_cpg_rb_filtered, adj.Rsquared > 0.3)
dim(cpgs_cvs_itu_lmrb_filtered)
cpgs_placenta_itu_lmrb_filtered <- subset(results_placenta_itu_lm_cpg_rb_filtered, adj.Rsquared > 0.3)
dim(cpgs_placenta_itu_lmrb_filtered)
cpgs_placenta_predo_lmrb_filtered <- subset(results_placenta_predo_lm_cpg_rb_filtered, adj.Rsquared > 0.3)
dim(cpgs_placenta_predo_lmrb_filtered)
cpgs_placenta_BET_lmrb_filtered <- subset(results_placenta_BET_lm_cpg_rb_filtered, adj.Rsquared > 0.3)
dim(cpgs_placenta_BET_lmrb_filtered)


cpgs_list_filtered <- list(rownames(cpgs_cvs_itu_lmrb_filtered), rownames(cpgs_placenta_itu_lmrb_filtered), rownames(cpgs_placenta_predo_lmrb_filtered), rownames(cpgs_placenta_BET_lmrb_filtered))
```

``` r
# which are the cpgs that overlap between the data sets with R2 > .30
overlap_cpgs_filtered <- Reduce(intersect,cpgs_list_filtered)
length(overlap_cpgs_filtered)
# 26.092
```

``` r
overlap_cpgs_genes_filtered <- genes_cpgs_fullinfo[genes_cpgs_fullinfo$CpG %in% overlap_cpgs_filtered, ]
dim(overlap_cpgs_genes_filtered)
```

``` r
write.table(overlap_cpgs_genes_filtered, file="Results/RData/overlap_cpgs_genes_rb_filtered.txt", row.names = FALSE,col.names = FALSE)
save(overlap_cpgs_genes_filtered, file="Results/RData/overlap_cpgs_genes_rb_filtered.Rdata")
# -> CpGs mapped to genes with R2 > .30 and in all data sets
```

``` r
genes_overlap_Symbol_filtered <- as.character(unique(overlap_cpgs_genes_filtered$gene)) 
length(genes_overlap_Symbol_filtered)
```

``` r
write.table(genes_overlap_Symbol_filtered, file="Results/RData/genes_overlap_Symbol_filtered.txt", row.names = FALSE,col.names = FALSE)
save(genes_overlap_Symbol_filtered, file="Results/RData/genes_overlap_Symbol_filtered.Rdata")
```

extract the CpGs with gene info etc. that are used as input for
enrichment (CpGs influenced by cell types)

``` r
identical(overlap_cpgs_genes_filtered$gene, as.character(overlap_cpgs_genes_filtered$name))
overlap_cpgs_genes_filtered$GeneSymbol <- overlap_cpgs_genes_filtered$gene
overlap_cpgs_genes_filtered$gene <- NULL
overlap_cpgs_genes_filtered$name <- NULL
```

``` r
CpGs_R30_Genes_refbased <- overlap_cpgs_genes_filtered[ ,c("CpG", "GeneSymbol")]
```

*Table S2 - reference-based*

``` r
write.xlsx2(CpGs_R30_Genes_refbased, "Results/List_CpGs_R30_Genes_refbased.xlsx", sheetName = "List", col.names = TRUE, row.names = FALSE, append = FALSE)
```

[to the top](#top)

## Tissue-specific Gene Enrichment

The TissueEnrich R package is used to calculate enrichment of
tissue-specific genes in a set of input genes.

The user can input genes to determine which tissue-specific genes are
enriched in those datasets.

The hypergeometric test is being used to determine if the
tissue-specific genes are enriched among the input genes. Function
teEnrichment: Given a gene list as input, this function calculates the
tissue-specific gene enrichment using tissue-specific genes from either
human or mouse RNA-Seq datasets.

### genes with R2 \> .30 that overlap between the data sets against all genes

``` r
load("Results/RData/genes_all_Symbol_filtered.Rdata") 
length(genes_all_Symbol_filtered)
```

``` r
genes_overlap_input_filtered <-GeneSet(geneIds=genes_overlap_Symbol_filtered,organism="Homo Sapiens",geneIdType=SymbolIdentifier()) 

genes_overlap_background_filtered <-GeneSet(geneIds=genes_all_Symbol_filtered,organism="Homo Sapiens",geneIdType=SymbolIdentifier())
```

When using teEnrichment, the user can specify the RNA-Seq dataset
(rnaSeqDataset) to be used for the tissue-specific gene enrichment
analysis. 1 for “Human Protein Atlas” (default) 2 for “GTEx” 3 for
“Mouse ENCODE”

``` r
output_GEnrich_filtered <- teEnrichment(inputGenes = genes_overlap_input_filtered, backgroundGenes = genes_overlap_background_filtered, rnaSeqDataset = 1)
```

The first object is the SummarizedExperiment object containing the
enrichment results, the second and the third object contains the
expression values and tissue-specificity information of the
tissue-specific genes for genes from the input gene set, and the fourth
is a GeneSet object containing genes that were not identified in the
tissue-specific gene data.

``` r
seEnrichmentOutput_filtered <-output_GEnrich_filtered[[1]]

enrichmentOutput_filtered<-setNames(data.frame(assay(seEnrichmentOutput_filtered),row.names = rowData(seEnrichmentOutput_filtered)[,1]), colData(seEnrichmentOutput_filtered)[,1])
enrichmentOutput_filtered$Tissue<-row.names(enrichmentOutput_filtered)
```

``` r
save(output_GEnrich_filtered, file="Results/RData/output_GEnrich_filtered.Rdata")
save(enrichmentOutput_filtered, file="Results/RData/enrichmentOutput_filtered.Rdata")
```

``` r
load("Results/RData/output_GEnrich_filtered.Rdata") 
# list object containing the enrichment results
load("Results/RData/enrichmentOutput_filtered.Rdata")
# containing the −Log10(P−Value) and fold-change, corresponding to the tissue-specific gene enrichment, along with the number of tissue-specific genes in the input gene set. This object can be used to visualize tissue-specific gene enrichment in the form of a bar chart using the −Log10(P−Value) values.
```

all genes mapped:

``` r
sum(enrichmentOutput_filtered$Tissue.Specific.Genes)
```

genes not identified in tissue-specific gene data

``` r
genes_not_ident_TissueEnrich_filtered <- geneIds(output_GEnrich_filtered[[4]]) #865
```

In TissueEnrich, they are dividing the genes into six different groups
which are specified in their paper
(<http://doi.org/10.1093/bioinformatics/bty890>). The missing genes
could be in the other three non-tissue specific gene groups (Not
Expressed, Expressed In all, or Mixed).

Tissue-specific genes are defined using the algorithm from the HPA
(Uhlén et al. 2015), and can be grouped as follows: Tissue Enriched:
Genes with an expression level greater than 1 (TPM or FPKM) that also
have at least five-fold higher expression levels in a particular tissue
compared to all other tissues. Group Enriched: Genes with an expression
level greater than 1 (TPM or FPKM) that also have at least five-fold
higher expression levels in a group of 2-7 tissues compared to all other
tissues, and that are not considered Tissue Enriched. Tissue Enhanced:
Genes with an expression level greater than 1 (TPM or FPKM) that also
have at least five-fold higher expression levels in a particular tissue
compared to the average levels in all other tissues, and that are not
considered Tissue Enriched or Group Enriched.

``` r
tissue_enrich_main_plot_refbased <- ggplot(enrichmentOutput_filtered,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-log10 (p-adjusted)')+
  scale_y_continuous(breaks = c(0, 2, 4, 6), limits = c(0,8), expand = c(0,0))+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank(), plot.margin=unit(c(1,1,1,2),"cm"))+
  geom_hline(yintercept=2, linetype="dashed", color = "black")

ggsave("Results/TissueEnrichment_CpGGenes.pdf",
       tissue_enrich_main_plot_refbased, width=84, height=50, units="mm", dpi=600, scale=2, device = cairo_pdf)

ggsave(tissue_enrich_main_plot_refbased, filename = "Results/TissueEnrichment_CpGGenes.png", dpi = 300, type = "cairo",
       width = 8, height = 5, units = "in")
```

generate a heatmap showing the expression of the placenta specific genes
across all the tissues.

``` r
seExp_filtered<-output_GEnrich_filtered[[2]][["Placenta"]]
exp_filtered<-setNames(data.frame(assay(seExp_filtered), row.names = rowData(seExp_filtered)[,1]), colData(seExp_filtered)[,1])
exp_filtered$Gene<-row.names(exp_filtered)
exp_filtered<-exp_filtered %>% gather(key = "Tissue", value = "expression",1:(ncol(exp_filtered)-1))

ggplot(exp_filtered, aes(Tissue, Gene)) + geom_tile(aes(fill = expression),
                                                    colour = "white") + scale_fill_gradient(low = "white",
                                                                                            high = "steelblue")+
  labs(x='', y = '')+
  theme_bw()+
  guides(fill = guide_legend(title = "Log2(TPM)"))+
  #theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1, size=5), panel.grid.major= element_blank(),panel.grid.minor = element_blank())
```

get p value for placenta

``` r
enrichmentOutput_filtered["Placenta", ]
options(scipen=999)
log10p_placenta_enrich_filtered <- enrichmentOutput_filtered["Placenta", "Log10PValue"]
p_placenta_enrich_filtered <- 10^-log10p_placenta_enrich_filtered

log10p_placenta_enrich_filtered
p_placenta_enrich_filtered
```

``` r
#list containing the tissue-specificity information for the input genes. The code below retrieves the tissue-specific genes along with the type of tissue-specificity in placenta tissue.
seGroupInf_Placenta_filtered <-output_GEnrich_filtered[[3]][["Placenta"]]
groupInf_Placenta_filtered <-data.frame(assay(seGroupInf_Placenta_filtered)) # 186
groupInf_Placenta_filtered
```

``` r
save(groupInf_Placenta_filtered, file="Results/RData/groupInf_Placenta_filtered.Rdata")
write.table(groupInf_Placenta_filtered$Gene, file="Results/RData/genes_cpgs_placenta_enriched_filtered.txt", row.names = F, col.names = F)
# 186 genes highly important / specific for placenta 
```

export the 186 placenta-specific genes in our input data to excel

*Table S3*

``` r
write.xlsx(groupInf_Placenta_filtered,"Results/genes_cpgs_placenta_enriched_filtered.xlsx", row.names = F, col.names = F)
# 186 genes highly important / specific for placenta 
```

[to the top](#top)

## Placenta cell enrichment

PlacentaCellEnrich enables users to provide background genes for
carrying out cell-specific gene enrichment. In this case, instead of
using all the genes in the dataset, a background gene set is being used
to carry out the enrichment analysis. It should be noted that the
background gene set must have all the genes of the input gene set. The
p-values are corrected for multiple hypothesis testing using the
Benjamini & Hochberg correction. If the background gene set is not
provided all the genes will be used as background.

How should I interpret my results? The recommended threshold value of
the adjusted p-values and fold-change for selecting the enriched cells
is 0.01 and 2 respectively. It is also recommended that the users should
look at cell-specific gene enrichment from all the sources and have the
highest confidence in results that are consistent across datasets.

INPUT:

-   gene set: genes with R2 > .30 that overlap between data sets
    (genes_overlap_Sympbol_filtered)  
-   background: genes that overlap between data sets
    (genes_all_Symbol_filtered)

``` r
cell_enrichment_results_refbased <- read_excel("Results/PlacentaCellEnrich_Tool/refbased/CellSpecificGeneEnrichment-VentoTormo.xlsx", col_names = T)
cell_enrichment_results_refbased$`-Log10PValue` <- as.numeric(cell_enrichment_results_refbased$`-Log10PValue`)
```

``` r
cell_enrichment_colors_cells <- cell_enrichment_results_refbased[,"Cell"]
```

``` r
cell_enrichment_results_refbased <- merge(cell_enrichment_results_refbased, cell_enrichment_colors_cells, by="Cell")
```

``` r
plot_cell_enrichment_refbased <- ggplot2::ggplot(cell_enrichment_results_refbased, aes(x = reorder(Cell, -`-Log10PValue`, sum), y= `-Log10PValue`, fill=Cell)) +
geom_bar(stat = 'identity')+ 
xlab("") + ylab("-log10 (p-adjusted)")+
scale_y_continuous(breaks = round(seq(0, 14, by = 2),1),limits = c(0,16), expand = c(0,0))+
theme_bw()+
theme(legend.position="none")+
theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title = element_text(size=10))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank(), plot.margin=unit(c(1,1,1,2),"cm"))+
geom_hline(yintercept=2, linetype="dashed", color = "black")
```

[to the top](#top)

# CpGs most influenced by cell types (reference-free)

## extract CpGs

load R2 (predicted by RPC cell types)

``` r
load("Results/RData/results_cvs_itu_lm_cpg_rf_filtered.Rdata")
load("Results/RData/results_placenta_itu_lm_cpg_rf_filtered.Rdata")
load("Results/RData/results_placenta_predo_lm_cpg_rf_filtered.Rdata")
load("Results/RData/results_placenta_BET_lm_cpg_rf_filtered.Rdata")

load("Input_Data_prepared/genes_cpgs_fullinfo.Rdata")
```

``` r
q9_cvs_rf <- quantile(results_cvs_itu_lm_cpg_rf_filtered$adj.Rsquared, probs=0.9)
q9_placenta_itu_rf <- quantile(results_placenta_itu_lm_cpg_rf_filtered$adj.Rsquared, probs=0.9)
q9_placenta_predo_rf <- quantile(results_placenta_predo_lm_cpg_rf_filtered$adj.Rsquared, probs=0.9)
q9_placenta_bet_rf <- quantile(results_placenta_BET_lm_cpg_rf_filtered$adj.Rsquared, probs=0.9)
```

mean of 90% quantiles in adjusted R2 of CpGs (refbased & reffree)

``` r
c9_rb_rf_r <- as.numeric(c(q9_cvs_rb, q9_placenta_itu_rb, q9_placenta_predo_rb, q9_placenta_bet_rb, q9_cvs_rf, q9_placenta_itu_rf, q9_placenta_predo_rf, q9_placenta_bet_rf))
mean(c9_rb_rf_r)
```

## extract CpGs with R2 \> 30% that are in all data sets (ref-free)

``` r
cpgs_cvs_itu_lmrf_filtered <- subset(results_cvs_itu_lm_cpg_rf_filtered, adj.Rsquared > 0.3)
dim(cpgs_cvs_itu_lmrf_filtered)
cpgs_placenta_itu_lmrf_filtered <- subset(results_placenta_itu_lm_cpg_rf_filtered, adj.Rsquared > 0.3)
dim(cpgs_placenta_itu_lmrf_filtered)
cpgs_placenta_predo_lmrf_filtered <- subset(results_placenta_predo_lm_cpg_rf_filtered, adj.Rsquared > 0.3)
dim(cpgs_placenta_predo_lmrf_filtered)
cpgs_placenta_BET_lmrf_filtered <- subset(results_placenta_BET_lm_cpg_rf_filtered, adj.Rsquared > 0.3)
dim(cpgs_placenta_BET_lmrf_filtered)
    
cpgs_list_rf_filtered <- list(rownames(cpgs_cvs_itu_lmrf_filtered), rownames(cpgs_placenta_itu_lmrf_filtered), rownames(cpgs_placenta_predo_lmrf_filtered), rownames(cpgs_placenta_BET_lmrf_filtered))
```

``` r
# which are the cpgs that overlap between the data sets
overlap_cpgs_rf_filtered <- Reduce(intersect,cpgs_list_rf_filtered)
length(overlap_cpgs_rf_filtered)
# 531
```

``` r
overlap_cpgs_genes_rf_filtered <- genes_cpgs_fullinfo[genes_cpgs_fullinfo$CpG %in% overlap_cpgs_rf_filtered, ]
dim(overlap_cpgs_genes_rf_filtered)
```

``` r
write.table(overlap_cpgs_genes_rf_filtered, file="Results/RData/overlap_cpgs_genes_rf_filtered.txt", row.names = FALSE,col.names = FALSE)
save(overlap_cpgs_genes_rf_filtered, file="Results/RData/overlap_cpgs_genes_rf_filtered.Rdata")
    # -> CpGs mapped to genes with R2 > .30 and in all data sets
```

``` r
genes_overlap_Symbol_rf_filtered <- as.character(unique(overlap_cpgs_genes_rf_filtered$gene)) 
length(genes_overlap_Symbol_rf_filtered)
```

``` r
write.table(genes_overlap_Symbol_rf_filtered, file="Results/RData/genes_overlap_Symbol_rf_filtered.txt", row.names = FALSE,col.names = FALSE)
    
save(genes_overlap_Symbol_rf_filtered, file="Results/RData/genes_overlap_Symbol_rf_filtered.Rdata")
```

[to the top](#top)

extract the CpGs with gene info etc. that are used as input for
enrichment (CpGs influenced by cell types)

``` r
identical(overlap_cpgs_genes_rf_filtered$gene, as.character(overlap_cpgs_genes_rf_filtered$name))
overlap_cpgs_genes_rf_filtered$GeneSymbol <- overlap_cpgs_genes_rf_filtered$gene
overlap_cpgs_genes_rf_filtered$gene <- NULL
overlap_cpgs_genes_rf_filtered$name <- NULL
```

``` r
CpGs_R30_Genes_reffree <- overlap_cpgs_genes_rf_filtered[ ,c("CpG", "GeneSymbol")]
```

``` r
write.xlsx2(CpGs_R30_Genes_reffree, "Results/List_CpGs_R30_Genes_reffree.xlsx", sheetName = "List", col.names = TRUE, row.names = FALSE, append = FALSE)
```

## Tissue-specific Gene Enrichment

### genes with R2 \> .30 that overlap between the data sets against all genes

``` r
load("Results/RData/genes_all_Symbol_filtered.Rdata")
load("Results/RData/genes_overlap_Symbol_rf_filtered.Rdata") 
```

``` r
genes_overlap_input_rf_filtered <-GeneSet(geneIds=genes_overlap_Symbol_rf_filtered,organism="Homo Sapiens",geneIdType=SymbolIdentifier()) 
genes_overlap_background_filtered <-GeneSet(geneIds=genes_all_Symbol_filtered,organism="Homo Sapiens",geneIdType=SymbolIdentifier())

output_GEnrich_rf_filtered <- teEnrichment(inputGenes = genes_overlap_input_rf_filtered, backgroundGenes = genes_overlap_background_filtered, rnaSeqDataset = 1)
```

``` r
seEnrichmentOutput_rf_filtered<-output_GEnrich_rf_filtered[[1]]
enrichmentOutput_rf_filtered<-setNames(data.frame(assay(seEnrichmentOutput_rf_filtered),row.names = rowData(seEnrichmentOutput_rf_filtered)[,1]), colData(seEnrichmentOutput_rf_filtered)[,1])
enrichmentOutput_rf_filtered$Tissue<-row.names(enrichmentOutput_rf_filtered)
```

``` r
save(output_GEnrich_rf_filtered, file="Results/RData/output_GEnrich_rf_filtered.Rdata")
save(enrichmentOutput_rf_filtered, file="Results/RData/enrichmentOutput_rf_filtered.Rdata")
```

``` r
load("Results/RData/output_GEnrich_rf_filtered.Rdata")
# ist object containing the enrichment results
load("Results/RData/enrichmentOutput_rf_filtered.Rdata")
# SummarizedExperiment object containing the −Log10(P−Value) and fold-change, corresponding to the tissue-specific gene enrichment, along with the number of tissue-specific genes in the input gene set. This object can be used to visualize tissue-specific gene enrichment in the form of a bar chart using the −Log10(P−Value) values.
```

all genes mapped:

``` r
sum(enrichmentOutput_rf_filtered$Tissue.Specific.Genes)
```

genes not identified in tissue-specific gene data

``` r
genes_not_ident_TissueEnrich_rf_filtered <- geneIds(output_GEnrich_rf_filtered[[4]])
```

``` r
tissueEnrichment_reffree_filtered <- ggplot(enrichmentOutput_rf_filtered,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-log10 (p-adjusted)')+
  scale_y_continuous(breaks = c(0, 2, 4, 6), limits = c(0,8), expand = c(0,0))+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=10),panel.grid.major= element_blank(),panel.grid.minor = element_blank(), plot.margin=unit(c(1,1,1,2),"cm"))+
  geom_hline(yintercept=2, linetype="dashed", color = "black")

ggsave("Results/TissueEnrichment_CpGGenes_reffree.pdf",
       tissueEnrichment_reffree_filtered, width=84, height=50, units="mm", dpi=600, scale=2, device = cairo_pdf)

ggsave(tissueEnrichment_reffree_filtered, filename = "Results/TissueEnrichment_CpGGenes_reffree.png", dpi = 300, type = "cairo",
       width = 8, height = 5, units = "in")
```

get p value for placenta

``` r
enrichmentOutput_rf_filtered["Placenta", ]
options(scipen=999)
log10p_placenta_enrich_rf_filtered <- enrichmentOutput_rf_filtered["Placenta", "Log10PValue"]
p_placenta_enrich_rf_filtered <- 10^-log10p_placenta_enrich_rf_filtered
```

see if tissue is significant (p \< .01 / -log10 p \>2)

``` r
enrichmentOutput_rf_filtered[(enrichmentOutput_rf_filtered[,1]>2),]
```

``` r
#list containing the tissue-specificity information for the input genes. The code below retrieves the tissue-specific genes along with the type of tissue-specificity in placenta tissue.
seGroupInf_Placenta_rf_filtered <-output_GEnrich_rf_filtered[[3]][["Placenta"]]
groupInf_Placenta_rf_filtered <-data.frame(assay(seGroupInf_Placenta_rf_filtered)) # 10
groupInf_Placenta_rf_filtered
```

``` r
save(groupInf_Placenta_rf_filtered, file="Results/RData/genes_groupInf_Placenta_rf_filtered.Rdata")
write.table(groupInf_Placenta_rf_filtered$Gene, file="Results/RData/genes_cpgs_placenta_enriched_rf_filtered.txt", row.names = F, col.names = F)
#10 genes highly important / specific for placenta: genes_cpgs_placenta_enriched_rf_filtered
```

``` r
#list of input genes that were not identified in the tissue-specific gene data.
unidentified_rf_filtered <- geneIds(output_GEnrich_rf_filtered[[4]]) 
```

## Placenta cell enrichment

INPUT:

-   gene set: genes with R2 > .30 that overlap between data sets
    (genes_overlap_Symbol_rf_filtered)  
-   background: genes that overlap between data sets
    (genes_all_Symbol_filtered)

``` r
cell_enrichment_results_reffree <- read_excel("Results/PlacentaCellEnrich_Tool/reffree/CellSpecificGeneEnrichment-VentoTormo.xlsx", col_names = T)
cell_enrichment_results_reffree$`-Log10PValue` <- as.numeric(cell_enrichment_results_reffree$`-Log10PValue`)
```

``` r
cell_enrichment_results_reffree <- merge(cell_enrichment_results_reffree, cell_enrichment_colors_cells, by="Cell")
```

``` r
plot_cell_enrichment_reffree <-
ggplot2::ggplot(cell_enrichment_results_reffree, aes(x = reorder(Cell, -`-Log10PValue`, sum), y= `-Log10PValue`, fill=Cell)) +
geom_bar(stat = 'identity')+ 
xlab("") + ylab("-log10 (p-adjusted)")+
scale_y_continuous(breaks = round(seq(0, 14, by = 2),1),limits = c(0,16), expand = c(0,0))+
theme_bw()+
theme(legend.position="none")+
theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title = element_text(size=10))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank(), plot.margin=unit(c(1,1,1,2),"cm"))+
geom_hline(yintercept=2, linetype="dashed", color = "black")
```

[to the top](#top)

# save Tissue Enrichment Plot

*Fig. 4*

``` r
ggsave("Results/Plots_tissue_Enrichment_combined.pdf",
       ggarrange(tissue_enrich_main_plot_refbased, tissueEnrichment_reffree_filtered, ncol=1, labels = c("a", "b")), width=84, height=160, units="mm", dpi=600, scale=2)
```

# save Cell Enrichment Plot

*Fig. 5*

``` r
ggsave("Results/Plots_Cell_Enrichment_combined.pdf",
       ggarrange(plot_cell_enrichment_refbased, plot_cell_enrichment_reffree, ncol=1, labels = c("a", "b")), width=84, height=160, units="mm", dpi=600, scale=2)
```

[to the top](#top)

# Plot Cell types (RPC, reference-based)

*Fig. 6*

RPC reference-based

``` r
CVS_ITU_RPC_Plot <- 
    as_tibble(Pheno_CVS_ITU_filtered) %>%
  pivot_longer(cols = all_of(names_refbased_cells),
               names_to = 'component',
               values_to = 'estimate') %>%
  ggplot(aes(x = Sample_Name, y = estimate, fill = component)) +
  geom_bar(stat = 'identity', width=1) +
  scale_fill_manual(values = plColors)+
  theme_minimal() +
  scale_y_continuous(limits = c(-0.1,1.1), breaks = c(0, 0.5, 1), labels = scales::percent) +
  theme(axis.text.x = element_blank(), 
        axis.text.y=element_text(size=10), 
        axis.title.y = element_text(size=10), 
        axis.title.x=element_blank(), 
        legend.text = element_text(size=10), 
        legend.key.size = unit(0.3, "cm"), 
        legend.margin=margin(t = 0, unit='cm'), 
        legend.title = element_blank(), 
        panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0,1))

Placenta_ITU_RPC_Plot <-
    as_tibble(Pheno_Placenta_ITU_reduced_filtered) %>%
  pivot_longer(cols = all_of(names_refbased_cells),
               names_to = 'component',
               values_to = 'estimate') %>%
  ggplot(aes(x = Sample_Name, y = estimate, fill = component)) +
  geom_bar(stat = 'identity', width=1) +
  scale_fill_manual(values = plColors)+
  theme_minimal() +
  scale_y_continuous(limits = c(-0.1,1.1), breaks = c(0, 0.5, 1), labels = scales::percent) +
  theme(axis.text.x = element_blank(), 
        axis.text.y=element_text(size=10), 
        axis.title.y = element_text(size=10), 
        axis.title.x=element_blank(), 
        legend.text = element_text(size=5), 
        legend.key.size = unit(0.3, "cm"), 
        legend.margin=margin(t = 0, unit='cm'), 
        legend.title = element_blank(), 
        panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0,1))

Placenta_PREDO_RPC_Plot <- 
    as_tibble(Pheno_Placenta_PREDO_filtered) %>%
  pivot_longer(cols = all_of(names_refbased_cells),
               names_to = 'component',
               values_to = 'estimate') %>%
  ggplot(aes(x = ID, y = estimate, fill = component)) +
  geom_bar(stat = 'identity', width=1) +
  scale_fill_manual(values = plColors)+
  theme_minimal() +
  scale_y_continuous(limits = c(-0.1,1.1), breaks = c(0, 0.5, 1), labels = scales::percent) +
  theme(axis.text.x = element_blank(), 
        axis.text.y=element_text(size=5), 
        axis.title.y = element_text(size=10), 
        axis.title.x=element_blank(), 
        legend.text = element_text(size=10), 
        legend.key.size = unit(0.3, "cm"), 
        legend.margin=margin(t = 0, unit='cm'), 
        legend.title = element_blank(), 
        panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0,1))

Placenta_BET_RPC_Plot <- 
    as_tibble(Pheno_Placenta_BET_filtered) %>%
  pivot_longer(cols = all_of(names_refbased_cells),
               names_to = 'component',
               values_to = 'estimate') %>%
  ggplot(aes(x = Sample_Name, y = estimate, fill = component)) +
  geom_bar(stat = 'identity', width=1) +
  scale_fill_manual(values = plColors)+
  theme_minimal() +
  scale_y_continuous(limits = c(-0.1,1.1), breaks = c(0, 0.5, 1), labels = scales::percent) +
  theme(axis.text.x = element_blank(), 
        axis.text.y=element_text(size=5), 
        axis.title.y = element_text(size=10), 
        axis.title.x=element_blank(), 
        legend.text = element_text(size=10), 
        legend.key.size = unit(0.3, "cm"), 
        legend.margin=margin(t = 0, unit='cm'), 
        legend.title = element_blank(), 
        panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0,1))
```

plot

``` r
CVS_ITU_RPC_Plot
Placenta_ITU_RPC_Plot
Placenta_PREDO_RPC_Plot
Placenta_BET_RPC_Plot
```

``` r
CVS_ITU_RPC_Plot_arrange <- 
    as_tibble(Pheno_CVS_ITU_filtered) %>%
  pivot_longer(cols = all_of(names_refbased_cells),
               names_to = 'component',
               values_to = 'estimate') %>%
  ggplot(aes(x = Sample_Name, y = estimate, fill = component)) +
  labs(tag = "")+
  geom_bar(stat = 'identity', width=1) +
  scale_fill_manual(values = plColors)+
  theme_minimal() +
  scale_y_continuous(limits = c(0,1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = scales::percent) +
  theme(axis.text.x = element_blank(), 
        axis.text.y=element_text(size=10), 
        axis.title.y = element_blank(), 
        axis.title.x=element_blank(), 
        legend.text = element_text(size=10), 
        legend.title = element_blank(),legend.position="bottom",
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  plot.title=element_text(size=10, hjust=0.5)) +
  coord_cartesian(ylim = c(0,1.1))


Placenta_ITU_RPC_Plot_arrange <-
    as_tibble(Pheno_Placenta_ITU_reduced_filtered) %>%
  pivot_longer(cols = all_of(names_refbased_cells),
               names_to = 'component',
               values_to = 'estimate') %>%
  ggplot(aes(x = Sample_Name, y = estimate, fill = component)) +
  labs(tag = "")+
  geom_bar(stat = 'identity', width=1) +
  scale_fill_manual(values = plColors)+
  theme_minimal() +
  scale_y_continuous(limits = c(0,1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = scales::percent) +
  theme(axis.text.x = element_blank(), 
        axis.text.y=element_text(size=9), 
        axis.title.y = element_blank(), 
        axis.title.x=element_blank(), 
        legend.text = element_text(size=10), 
        legend.title = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), plot.title=element_text(size=10, hjust=0.5)) +
  coord_cartesian(ylim = c(0,1.1))

Placenta_PREDO_RPC_Plot_arrange <- 
    as_tibble(Pheno_Placenta_PREDO_filtered) %>%
  pivot_longer(cols = all_of(names_refbased_cells),
               names_to = 'component',
               values_to = 'estimate') %>%
  ggplot(aes(x = ID, y = estimate, fill = component)) +
  labs(tag = "")+
  geom_bar(stat = 'identity', width=1) +
  scale_fill_manual(values = plColors)+
  theme_minimal() +
  scale_y_continuous(limits = c(0,1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = scales::percent) +
  theme(axis.text.x = element_blank(), 
        axis.text.y=element_text(size=10), 
        axis.title.y = element_blank(), 
        axis.title.x=element_blank(), 
        legend.text = element_text(size=10), 
        legend.title = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), plot.title=element_text(size=10, hjust=0.5)) +
  coord_cartesian(ylim = c(0,1.1))

Placenta_BET_RPC_Plot_arrange <- 
    as_tibble(Pheno_Placenta_BET_filtered) %>%
  pivot_longer(cols =all_of(names_refbased_cells),
               names_to = 'component',
               values_to = 'estimate') %>%
  ggplot(aes(x = Sample_Name, y = estimate, fill = component)) +
  labs(tag = "")+
  geom_bar(stat = 'identity', width=1) +
  scale_fill_manual(values = plColors)+
  theme_minimal() +
  scale_y_continuous(limits = c(0,1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = scales::percent) +
  theme(axis.text.x = element_blank(), 
        axis.text.y=element_text(size=10), 
        axis.title.y = element_blank(), 
        axis.title.x=element_blank(), 
        legend.text = element_text(size=10), 
        legend.title = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), plot.title=element_text(size=10, hjust=0.5)) +
  coord_cartesian(ylim = c(0,1.1))
```

``` r
RPC_Plot_Combined <- 
ggarrange(
          CVS_ITU_RPC_Plot_arrange +
           theme(legend.position="none", plot.margin = margin(r = 0.2, t= -0.2) ),
          Placenta_ITU_RPC_Plot_arrange +
               theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_blank(), plot.margin = margin(r = 0.2, l = 0.2)),
          Placenta_PREDO_RPC_Plot_arrange +
               theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_blank(), plot.margin = margin(l = 0.2)),
          Placenta_BET_RPC_Plot_arrange +
               theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_blank(), plot.margin = margin(l = 0.2)),
          common.legend = T,
          labels = c("a", "b", "c","d", size=10),
          legend="bottom",
          nrow = 1,
          align = "hv")

RPC_Plot_Combined
```

``` r
ggsave("Results/RPC_Cells_Combined.pdf",
RPC_Plot_Combined, width=84, height=40, units="mm", dpi=600, scale=2, device = cairo_pdf)

ggsave(RPC_Plot_Combined, filename = "Results/RPC_Cells_Combined.png", dpi = 300, type = "cairo",
       width = 6, height = 3, units = "in")
```

plots with mean and sd

``` r
cells_cvs <- data.frame(psych::describe(Pheno_CVS_ITU_filtered[ ,names_refbased_cells]))
cells_cvs <- cells_cvs[ ,c("mean", "sd")]
rownames(cells_cvs) <- names_refbased_cells

msd_plot_cells_cvs <- ggplot(cells_cvs, aes(x=as.factor(rownames(cells_cvs)), y=mean, fill = as.factor(rownames(cells_cvs)))) +
  geom_bar(position=position_dodge(), stat="identity") +
  scale_fill_manual(values = plColors)+
  theme_minimal()+
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,position=position_dodge(.9))+
  labs(tag = "a", title="CVS (ITU)")+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.y = element_text(size=10), 
        legend.text = element_text(size=10), 
        legend.title = element_blank(), 
        plot.title=element_text(size=10))+
  scale_y_continuous(limits = c(0,1.1), breaks = c(0, 0.2,0.4, 0.6, 0.8, 1), labels = scales::percent)+
  guides(fill = guide_legend(label.position = "bottom", nrow=1))
  coord_cartesian(ylim = c(0,1.1))

cells_placenta_itu <- data.frame(psych::describe(Pheno_Placenta_ITU_reduced_filtered[ ,names_refbased_cells]))
cells_placenta_itu <- cells_placenta_itu[ ,c("mean", "sd")]
rownames(cells_placenta_itu) <- names_refbased_cells

msd_plot_cells_placenta_itu <- ggplot(cells_placenta_itu, aes(x=as.factor(rownames(cells_placenta_itu)), y=mean, fill = as.factor(rownames(cells_cvs)))) +
  geom_bar(position=position_dodge(), stat="identity") +
  scale_fill_manual(values = plColors)+
  theme_minimal()+
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,position=position_dodge(.9))+
  labs(tag = "b", title="Placenta (ITU)")+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.y = element_text(size=10), 
        legend.text = element_text(size=10), 
        legend.title = element_blank(), 
        plot.title=element_text(size=10))+
   scale_y_continuous(limits = c(0,1.1), breaks = c(0, 0.2,0.4, 0.6, 0.8, 1), labels = scales::percent)
  coord_cartesian(ylim = c(0,1.1))

cells_placenta_predo <- data.frame(psych::describe(Pheno_Placenta_PREDO_filtered[ ,names_refbased_cells]))
cells_placenta_predo <- cells_placenta_predo[ ,c("mean", "sd")]
rownames(cells_placenta_predo) <- names_refbased_cells

msd_plot_cells_placenta_predo <- ggplot(cells_placenta_predo, aes(x=as.factor(rownames(cells_placenta_predo)), y=mean, fill = as.factor(rownames(cells_cvs)))) +
  geom_bar(position=position_dodge(), stat="identity") +
  scale_fill_manual(values = plColors)+
  theme_minimal()+
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,position=position_dodge(.9))+
  labs(tag = "c", title="Placenta (PREDO)")+
  scale_y_continuous(limits = c(0,1.1), breaks = c(0, 0.2,0.4, 0.6, 0.8, 1), labels = scales::percent)+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.y = element_text(size=10), 
        legend.text = element_text(size=10), 
        legend.title = element_blank(), 
        plot.title=element_text(size=10))+
  coord_cartesian(ylim = c(0,1.1))

cells_placenta_BET <- data.frame(psych::describe(Pheno_Placenta_BET_filtered[ ,names_refbased_cells]))
cells_placenta_BET <- cells_placenta_BET[ ,c("mean", "sd")]
rownames(cells_placenta_BET) <- names_refbased_cells

msd_plot_cells_placenta_BET <- ggplot(cells_placenta_BET, aes(x=as.factor(rownames(cells_placenta_BET)), y=mean, fill = as.factor(rownames(cells_cvs)))) +
  geom_bar(position=position_dodge(), stat="identity") +
  scale_fill_manual(values = plColors)+
  theme_minimal()+
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,position=position_dodge(.9))+
  labs(tag = "d", title="Placenta (BET)")+
  scale_y_continuous(limits = c(0,1.1), breaks = c(0, 0.2,0.4, 0.6, 0.8, 1), labels = scales::percent)+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.y = element_text(size=10), 
        legend.text = element_text(size=10), 
        legend.title = element_blank(), 
        plot.title=element_text(size=10))+
  coord_cartesian(ylim = c(0,1.1))

msd_plot_cells_cvs
msd_plot_cells_placenta_itu
msd_plot_cells_placenta_predo
msd_plot_cells_placenta_BET
```

``` r
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
```

``` r
gglegend<-g_legend(msd_plot_cells_cvs)
```

``` r
grid.arrange(arrangeGrob(msd_plot_cells_cvs + theme(legend.position="none"),
              CVS_ITU_RPC_Plot_arrange + theme(legend.position="none", axis.text.y = element_blank()),          
              msd_plot_cells_placenta_itu + theme(legend.position="none"),
              Placenta_ITU_RPC_Plot_arrange + theme(legend.position="none", axis.text.y = element_blank()),
              msd_plot_cells_placenta_predo + theme(legend.position="none"),
              Placenta_PREDO_RPC_Plot_arrange + theme(legend.position="none", axis.text.y = element_blank()),
              msd_plot_cells_placenta_BET + theme(legend.position="none"),
              Placenta_BET_RPC_Plot_arrange + theme(legend.position="none", axis.text.y = element_blank()),
              nrow=4),
             gglegend, nrow=2, heights=c(15, 2))
```

``` r
ggsave("Results/RPC_Placenta_combined_cells_estimate.pdf",
grid.arrange(arrangeGrob(msd_plot_cells_cvs + theme(legend.position="none"),
              CVS_ITU_RPC_Plot_arrange + theme(legend.position="none", axis.text.y = element_blank()),          
              msd_plot_cells_placenta_itu + theme(legend.position="none"),
              Placenta_ITU_RPC_Plot_arrange + theme(legend.position="none", axis.text.y = element_blank()),
              msd_plot_cells_placenta_predo + theme(legend.position="none"),
              Placenta_PREDO_RPC_Plot_arrange + theme(legend.position="none", axis.text.y = element_blank()),
              msd_plot_cells_placenta_BET + theme(legend.position="none"),
              Placenta_BET_RPC_Plot_arrange + theme(legend.position="none", axis.text.y = element_blank()),
              nrow=4),
             gglegend, nrow=2, heights=c(15, 2)), width=84, height=150, units="mm", dpi=600, scale=2, device = cairo_pdf)
```

[to the top](#top)

# Compare data sets in RPC cell type proportions

``` r
colwise(mean)(Pheno_CVS_ITU_filtered[c("Trophoblasts")])
colwise(mean)(Pheno_CVS_ITU_filtered[c("Stromal")])
colwise(mean)(Pheno_CVS_ITU_filtered[c("Hofbauer")])
colwise(mean)(Pheno_CVS_ITU_filtered[c("Endothelial")])
colwise(mean)(Pheno_CVS_ITU_filtered[c("nRBC")])
colwise(mean)(Pheno_CVS_ITU_filtered[c("Syncytiotrophoblast")])
```

``` r
colwise(median)(Pheno_CVS_ITU_filtered[c("Trophoblasts")])
colwise(median)(Pheno_CVS_ITU_filtered[c("Stromal")])
colwise(median)(Pheno_CVS_ITU_filtered[c("Hofbauer")])
colwise(median)(Pheno_CVS_ITU_filtered[c("Endothelial")])
colwise(median)(Pheno_CVS_ITU_filtered[c("nRBC")])
colwise(median)(Pheno_CVS_ITU_filtered[c("Syncytiotrophoblast")])
```

## Placenta: CVS vs. term Placenta in ITU

``` r
Pheno_Placenta_CVS_ITU_filtered <- merge(Pheno_CVS_ITU_filtered, Pheno_Placenta_ITU_reduced_filtered, by="Sample_Name") # n = 85
n_Placenta_CVS_filtered <- nrow(Pheno_Placenta_CVS_ITU_filtered)
```

``` r
colwise(mean)(Pheno_Placenta_CVS_ITU_filtered[c("Trophoblasts.x", "Trophoblasts.y")])
colwise(mean)(Pheno_Placenta_CVS_ITU_filtered[c("Stromal.x", "Stromal.y")])
colwise(mean)(Pheno_Placenta_CVS_ITU_filtered[c("Hofbauer.x", "Hofbauer.y")])
colwise(mean)(Pheno_Placenta_CVS_ITU_filtered[c("Endothelial.x", "Endothelial.y")])
colwise(mean)(Pheno_Placenta_CVS_ITU_filtered[c("nRBC.x", "nRBC.y")])
colwise(mean)(Pheno_Placenta_CVS_ITU_filtered[c("Syncytiotrophoblast.x", "Syncytiotrophoblast.y")])
```

``` r
colwise(median)(Pheno_Placenta_CVS_ITU_filtered[c("Trophoblasts.x", "Trophoblasts.y")])
colwise(median)(Pheno_Placenta_CVS_ITU_filtered[c("Stromal.x", "Stromal.y")])
colwise(median)(Pheno_Placenta_CVS_ITU_filtered[c("Hofbauer.x", "Hofbauer.y")])
colwise(median)(Pheno_Placenta_CVS_ITU_filtered[c("Endothelial.x", "Endothelial.y")])
colwise(median)(Pheno_Placenta_CVS_ITU_filtered[c("nRBC.x", "nRBC.y")])
colwise(median)(Pheno_Placenta_CVS_ITU_filtered[c("Syncytiotrophoblast.x", "Syncytiotrophoblast.y")])
```

use Wilcoxon (paired)

``` r
# trophoblasts
wilcox.test(Pheno_Placenta_CVS_ITU_filtered$Trophoblasts.x, Pheno_Placenta_CVS_ITU_filtered$Trophoblasts.y, paired=T)
t.test(Pheno_Placenta_CVS_ITU_filtered$Trophoblasts.x, Pheno_Placenta_CVS_ITU_filtered$Trophoblasts.y, paired=T) # with t
pi_T <- wilcox.test(Pheno_Placenta_CVS_ITU_filtered$Trophoblasts.x, Pheno_Placenta_CVS_ITU_filtered$Trophoblasts.y, paired=T)$p.value # p from Wilxocon
wilcoxonZ(Pheno_Placenta_CVS_ITU_filtered$Trophoblasts.x, Pheno_Placenta_CVS_ITU_filtered$Trophoblasts.y, paired=T) # compute z
wilcoxonZ(Pheno_Placenta_CVS_ITU_filtered$Trophoblasts.x, Pheno_Placenta_CVS_ITU_filtered$Trophoblasts.y, paired=T) / sqrt(n_Placenta_CVS_filtered) # r

# stromal
wilcox.test(Pheno_Placenta_CVS_ITU_filtered$Stromal.x, Pheno_Placenta_CVS_ITU_filtered$Stromal.y, paired=T)
t.test(Pheno_Placenta_CVS_ITU_filtered$Stromal.x, Pheno_Placenta_CVS_ITU_filtered$Stromal.y, paired=T)
pi_ST <- wilcox.test(Pheno_Placenta_CVS_ITU_filtered$Stromal.x,Pheno_Placenta_CVS_ITU_filtered$Stromal.y, paired=T)$p.value
wilcoxonZ(Pheno_Placenta_CVS_ITU_filtered$Stromal.x, Pheno_Placenta_CVS_ITU_filtered$Stromal.y, paired=T) # compute z
wilcoxonZ(Pheno_Placenta_CVS_ITU_filtered$Stromal.x, Pheno_Placenta_CVS_ITU_filtered$Stromal.y, paired=T) / sqrt(n_Placenta_CVS_filtered) # r

# Hofbauer
wilcox.test(Pheno_Placenta_CVS_ITU_filtered$Hofbauer.x, Pheno_Placenta_CVS_ITU_filtered$Hofbauer.y, paired=T)
t.test(Pheno_Placenta_CVS_ITU_filtered$Hofbauer.x, Pheno_Placenta_CVS_ITU_filtered$Hofbauer.y, paired=T)
pi_H <- wilcox.test(Pheno_Placenta_CVS_ITU_filtered$Hofbauer.x, Pheno_Placenta_CVS_ITU_filtered$Hofbauer.y, paired=T)$p.value
wilcoxonZ(Pheno_Placenta_CVS_ITU_filtered$Hofbauer.x, Pheno_Placenta_CVS_ITU_filtered$Hofbauer.y, paired=T) # compute z
wilcoxonZ(Pheno_Placenta_CVS_ITU_filtered$Hofbauer.x, Pheno_Placenta_CVS_ITU_filtered$Hofbauer.y, paired=T) / sqrt(n_Placenta_CVS_filtered) # r

# endothelial
wilcox.test(Pheno_Placenta_CVS_ITU_filtered$Endothelial.x, Pheno_Placenta_CVS_ITU_filtered$Endothelial.y, paired=T)
t.test(Pheno_Placenta_CVS_ITU_filtered$Endothelial.x, Pheno_Placenta_CVS_ITU_filtered$Endothelial.y, paired=T)
pi_E<- wilcox.test(Pheno_Placenta_CVS_ITU_filtered$Endothelial.x,Pheno_Placenta_CVS_ITU_filtered$Endothelial.y, paired=T)$p.value
wilcoxonZ(Pheno_Placenta_CVS_ITU_filtered$Endothelial.x, Pheno_Placenta_CVS_ITU_filtered$Endothelial.y, paired=T) # compute z
wilcoxonZ(Pheno_Placenta_CVS_ITU_filtered$Endothelial.x, Pheno_Placenta_CVS_ITU_filtered$Endothelial.y, paired=T) / sqrt(n_Placenta_CVS_filtered) # r

# nRBC
wilcox.test(Pheno_Placenta_CVS_ITU_filtered$nRBC.x, Pheno_Placenta_CVS_ITU_filtered$nRBC.y, paired=T)
t.test(Pheno_Placenta_CVS_ITU_filtered$nRBC.x, Pheno_Placenta_CVS_ITU_filtered$nRBC.y, paired=T)
pi_n<- wilcox.test(Pheno_Placenta_CVS_ITU_filtered$nRBC.x, Pheno_Placenta_CVS_ITU_filtered$nRBC.y, paired=T)$p.value
wilcoxonZ(Pheno_Placenta_CVS_ITU_filtered$nRBC.x, Pheno_Placenta_CVS_ITU_filtered$nRBC.y, paired=T) # compute z
wilcoxonZ(Pheno_Placenta_CVS_ITU_filtered$nRBC.x, Pheno_Placenta_CVS_ITU_filtered$nRBC.y, paired=T) / sqrt(n_Placenta_CVS_filtered) # r

# syncytiotrophoblast
wilcox.test(Pheno_Placenta_CVS_ITU_filtered$Syncytiotrophoblast.x, Pheno_Placenta_CVS_ITU_filtered$Syncytiotrophoblast.y, paired=T)
t.test(Pheno_Placenta_CVS_ITU_filtered$Syncytiotrophoblast.x, Pheno_Placenta_CVS_ITU_filtered$Syncytiotrophoblast.y, paired=T)
pi_S<- wilcox.test(Pheno_Placenta_CVS_ITU_filtered$Syncytiotrophoblast.x, Pheno_Placenta_CVS_ITU_filtered$Syncytiotrophoblast.y, paired=T)$p.value
wilcoxonZ(Pheno_Placenta_CVS_ITU_filtered$Syncytiotrophoblast.x, Pheno_Placenta_CVS_ITU_filtered$Syncytiotrophoblast.y, paired=T) # compute z
wilcoxonZ(Pheno_Placenta_CVS_ITU_filtered$Syncytiotrophoblast.x, Pheno_Placenta_CVS_ITU_filtered$Syncytiotrophoblast.y, paired=T) / sqrt(n_Placenta_CVS_filtered) # r

# 
pci_placenta <- c(pi_T, pi_ST, pi_H, pi_E, pi_n, pi_S) # all p values
p.adjust(pci_placenta, method = "bonferroni", n = length(pci_placenta)) # multiple test correction with Bonferroni
```

[to the top](#top)

## Term Placentas

make one data frame

``` r
RPC_Cells_ITU <- Pheno_Placenta_ITU_reduced_filtered[ ,names_refbased_cells]
RPC_Cells_ITU$group <- rep("ITU", nrow(RPC_Cells_ITU))

RPC_Cells_PREDO <- Pheno_Placenta_PREDO_filtered[ ,names_refbased_cells]
RPC_Cells_PREDO$group <- rep("PREDO", nrow(RPC_Cells_PREDO))

RPC_Cells_BET <- Pheno_Placenta_BET_filtered[ ,names_refbased_cells]
RPC_Cells_BET$group <- rep("BET", nrow(RPC_Cells_BET))

RPC_Cells_term <- rbind(RPC_Cells_ITU, RPC_Cells_PREDO, RPC_Cells_BET)
RPC_Cells_term$group <- as.factor(RPC_Cells_term$group)
```

``` r
levels(RPC_Cells_term$group)
RPC_Cells_term$group <- factor(RPC_Cells_term$group, levels = c("ITU", "PREDO", "BET"))
levels(RPC_Cells_term$group)
```

``` r
save(RPC_Cells_term, file="Input_Data_prepared/RPC_Cells_term.Rdata")
```

*nonparametric inference for multivariate data* We use the npmv package
Unlike in classical multivariate analysis of variance, multivariate
normality is not required for the data. In fact, the different response
variables may even be measured on different scales (binary, ordinal,
quantitative). p values are calculated for overall tests (permutation
tests and F approximations), and, using multiple testing algorithms
which control the familywise error rate, significant subsets of response
variables and factor levels are identified. The package may be used for
low- or high- dimensional data with small or with large sample sizes and
many or few factor levels.

Typical global statistical hypotheses in this context are the following.
“Are the a samples from the same population (multivariate
distribution)?”

``` r
nonpartest(Trophoblasts|Stromal|Syncytiotrophoblast|Endothelial|nRBC|Hofbauer ~ group,RPC_Cells_term)
```

In addition to the F -distribution approximations that are provided,
each of the four test statistics is also used as the basis for a
multivariate permutation or randomization test. To this end, the N data
vectors are permuted, and the multivariate test statistics recalculated
each time. For each of the four tests, these resulting values form the
respective distribution whose quantiles are used to determine the p
value of the corresponding permutation test (if all N! permutations are
performed) or randomization test (if a predetermined number of random
permutations is performed). The relative effects quantify the tendencies
observed in the data in term of probabilities. For example, the samples
from BET tend to have larger proportions of trophoblasts compared to the
other cohorts. The probability that a randomly chosen sample from BET
exhibits a larger proportion of trophoblasts than a randomly chosen
sample from the other group is 0.856.

follow-up test:

``` r
ssnonpartest(Trophoblasts|Stromal|Syncytiotrophoblast|Endothelial|nRBC|Hofbauer ~group,RPC_Cells_term, alpha = 0.01, factors.and.variables = TRUE)
```

Regarding the group (factor levels), all pairwise comparisons are
significant. Regarding the variables, every one of the six variables
alone shows a significant difference between the groups and so does the
combination of variables

``` r
colwise(median)(Pheno_Placenta_ITU_reduced_filtered[c("Trophoblasts")])
colwise(median)(Pheno_Placenta_ITU_reduced_filtered[c("Stromal")])
colwise(median)(Pheno_Placenta_ITU_reduced_filtered[c("Hofbauer")])
colwise(median)(Pheno_Placenta_ITU_reduced_filtered[c("Endothelial")])
colwise(median)(Pheno_Placenta_ITU_reduced_filtered[c("nRBC")])
colwise(median)(Pheno_Placenta_ITU_reduced_filtered[c("Syncytiotrophoblast")])
```

``` r
colwise(mean)(Pheno_Placenta_ITU_reduced_filtered[c("Trophoblasts")])
colwise(mean)(Pheno_Placenta_ITU_reduced_filtered[c("Stromal")])
colwise(mean)(Pheno_Placenta_ITU_reduced_filtered[c("Hofbauer")])
colwise(mean)(Pheno_Placenta_ITU_reduced_filtered[c("Endothelial")])
colwise(mean)(Pheno_Placenta_ITU_reduced_filtered[c("nRBC")])
colwise(mean)(Pheno_Placenta_ITU_reduced_filtered[c("Syncytiotrophoblast")])
```

``` r
colwise(median)(Pheno_Placenta_PREDO_filtered[c("Trophoblasts")])
colwise(median)(Pheno_Placenta_PREDO_filtered[c("Stromal")])
colwise(median)(Pheno_Placenta_PREDO_filtered[c("Hofbauer")])
colwise(median)(Pheno_Placenta_PREDO_filtered[c("Endothelial")])
colwise(median)(Pheno_Placenta_PREDO_filtered[c("nRBC")])
colwise(median)(Pheno_Placenta_PREDO_filtered[c("Syncytiotrophoblast")])
```

``` r
colwise(mean)(Pheno_Placenta_PREDO_filtered[c("Trophoblasts")])
colwise(mean)(Pheno_Placenta_PREDO_filtered[c("Stromal")])
colwise(mean)(Pheno_Placenta_PREDO_filtered[c("Hofbauer")])
colwise(mean)(Pheno_Placenta_PREDO_filtered[c("Endothelial")])
colwise(mean)(Pheno_Placenta_PREDO_filtered[c("nRBC")])
colwise(mean)(Pheno_Placenta_PREDO_filtered[c("Syncytiotrophoblast")])
```

``` r
colwise(median)(Pheno_Placenta_BET_filtered[c("Trophoblasts")])
colwise(median)(Pheno_Placenta_BET_filtered[c("Stromal")])
colwise(median)(Pheno_Placenta_BET_filtered[c("Hofbauer")])
colwise(median)(Pheno_Placenta_BET_filtered[c("Endothelial")])
colwise(median)(Pheno_Placenta_BET_filtered[c("nRBC")])
colwise(median)(Pheno_Placenta_BET_filtered[c("Syncytiotrophoblast")])
```

``` r
colwise(mean)(Pheno_Placenta_BET_filtered[c("Trophoblasts")])
colwise(mean)(Pheno_Placenta_BET_filtered[c("Stromal")])
colwise(mean)(Pheno_Placenta_BET_filtered[c("Hofbauer")])
colwise(mean)(Pheno_Placenta_BET_filtered[c("Endothelial")])
colwise(mean)(Pheno_Placenta_BET_filtered[c("nRBC")])
colwise(mean)(Pheno_Placenta_BET_filtered[c("Syncytiotrophoblast")])
```

[to the top](#top)

# Phenotype relationships

## CVS (ITU)

``` r
png(file="Results/Pheno_cor_CVS_cells.png", width= 12000, height=6000, res=600)
ggpairs(Pheno_CVS_ITU_filtered[,c(9:19, 42, 43,45,46, 66, 108)])
dev.off()
```

correlation with gestational age

``` r
# correlation cell types with gestational age at sampling
r_gestage_cvs_cell_cors <- sapply(Pheno_CVS_ITU_filtered[,names_refbased_cells], FUN=function(x, y) cor.test(x, y, method="spearman", exact=F)$estimate, y=Pheno_CVS_ITU_filtered$gestage_at_CVS_weeks)
p_gestage_cvs_cell_cors <- sapply(Pheno_CVS_ITU_filtered[,names_refbased_cells], FUN=function(x, y) cor.test(x, y, method="spearman", exact=F)$p.value, y=Pheno_CVS_ITU_filtered$gestage_at_CVS_weeks)

r_gestage_cvs_cell_cors
p_gestage_cvs_cell_cors
```

``` r
# correct for multiple tests
p.adjust(p_gestage_cvs_cell_cors,  method="bonferroni", n=6)
```

GA: Trohpoblasts: r = -0.315 p \< 0.001 Synctiotrophoblast: r = 0.362 p
\< 0.001

``` r
cor_GA_trophoblasts_cvs <- ggscatter(Pheno_CVS_ITU_filtered, x = "gestage_at_CVS_weeks", y = "Trophoblasts", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "gestational age (weeks)\n at CVS sampling", ylab = "proportion \n Trophoblast cells")+
          stat_cor(method = "spearman", label.x = 10, label.y = 0.5,r.digits = 2, aes(label = paste("'r = '", ..r..,sep = "~", "'**'")))+
          scale_x_continuous(breaks=c(10,12,14,16))+
          scale_y_continuous(breaks=c(0.2,0.3,0.4,0.5))+
  labs(tag="a", size=10)+
     theme(axis.text.x = element_text(size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

cor_GA_syncytiotrophoblasts_cvs <- ggscatter(Pheno_CVS_ITU_filtered, x = "gestage_at_CVS_weeks", y = "Syncytiotrophoblast",
          add = "reg.line", conf.int = TRUE, 
          xlab = "gestational age (weeks)\n at CVS sampling", ylab = "proportion \n Syncytiotrophoblast cells")+
          stat_cor(method = "spearman", label.x = 10, label.y = 0.7,r.digits = 2, aes(label = paste("'r = '", ..r..,sep = "~", "'**'")))+
          scale_x_continuous(breaks=c(10,12,14,16))+
          scale_y_continuous(breaks=c(0.5,0.6,0.7))+
    theme(axis.text.x = element_text(size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

GA_CVS <-ggarrange(cor_GA_trophoblasts_cvs, cor_GA_syncytiotrophoblasts_cvs,align="hv", nrow=1)
GA_CVS_a <- annotate_figure(GA_CVS, top = text_grob("CVS (ITU)", size = 10))

GA_CVS_a
```

``` r
png(file="Results/GA_CVS_trophoblasts_cor.png", width= 3000, height=2000, res=600)
cor_GA_trophoblasts_cvs
dev.off()

png(file="Results/GA_CVS_synctiotrophoblast_cor.png", width= 3000, height=2000, res=600)
cor_GA_syncytiotrophoblasts_cvs
dev.off()
```

differences between sexes

``` r
# Trophoblast
wilcox.test(Pheno_CVS_ITU_filtered$Trophoblasts ~ Pheno_CVS_ITU_filtered$Child_Sex)
t.test(Pheno_CVS_ITU_filtered$Trophoblasts ~ Pheno_CVS_ITU_filtered$Child_Sex) # with t
pi_T_sex <- wilcox.test(Pheno_CVS_ITU_filtered$Trophoblasts ~ Pheno_CVS_ITU_filtered$Child_Sex)$p.value # p from Wilxocon

# Hofbauer
wilcox.test(Pheno_CVS_ITU_filtered$Hofbauer ~ Pheno_CVS_ITU_filtered$Child_Sex)
t.test(Pheno_CVS_ITU_filtered$Hofbauer ~ Pheno_CVS_ITU_filtered$Child_Sex) # with t
pi_H_sex <- wilcox.test(Pheno_CVS_ITU_filtered$Hofbauer ~ Pheno_CVS_ITU_filtered$Child_Sex)$p.value # p from Wilxocon

# endothelial
wilcox.test(Pheno_CVS_ITU_filtered$Endothelial ~ Pheno_CVS_ITU_filtered$Child_Sex)
t.test(Pheno_CVS_ITU_filtered$Endothelial ~ Pheno_CVS_ITU_filtered$Child_Sex) # with t
pi_E_sex <- wilcox.test(Pheno_CVS_ITU_filtered$Endothelial ~ Pheno_CVS_ITU_filtered$Child_Sex)$p.value # p from Wilxocon

# stromal
wilcox.test(Pheno_CVS_ITU_filtered$Stromal ~ Pheno_CVS_ITU_filtered$Child_Sex)
t.test(Pheno_CVS_ITU_filtered$Stromal ~ Pheno_CVS_ITU_filtered$Child_Sex) # with t
pi_St_sex <- wilcox.test(Pheno_CVS_ITU_filtered$Stromal ~ Pheno_CVS_ITU_filtered$Child_Sex)$p.value # p from Wilxocon

# syncytiotrophoblast
wilcox.test(Pheno_CVS_ITU_filtered$Syncytiotrophoblast ~ Pheno_CVS_ITU_filtered$Child_Sex)
t.test(Pheno_CVS_ITU_filtered$Syncytiotrophoblast ~ Pheno_CVS_ITU_filtered$Child_Sex) # with t
pi_S_sex <- wilcox.test(Pheno_CVS_ITU_filtered$Syncytiotrophoblast ~ Pheno_CVS_ITU_filtered$Child_Sex)$p.value # p from Wilxocon

# nRBC
wilcox.test(Pheno_CVS_ITU_filtered$nRBC ~ Pheno_CVS_ITU_filtered$Child_Sex)
t.test(Pheno_CVS_ITU_filtered$nRBC ~ Pheno_CVS_ITU_filtered$Child_Sex) # with t
pi_nR_sex <- wilcox.test(Pheno_CVS_ITU_filtered$nRBC ~ Pheno_CVS_ITU_filtered$Child_Sex)$p.value # p from Wilxocon
```

``` r
# correct for multiple tests
p.adjust(c(pi_T_sex, pi_S_sex, pi_H_sex, pi_E_sex, pi_nR_sex, pi_St_sex), method="bonferroni", n=6)
```

no significant differences between sexes

[to the top](#top)

## Placenta (ITU)

``` r
png(file="Results/Pheno_cor_Placenta_ITU_cells.png", width= 12000, height=6000, res=600)
ggpairs(Pheno_Placenta_ITU_reduced_filtered[,c(9:22,45, 46,48,49, 111)])
dev.off()
```

gestational age

``` r
# correlation cell types with gestational age at sampling
r_gestage_feplacenta_cell_cors <- sapply(Pheno_Placenta_ITU_reduced_filtered[,names_refbased_cells], FUN=function(x, y) cor.test(x, y, method="spearman", exact=F)$estimate, y=Pheno_Placenta_ITU_reduced_filtered$Gestational_Age_Weeks)
p_gestage_feplacenta_cell_cors <- sapply(Pheno_Placenta_ITU_reduced_filtered[,names_refbased_cells], FUN=function(x, y) cor.test(x, y, method="spearman", exact=F)$p.value, y=Pheno_Placenta_ITU_reduced_filtered$Gestational_Age_Weeks)

r_gestage_feplacenta_cell_cors
p_gestage_feplacenta_cell_cors
```

``` r
# correct for multiple tests
p.adjust(p_gestage_feplacenta_cell_cors, method="bonferroni", n=6)
```

GA: nothing significant

``` r
cor_GA_trophoblasts_itu <- ggscatter(Pheno_Placenta_ITU_reduced_filtered, x = "Gestational_Age_Weeks", y = "Trophoblasts", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "gestational age (weeks)", ylab = "proportion \n Trophoblast cells")+
          stat_cor(method = "spearman", label.x = 29, label.y = 0.3, r.digits = 2, aes(label = paste("'r = '", ..r..,sep = "~", "''")))+
  labs(tag="b")+
          scale_x_continuous(breaks=c(28,30,32,34,36,38,40,42))+
          scale_y_continuous(breaks=c(0.1,0.2,0.3))+
            theme(axis.text.x = element_text(size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

cor_GA_syncytiotrophoblasts_itu <- ggscatter(Pheno_Placenta_ITU_reduced_filtered, x = "Gestational_Age_Weeks", y = "Syncytiotrophoblast", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "gestational age (weeks)", ylab = "proportion \n Syncytiotrophoblast cells")+
          stat_cor(method = "spearman", label.x = 29, label.y = 1, r.digits = 1, aes(label = paste("'r = '", ..r..,sep = "~", "''")))+
          scale_x_continuous(breaks=c(28,30,32,34,36,38,40,42))+
          scale_y_continuous(breaks=c(0.6,0.7,0.8,0.9,1.0))+
        theme(axis.text.x = element_text(size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

GA_ITU <-ggarrange(cor_GA_trophoblasts_itu, cor_GA_syncytiotrophoblasts_itu,align="hv", nrow=1)
GA_ITU_a <- annotate_figure(GA_ITU, top = text_grob("Placenta (ITU)", size = 10))

GA_ITU_a
```

differences between sexes

``` r
# Trophoblast
wilcox.test(Pheno_Placenta_ITU_reduced_filtered$Trophoblasts ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex)
t.test(Pheno_Placenta_ITU_reduced_filtered$Trophoblasts ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex) # with t
p_T_sex <- wilcox.test(Pheno_Placenta_ITU_reduced_filtered$Trophoblasts ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex)$p.value # p from Wilxocon

# Hofbauer
wilcox.test(Pheno_Placenta_ITU_reduced_filtered$Hofbauer ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex)
t.test(Pheno_Placenta_ITU_reduced_filtered$Hofbauer ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex) # with t
p_H_sex <- wilcox.test(Pheno_Placenta_ITU_reduced_filtered$Hofbauer ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex)$p.value # p from Wilxocon

# endothelial
wilcox.test(Pheno_Placenta_ITU_reduced_filtered$Endothelial ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex)
t.test(Pheno_Placenta_ITU_reduced_filtered$Endothelial ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex) # with t
p_E_sex <- wilcox.test(Pheno_Placenta_ITU_reduced_filtered$Endothelial ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex)$p.value # p from Wilxocon

# syncytiotrophoblast
wilcox.test(Pheno_Placenta_ITU_reduced_filtered$Syncytiotrophoblast ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex)
t.test(Pheno_Placenta_ITU_reduced_filtered$Syncytiotrophoblast ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex) # with t
p_S_sex <- wilcox.test(Pheno_Placenta_ITU_reduced_filtered$Syncytiotrophoblast ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex)$p.value # p from Wilxocon

# stromal
wilcox.test(Pheno_Placenta_ITU_reduced_filtered$Stromal ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex)
t.test(Pheno_Placenta_ITU_reduced_filtered$Stromal ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex) # with t
p_St_sex <- wilcox.test(Pheno_Placenta_ITU_reduced_filtered$Stromal ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex)$p.value # p from Wilxocon

# nRBC
wilcox.test(Pheno_Placenta_ITU_reduced_filtered$nRBC ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex)
t.test(Pheno_Placenta_ITU_reduced_filtered$nRBC ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex) # with t
p_nR_sex <- wilcox.test(Pheno_Placenta_ITU_reduced_filtered$nRBC ~ Pheno_Placenta_ITU_reduced_filtered$Child_Sex)$p.value # p from Wilxocon
```

no significant differences

[to the top](#top)

## Placenta (PREDO)

``` r
png(file="Results/Pheno_cor_Placenta_PREDO_cells.png", width= 12000, height=6000, res=600)
ggpairs(Pheno_Placenta_PREDO_filtered[,c(13:20, 22,24,25,27, 102)])
dev.off()
```

``` r
# correlation cell types with gestational age at sampling
r_gestage_deplacenta_cell_cors <- sapply(Pheno_Placenta_PREDO_filtered[,names_refbased_cells], FUN=function(x, y) cor.test(x, y, method="spearman", exact=F)$estimate, y=Pheno_Placenta_PREDO_filtered$Gestational_Age)
p_gestage_deplacenta_cell_cors <- sapply(Pheno_Placenta_PREDO_filtered[,names_refbased_cells], FUN=function(x, y) cor.test(x, y, method="spearman", exact=F)$p.value, y=Pheno_Placenta_PREDO_filtered$Gestational_Age)

r_gestage_deplacenta_cell_cors
p_gestage_deplacenta_cell_cors
```

``` r
# correct for multiple tests
p.adjust(p_gestage_deplacenta_cell_cors, method="bonferroni", n=6)
```

Gestational Age: trophoblasts r = -0.254, p = 0.016 syncytiotrophoblast
r = 0.246, p = 0.021

``` r
cor_GA_trophoblasts_predo <- ggscatter(Pheno_Placenta_PREDO_filtered, x = "Gestational_Age", y = "Trophoblasts", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "gestational age (weeks)", ylab = "proportion \n Trophoblast cells")+
          stat_cor(method = "spearman", label.x = 32, label.y = 0.3,r.digits = 2, aes(label = paste("'r = '", ..r..,sep = "~", "''")))+
  labs(tag="c", size=10)+
          scale_x_continuous(breaks=c(32,34,36,38,40,42))+
          scale_y_continuous(breaks=c(0.0,0.1,0.2,0.3, 0.3))+
  theme(axis.text.x = element_text(size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())
  
cor_GA_syncytiotrophoblasts_predo <- ggscatter(Pheno_Placenta_PREDO_filtered, x = "Gestational_Age", y = "Syncytiotrophoblast", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "gestational age (weeks)", ylab = "proportion \n Syncytiotrophoblast cells")+
          stat_cor(method = "spearman", label.x = 32, label.y = 1.0,r.digits = 2, aes(label = paste("'r = '", ..r..,sep = "~", "''")))+
          scale_x_continuous(breaks=c(32,34,36,38,40,42))+
          scale_y_continuous(breaks=c(0.6,0.7,0.8,0.9,1.0))+
  theme(axis.text.x = element_text(size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

GA_PREDO <-ggarrange(cor_GA_trophoblasts_predo, cor_GA_syncytiotrophoblasts_predo,align="hv", nrow=1)
GA_PREDO_a <- annotate_figure(GA_PREDO, top = text_grob("Placenta (PREDO)", size = 10))

GA_PREDO_a
```

``` r
png(file="Results/GA_PlacentaPredo_Trophoblasts_cor.png", width= 3000, height=2000, res=600)
cor_GA_trophoblasts_predo
dev.off()

png(file="Results/GA_PlacentaPredo_Synctiotrophoblast_cor.png", width= 3000, height=2000, res=600)
cor_GA_syncytiotrophoblasts_predo
dev.off()
```

differences between sexes

``` r
# Trophoblast
wilcox.test(Pheno_Placenta_PREDO_filtered$Trophoblasts ~ Pheno_Placenta_PREDO_filtered$Child_Sex)
t.test(Pheno_Placenta_PREDO_filtered$Trophoblasts ~ Pheno_Placenta_PREDO_filtered$Child_Sex) # with t
pp_T_sex <- wilcox.test(Pheno_Placenta_PREDO_filtered$Trophoblasts ~ Pheno_Placenta_PREDO_filtered$Child_Sex)$p.value # p from Wilxocon

# Hofbauer
wilcox.test(Pheno_Placenta_PREDO_filtered$Hofbauer ~ Pheno_Placenta_PREDO_filtered$Child_Sex)
t.test(Pheno_Placenta_PREDO_filtered$Hofbauer ~ Pheno_Placenta_PREDO_filtered$Child_Sex) # with t
pp_H_sex <- wilcox.test(Pheno_Placenta_PREDO_filtered$Hofbauer ~ Pheno_Placenta_PREDO_filtered$Child_Sex)$p.value # p from Wilxocon

# endothelial
wilcox.test(Pheno_Placenta_PREDO_filtered$Endothelial ~ Pheno_Placenta_PREDO_filtered$Child_Sex)
t.test(Pheno_Placenta_PREDO_filtered$Endothelial ~ Pheno_Placenta_PREDO_filtered$Child_Sex) # with t
pp_E_sex <- wilcox.test(Pheno_Placenta_PREDO_filtered$Endothelial ~ Pheno_Placenta_PREDO_filtered$Child_Sex)$p.value # p from Wilxocon

# syncytiotrophoblast
wilcox.test(Pheno_Placenta_PREDO_filtered$Syncytiotrophoblast ~ Pheno_Placenta_PREDO_filtered$Child_Sex)
t.test(Pheno_Placenta_PREDO_filtered$Syncytiotrophoblast ~ Pheno_Placenta_PREDO_filtered$Child_Sex) # with t
pp_S_sex <- wilcox.test(Pheno_Placenta_PREDO_filtered$Syncytiotrophoblast ~ Pheno_Placenta_PREDO_filtered$Child_Sex)$p.value # p from Wilxocon

# stromal
wilcox.test(Pheno_Placenta_PREDO_filtered$Stromal ~ Pheno_Placenta_PREDO_filtered$Child_Sex)
t.test(Pheno_Placenta_PREDO_filtered$Stromal ~ Pheno_Placenta_PREDO_filtered$Child_Sex) # with t
pp_St_sex <- wilcox.test(Pheno_Placenta_PREDO_filtered$Stromal ~ Pheno_Placenta_PREDO_filtered$Child_Sex)$p.value # p from Wilxocon

# nRBC
wilcox.test(Pheno_Placenta_PREDO_filtered$nRBC ~ Pheno_Placenta_PREDO_filtered$Child_Sex)
t.test(Pheno_Placenta_PREDO_filtered$nRBC ~ Pheno_Placenta_PREDO_filtered$Child_Sex) # with t
pp_nR_sex <- wilcox.test(Pheno_Placenta_PREDO_filtered$nRBC ~ Pheno_Placenta_PREDO_filtered$Child_Sex)$p.value # p from Wilxocon
```

small difference in trophoblasts (higher in females), but not after
multiple test correction

``` r
# correct for multiple tests
p.adjust(c(pi_S_sex, pi_H_sex, pi_T_sex, pi_E_sex, pi_St_sex, pi_nR_sex), method="bonferroni", n=6)
```

[to the top](#top)

## Placenta (BET)

``` r
png(file="Results/Pheno_cor_Placenta_BET_cells.png", width= 12000, height=6000, res=600)
ggpairs(Pheno_Placenta_BET_filtered[,c(35:44,50,51)])
dev.off()
```

correlation with gestational age

``` r
# correlation cell types with gestational age at sampling
r_gestage_Placenta_BET_cell_cors <- sapply(Pheno_Placenta_BET_filtered[,names_refbased_cells], FUN=function(x, y) cor.test(x, y, method="spearman", exact=F)$estimate, y=Pheno_Placenta_BET_filtered$gestage_weeks)
p_gestage_Placenta_BET_cell_cors <- sapply(Pheno_Placenta_BET_filtered[,names_refbased_cells], FUN=function(x, y) cor.test(x, y, method="spearman", exact=F)$p.value, y=Pheno_Placenta_BET_filtered$gestage_weeks)

r_gestage_Placenta_BET_cell_cors
p_gestage_Placenta_BET_cell_cors
```

``` r
# correct for multiple tests
p.adjust(p_gestage_Placenta_BET_cell_cors, method="bonferroni", n=6)
```

GA: Trohpoblasts: r = -0.415 p \< 0.001 Synctiotrophoblast: r = 0.367 p
\< 0.001

``` r
cor_GA_trophoblasts_BET <- ggscatter(Pheno_Placenta_BET_filtered, x = "gestage_weeks", y = "Trophoblasts", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "gestational age (weeks)", ylab = "proportion \n Trophoblast cells")+
          stat_cor(method = "spearman", label.x = 34, label.y = 0.4, r.digits = 2, aes(label = paste("'r = '", ..r..,sep = "~", "'**'")))+
  labs(tag="d")+
          scale_x_continuous(breaks=c(34,36,38,40,42))+
          scale_y_continuous(breaks=c(0.1,0.2,0.3,0.4))+
            theme(axis.text.x = element_text(size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

cor_GA_syncytiotrophoblasts_BET <- ggscatter(Pheno_Placenta_BET_filtered, x = "gestage_weeks", y = "Syncytiotrophoblast", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "gestational age (weeks)", ylab = "proportion \n Syncytiotrophoblast cells")+
          stat_cor(method = "spearman", label.x = 34, label.y = 0.8, r.digits = 2, aes(label = paste("'r = '", ..r..,sep = "~", "'**'")))+
          scale_x_continuous(breaks=c(34,36,38,40,42))+
          scale_y_continuous(breaks=c(0.5,0.6,0.7,0.8))+
        theme(axis.text.x = element_text(size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=10), axis.title.x=element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

GA_BET <-ggarrange(cor_GA_trophoblasts_BET, cor_GA_syncytiotrophoblasts_BET,align="hv", nrow=1)
GA_BET_a <- annotate_figure(GA_BET, top = text_grob("Placenta (BET)", size = 10))

GA_BET_a
```

``` r
png(file="Results/GA_Placenta_BET_trophoblasts_cor.png", width= 3000, height=2000, res=600)
cor_GA_trophoblasts_BET
dev.off()

png(file="Results/GA_Placenta_BET_synctiotrophoblast_cor.png", width= 3000, height=2000, res=600)
cor_GA_syncytiotrophoblasts_BET
dev.off()
```

differences between sexes

``` r
# Trophoblast
wilcox.test(Pheno_Placenta_BET_filtered$Trophoblasts ~ Pheno_Placenta_BET_filtered$Sex)
t.test(Pheno_Placenta_BET_filtered$Trophoblasts ~ Pheno_Placenta_BET_filtered$Sex) # with t
pie_T_sex <- wilcox.test(Pheno_Placenta_BET_filtered$Trophoblasts ~ Pheno_Placenta_BET_filtered$Sex)$p.value # p from Wilxocon

# Hofbauer
wilcox.test(Pheno_Placenta_BET_filtered$Hofbauer ~ Pheno_Placenta_BET_filtered$Sex)
t.test(Pheno_Placenta_BET_filtered$Hofbauer ~ Pheno_Placenta_BET_filtered$Sex) # with t
pie_H_sex <- wilcox.test(Pheno_Placenta_BET_filtered$Hofbauer ~ Pheno_Placenta_BET_filtered$Sex)$p.value # p from Wilxocon

# endothelial
wilcox.test(Pheno_Placenta_BET_filtered$Endothelial ~ Pheno_Placenta_BET_filtered$Sex)
t.test(Pheno_Placenta_BET_filtered$Endothelial ~ Pheno_Placenta_BET_filtered$Sex) # with t
pie_E_sex <- wilcox.test(Pheno_Placenta_BET_filtered$Endothelial ~ Pheno_Placenta_BET_filtered$Sex)$p.value # p from Wilxocon

# stromal
wilcox.test(Pheno_Placenta_BET_filtered$Stromal ~ Pheno_Placenta_BET_filtered$Sex)
t.test(Pheno_Placenta_BET_filtered$Stromal ~ Pheno_Placenta_BET_filtered$Sex) # with t
pie_St_sex <- wilcox.test(Pheno_Placenta_BET_filtered$Stromal ~ Pheno_Placenta_BET_filtered$Sex)$p.value # p from Wilxocon

# syncytiotrophoblast
wilcox.test(Pheno_Placenta_BET_filtered$Syncytiotrophoblast ~ Pheno_Placenta_BET_filtered$Sex)
t.test(Pheno_Placenta_BET_filtered$Syncytiotrophoblast ~ Pheno_Placenta_BET_filtered$Sex) # with t
pie_S_sex <- wilcox.test(Pheno_Placenta_BET_filtered$Syncytiotrophoblast ~ Pheno_Placenta_BET_filtered$Sex)$p.value # p from Wilxocon

# nRBC
wilcox.test(Pheno_Placenta_BET_filtered$nRBC ~ Pheno_Placenta_BET_filtered$Sex)
t.test(Pheno_Placenta_BET_filtered$nRBC ~ Pheno_Placenta_BET_filtered$Sex) # with t
pie_nR_sex <- wilcox.test(Pheno_Placenta_BET_filtered$nRBC ~ Pheno_Placenta_BET_filtered$Sex)$p.value # p from Wilxocon
```

``` r
# correct for multiple tests
p.adjust(c(pie_T_sex, pie_S_sex, pie_H_sex, pie_E_sex, pie_nR_sex, pie_St_sex), method="bonferroni", n=6)
```

no significant differences between sexes

[to the top](#top)

## arrange plots together

*Fig. 7*

``` r
ggsave("Results/GA_cor.pdf",
grid.arrange(GA_CVS_a,
             GA_ITU_a,
              GA_PREDO_a,          
              GA_BET_a, nrow=4),
              width=84, height=120, units="mm", dpi=600, scale=2)
```

[to the top](#top)
