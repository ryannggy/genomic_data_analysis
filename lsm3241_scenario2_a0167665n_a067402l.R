######Obtain Raw Microarray Data from GEO######
library(GEOquery)
gse <- getGEO('GSE143150',GSEMatrix=FALSE)
gsm <- GSMList(gse)[[1]]

# subset of metadata present in GSM but NOT in GSE
names(Meta(gsm))[!(names(Meta(gsm)) %in% names(Meta(gse)))]
Meta(gsm)[!(names(Meta(gsm)) %in% names(Meta(gse)))]

# characteristics_Ch1
for (gsm in GSMList(gse)) { 
  print(Meta(gsm)[['characteristics_ch1']])
}
genotype <- function(gsm) {
  Meta(gsm)[['characteristics_ch1']][2]
}
treatment <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][3]
}
sapply(GSMList(gse),treatment)
sapply(GSMList(gse),genotype)

# creating data frame of treatment and genotype
pd <- data.frame(treatment=as.factor(sapply(GSMList(gse),treatment)))
pd2 <- data.frame(genotype=as.factor(sapply(GSMList(gse),genotype)))
pd$genotype <- pd2$genotype
pd$genotype <- as.factor(pd$genotype)
pd$treatment <- as.factor(pd$treatment)
levels(pd$treatment) <- c("DM", "FM")
levels(pd$genotype) <- c("JMJD2B_Knockdown", "Wild_Type")

# reading in CEL files with the phenoData
library(oligo)
celfiles <- paste0('genomic_data_analysis/GSE143150/celfiles/',rownames(pd),'.CEL')
data <- read.celfiles(celfiles, phenoData=new("AnnotatedDataFrame", pd))
pData(data)

# rerforming RMA in R
eset <- oligo::rma(data)

library(affy)
par(mfrow=c(1,2))
affy::plotDensity(exprs(data),xlab='log intensity',main="feature level densities before RMA",lwd=2)
affy::plotDensity(exprs(eset),xlab='log intensity',main="feature level densities after RMA",lwd=2)

#2x2 Factorial using LIMMA
library(limma)
(t <- factor(eset$treatment, levels=c("DM","FM")))
(g <- factor(eset$genotype, levels=c("Wild_Type", "JMJD2B_Knockdown")))
(model <- model.matrix( ~ 0 + t:g))
colnames(model) <- c("DM_WT","FM_WT","DM_KD","FM_KD")
contrasts <- makeContrasts(DM_WT-FM_WT,
                           DM_KD-FM_KD,
                           DM_WT-DM_KD,
                           FM_WT-FM_KD,
                           levels=model)
fit <- lmFit(eset, model)
fitted.contrast <- contrasts.fit(fit,contrasts)
(fitted.ebayes <- eBayes(fitted.contrast))

# Annotation of genomic data in AnnotationDb
library(huex10sttranscriptcluster.db)
columns(huex10sttranscriptcluster.db)

# Basic downstream analysis of microarray data: DM_WT-FM_WT
(ps1 <- topTable(fitted.ebayes, coef=1, number=Inf,p.value = 0.05,lfc=1))
#write.csv(ps1, 'ps1.csv')
(interesting_genes1 <- rownames(ps1))
(ig1 <- AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes1,c("SYMBOL","GENENAME"),keytype="PROBEID"))
#write.csv(ig1, 'ig1.csv')
volcanoplot(fitted.ebayes, main=sprintf("DM_WT-FM_WT: %d features pass our cutoffs",nrow(ps1)))
points(ps1[['logFC']],-log10(ps1[['P.Value']]),col='red')

# upregulated vs downregulated genes
(ps1_up <- rownames(ps1[ps1$logFC>0,]))
(ig1_up <- AnnotationDbi::select(huex10sttranscriptcluster.db,ps1_up,c("SYMBOL","GENENAME"),keytype="PROBEID"))
#write.csv(ig1_up, 'ig1_up.csv')
(ps1_down <- rownames(ps1[ps1$logFC<0,]))
(ig1_down <- AnnotationDbi::select(huex10sttranscriptcluster.db,ps1_down,c("SYMBOL","GENENAME"),keytype="PROBEID"))
#write.csv(ig1_down, 'ig1_down.csv')

# heatmap of gene expression
library(RColorBrewer)
eset_of_interest1 <- eset[rownames(ps1),]
(hm1 <- heatmap(exprs(eset_of_interest1),
          labCol=eset$culture,labRow=NA,
          col       = rev(brewer.pal(10, "RdBu")),
          distfun   = function(x) as.dist(1-cor(t(x)))))
jpeg(file="hm1.jpg")


# Basic downstream analysis of microarray data: DM_KD-FM_KD
ps2 <- topTable(fitted.ebayes, coef=2, number=Inf,p.value = 0.05,lfc=1)
#write.csv(ps2, 'ps2.csv')
(interesting_genes2 <- rownames(ps2))
(ig2 <- AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes2,c("SYMBOL","GENENAME"),keytype="PROBEID"))
#write.csv(ig2, 'ig2.csv')
volcanoplot(fitted.ebayes, coef=2, main=sprintf("DM_KD-FM_KD: %d features pass our cutoffs",nrow(ps2)))
points(ps2[['logFC']],-log10(ps2[['P.Value']]),col='red')

# upregulated vs downregulated genes
(ps2_up <- rownames(ps2[ps2$logFC>0,]))
(ig2_up <- AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up,c("SYMBOL","GENENAME"),keytype="PROBEID"))
#write.csv(ig2_up, 'ig2_up.csv')
(ps2_down <- rownames(ps2[ps2$logFC<0,]))
(ig2_down <- AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_down,c("SYMBOL","GENENAME"),keytype="PROBEID"))
#write.csv(ig2_down, 'ig2_down.csv')

# heatmap of gene expression
eset_of_interest2 <- eset[rownames(ps2),]
heatmap(exprs(eset_of_interest2),
        labCol=eset$culture,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

##### unable to obtain results for these two following contrasts, uncomment to try #####
#Basic downstream analysis of microarray data: DM_WT-DM_KD
#ps3 <- topTable(fitted.ebayes, coef=3, number=Inf,p.value = 0.05,lfc=1)
#(interesting_genes3 <- rownames(ps3))
#AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes3,c("SYMBOL","GENENAME"),keytype="PROBEID")
#volcanoplot(fitted.ebayes, coef = 3, main=sprintf("DM_WT-DM_KD: %d features pass our cutoffs",nrow(ps3)))

#Basic downstream analysis of microarray data: FM_WT-FM_KD
#ps4 <- topTable(fitted.ebayes, coef=4, number=Inf,p.value = 0.05,lfc=1)
#(interesting_genes4 <- rownames(ps4))
#AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes4,c("SYMBOL","GENENAME"),keytype="PROBEID")
#volcanoplot(fitted.ebayes, coef = 4, main=sprintf("FM_WT-FM_KD: %d features pass our cutoffs",nrow(ps4)))
#####

# Basic downstream analysis of microarray data: (DM_WT-FM_WT)-(DM_KD-FM_KD)
contrasts2 <- makeContrasts((DM_WT-FM_WT)+(DM_KD-FM_KD),
                            levels=model)
fitted.contrast2 <- contrasts.fit(fit,contrasts2)
fitted.ebayes2 <- eBayes(fitted.contrast2)
topTable(fitted.ebayes2)
fitted.ebayes2

#getting a limited number of probesets
(ps5 <- topTable(fitted.ebayes2,number=Inf,p.value = 0.05,lfc=1, coef=1))
(interesting_genes5 <- rownames(ps5))
(ig5 <-AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes5,c("SYMBOL","GENENAME"),keytype="PROBEID"))
write.csv(ps5, 'ps5.csv')
write.csv(ig5, 'ig5.csv')
volcanoplot(fitted.ebayes2, coef=1, main=sprintf("(DM_WT-FM_WT)+(DM_KD-FM_KD): %d features pass our cutoffs",nrow(ps5)))
points(ps5[['logFC']],-log10(ps5[['P.Value']]),col='red')

(ps5_up <- rownames(ps5[ps5$logFC>0,]))
(ig5_up <- AnnotationDbi::select(huex10sttranscriptcluster.db,ps5_up,c("SYMBOL","GENENAME"),keytype="PROBEID"))
#write.csv(ig5_up, 'ig5_up.csv')
(ps5_down <- rownames(ps5[ps5$logFC<0,]))
(ig5_down <- AnnotationDbi::select(huex10sttranscriptcluster.db,ps5_down,c("SYMBOL","GENENAME"),keytype="PROBEID"))
write.csv(ig5_down, 'ig5_down.csv')

# heatmap of gene expression
eset_of_interest2 <- eset[rownames(ps5),]
heatmap(exprs(eset_of_interest2),
        labCol=eset$culture,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))
