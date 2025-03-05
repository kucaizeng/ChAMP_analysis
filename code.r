library(data.table)
library(tidyverse)
library(ChAMP)
library(conumee)
library(openxlsx)
library(GenomicRanges)
library(minfi)
library(ComplexHeatmap)
library(unixtools)
set.tempdir("./tmp")


###  ChAMP R

data <- chapm.load(directory="/data/",arraytype="EPIC")
CpG.GUI(CpG=rownames(data$beta),arraytype="EPIC")
champ.QC(
    beta= data$beta,
    pheno= data$pd$Sample_Group,
    mdsPlot= TRUE,
    densityPlot= TRUE,
    PDFplot= TRUE,
    Rplot= TRUE,
    Feature.sel= "None",
    resultsDir="./CHAMP_QCimages")

QC.GUI(
    beta= data$beta,
    pheno= data$pd$Sample_Group,
    arraytype="EPIC")

data_norm <- champ.norm(
    beta= data$beta,
    rgSet= data$rgSet,
    mset= data$ mset,
    resultsDir="./CHAMP_Normalization/",
    method="BMIQ",
    plotBMIQ=FALSE,
    arraytype="EPIC",
    cores=3
)

champ.SVD(
    beta=data_norm,
    rgSet=NULL,
    pd=data$pd,
    RGEffect=FALSE,
    PDFplot=TRUE,
    Rplot=TRUE,
    resultsDir="./CHAMP_SVDimmages"
)

data_Combat <- champ.rumCombat(beta=data_norm,pd=data$pd,batchname=c('Slide'))

data_DMP <- champ.DMP(
    beta=data_norm,
    pheno= data$pd$Sample_Group,
    compare.group=NULL,
    adjPVal=0.05,
    adjust.method="BH",
    arraytype="EPIC"
)

data_DMR <- champ.DMR(
    bata=data_norm,
    pheno= data$pd$Sample_Group,
    compare.group=NULL,
    arraytype="EPIC",
    method="Bumphunter",
    minProbes=7,
    adjPvalDmr=0.05,
    cores=3,
    maxGap=300,
    cutoff=NULL,
    pickCutoff=TRUE,
    smooth=TRUE,
    smoothFunction=loessByCluster,
    useWeights=FALSE,
    permutations=NULL,
    B=250,
    nullMethod="bootstrap",
    meanLassoRadius=375,
    minDmrSeq=1000,
    minDmrSize=50,
    adjPvalProbe=0.05,
    Rplot=T,
    PDFplot=T,
    resultsDir="./CHAMP_ProbeLasso/",
    rmSNPCH=T, 
    fdr=0.05,
    dist=2,
    mafuct=0.05,
    lambda=1000,
    C=2
)

DMR.GUI(DMR=data_DMR)


###   clusterProfiler



#ORA mat
DEG_mat <- read.csv("DEG.csv")
DEGdata <- subset(DEG_mat,DEG_mat$logFC>0&DEG_mat$pvalue<0.05)

DEgenelist <- as.character(DEGdata$gene_name)


#GSEA  sort mat
FCgenelist <- DEGdata$logFC #numeric vector
names(FCgenelist) <- as.character(DEGdata$gene_name) #named vector
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order

library(msigdbr)
Dm_msigdbr <- msigdbr(species="Homo sapiens")
DmGO <- msigdbr(species="Homo sapiens",category="C5") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)
em <- enricher(DEgenelist,TERM2GENE=DmGO[,c(1,3)])
em1 <- GSEA(FCgenelist,TERM2GENE=DmGO[,c(1,3)])


ggo <- groupGO(gene=DEgenelist, OrgDb=org.Hs.eg.db,ont="BP", level=2, readable=F)
egoMF <- enrichGO(DEgenelist, OrgDb=org.Hs.eg.db, ont='MF',pAdjustMethod='BH', pvalueCutoff=0.05, qvalueCutoff=0.2, keyType='SYMBOL')
egoall <- enrichGO(DEgenelist, OrgDb=org.Hs.eg.db, ont='ALL',pAdjustMethod='BH', pvalueCutoff=0.05, qvalueCutoff=0.2, keyType='SYMBOL')

gene.df <- bitr(DEgenelist,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Hs.eg.db)
gene.kegg <- bitr_kegg(gene.df$ENTREZID,fromType="ncbi-geneid",toType="kegg",organism='dme')
ekegg <- enrichKEGG(gene.kegg$kegg, organism='dme',keyType="kegg",pvalueCutoff=0.5,pAdjustMethod='BH',qvalueCutoff=0.5,minGSSize=10,maxGSSize=500,use_internal_data=T)

barplot(egoMF,showCategory=10)
dotplot(egoall,showCategory=10)
dotplot(egoall,title='Top5 GO terms of each sub-class',showCategory=5,split='ONTOLOGY')+ facet_grid(ONTOLOGY~.,scale="free")
emapplot(egoMF,showCategory=10) 
cnetplot(egoMF,showCategory=5)
cnetplot(egoall,showCategory=10,foldChange=FClist,circular=TRUE,colorEdge=TRUE)
goplot(egoMF,showCategory=5)
ridgeplot(egseGO)

library(oncoPredict)
library(gtools)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggplot2)

data<- readRDS('data.RDS')

exp<- assay(data,2) %>% as.data.frame

exp<-fpkmToTpm(exp)

# CTRP2
CTR2_exp<- readRDS('oncopredict/Training Data/CTRP2_Expr (TPM, not log transformed).rds')
CTR2_res<- readRDS('oncopredict/Training Data/CTRP2_Res.rds')
CTRP_tissue<- read.xlsx('oncopredict/GLDS/CTRPv2/Harmonized_CCL_Data_v1.0.xlsx')

CTRP_tissue<- CTRP_tissue %>% filter(cancer.type == 'kidney cancer')
kidney_ctrp<- CTRP_tissue$cvcl.cell.line.name

CTR2_res<- CTR2_res %>% as.data.frame() %>% filter(rownames(CTR2_res) %in% kidney_ctrp) %>% as.matrix()
CTR2_res<- exp(CTR2_res)

calcPhenotype(trainingExprData = CTR2_exp,
              trainingPtype = CTR2_res,
              testExprData = exp,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              cc=TRUE, #biomarker
              rsq = TRUE, #R^2
              pcr = FALSE, #pcr must be FALSE if cc is TRUE
              report_pc = FALSE,  
              removeLowVaringGenesFrom = 'rawData' )

#GDSC2
GDSC2_Expr <- readRDS('oncopredict/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res <- readRDS("oncopredict/Training Data/GDSC2_Res.rds")
GDSC_tissue<- read.xlsx('oncopredict/GLDS/GDSCv2/Cell_Lines_Details.xlsx',2)
GDSC_tissue<- GDSC_tissue %>% filter(Site=='kidney')
GDSC_tissue$COSMIC_ID<- paste('COSMIC',GDSC_tissue$COSMIC_ID,sep='_')
kidney_gdsc<- GDSC_tissue$COSMIC_ID
GDSC2_Res<- GDSC2_Res %>% as.data.frame() %>% filter(rownames(GDSC2_Res) %in% kidney_gdsc) %>% as.matrix
GDSC2_Res<- exp(GDSC2_Res)
exp_log<- log(exp+1)
 calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = exp_log,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              cc=TRUE, #biomarker
              rsq = TRUE, #
              pcr = FALSE, #pcr must be FALSE if cc is TRUE
              report_pc = FALSE,  
              removeLowVaringGenesFrom = 'rawData' )
#
ctrp_drug_prediction<- read.csv('CTRP/DrugPredictions.csv',row.names = 1)
drugRelatedness <- read.xlsx("GLDS/CTRPv2/Harmonized_Compound_Data_v1.0.xlsx")

drugRelatedness<- drugRelatedness %>% filter(Dataset=='GDSC1_v8.2') %>% select(Compound_Name_in_Dataset,Compound_Molecular_Targets)

mut_infor<- read.xlsx('primary_samples_somatic_mutations01_matrix.xlsx',rowNames = TRUE) %>% as.matrix()

glds(drugMat=ctrp_drug_prediction,
     drugRelatedness=drugRelatedness,
     markerMat=mut_infor,
     minMuts=5,
     additionalCovariateMatrix=NULL,
     threshold=0.7)

#plot
glds_p<- read.csv('GDSC_mut/gldsPs.csv',row.names = 1)
glds_p1<- na.omit(glds_p)
glds_p1$gene<- rownames(glds_p1)
glds_p2<- reshape2::melt(glds_p1, id.vars = c("gene"), 
                  variable.name = c('drug'),
                  value.name = 'pvalue')

glds_p3<- glds_p2[order(glds_p2$pvalue),] 
glds_p3<- glds_p3[1:30,] 
glds_p3$name<- paste(glds_p3$drug,glds_p3$gene,sep='_')
glds_p3$type<-'glds'
glds_p3<- glds_p3[,-c(1,2)]

drug_prediction<- read.table('GDSC.V2_all_calcPhenotype_Output/pvalues.txt')
colnames(drug_prediction)<-as.vector(fix)
drug_prediction$gene<- rownames(drug_prediction)
drug_prediction_1<- reshape2::melt(drug_prediction, id.vars = c("gene"), 
                  variable.name = c('drug'),
                  value.name = 'pvalue')

drug_prediction_1$name<-paste(drug_prediction_1$drug,drug_prediction_1$gene,sep='_')
drug_prediction_1$type<-'prediction'
drug_prediction_1<- drug_prediction_1 %>% filter(name %in% glds_p3$name) %>% select (pvalue,name,type)

glds_predict<- rbind(glds_p3,drug_prediction_1)
glds_predict$`-log10(pvalue)`<- -log10(glds_predict$pvalue)

fig01<- ggplot(glds_predict)+ geom_point(aes(x=-log10(pvalue),y=name,color=type))+theme_classic() 
fig01