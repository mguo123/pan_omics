library(tidyverse)
library(motifbreakR)
library(MotifDb)
library(rtracklayer)
ibrary(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

variants = snps.from.rsid(rsid='rs133015',dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37,search.genome = BSgenome.Hsapiens.UCSC.hg19)

MotifDb

# data(motifbreakR_motif)
# motifbreakR_motif
data( hocomoco)
 hocomoco
 
 motifs <- query(MotifDb, andStrings=c( "hsapiens"),
                orStrings=c("jaspar2018", "hocomocov11-core-A", "hocomocov11-core-B", "hocomocov11-core-C","HOCOMOCOv11-secondary-A","HOCOMOCOv11-secondary-B","HOCOMOCOv11-secondary-C","HOCOMOCOv11-secondary-D"))
length(motifs)

results = motifbreakR(snpList=variants,filterp=TRUE,pwmList =  motifs , threshold=1e-4,method='ic',bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),BPPARAM = BiocParallel::bpparam())


write.table(data.frame(results), file='/Users/mguo123/Google Drive/1_khavari/omics_project-LD/pan_omics/data/processed/snp_motifs/rs133015_motifbreakr.txt', sep="'\t", col.names=TRUE, row.names=FALSE, quote=FALSE )


pdf('/Users/mguo123/Google Drive/1_khavari/omics_project-LD/pan_omics/data/processed/snp_motifs/rs133015_motifplot.pdf')
plotMB(results = results, rsid = 'rs133015',effect=c('strong','weak'))
dev.off()
